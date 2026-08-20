#!/usr/bin/env python3
"""Execute locked Phase 5B continuous basal-classical axis analyses."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("05_results/.matplotlib").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

try:
    import diptest as py_diptest
except Exception:  # pragma: no cover - availability is environment-specific.
    py_diptest = None


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results" / "tables"
FIG_DIR = ROOT / "05_results" / "figures"
REPORT_DIR = ROOT / "04_analysis" / "07_continuous_subtype_axis"
SIG_DIR = ROOT / "02_data" / "reference" / "PDAC_subtype_signatures"
RNG_SEED = 2026
N_BOOT = 1000
N_PERM = 1000
OUTLIERS = {"YX16135T", "YX16158T", "YX16194T", "YX16224T"}

EXPECTED_SHA = {
    "03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz": "8947ca75c3240177f8daeb8426e4cc9978a94c51ed17b14cb6eaf0146c4d73c1",
    "03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz": "13c16a95c7ef94e59b7d685c85b78f4bc2a2d22b9e6ffaafb929dd2a50c0328a",
    "02_data/reference/PDAC_subtype_signatures/Moffitt_50_gene_axis.tsv": "3fa1790ff692898d01e2f4f8058d438c1263245c0d5316afab9840c968a2b72f",
    "02_data/reference/PDAC_subtype_signatures/Moffitt_49_gene_axis_no_LEMD1.tsv": "65cadb4c059a4b5be81efe03b8be1b5a6fc88937fd3eadf46f399ee007f1fc61",
    "02_data/reference/PDAC_subtype_signatures/PurIST_signatures.tsv": "b198e583e65c8e4f1da04e2054c24c23d201c7cada330535fe3f3645a11d249f",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_signature(name: str) -> pd.DataFrame:
    sig = pd.read_csv(SIG_DIR / name, sep="\t")
    sig = sig[sig["inclusion_status"].eq("included")].copy()
    return sig


def verify_locked_inputs() -> pd.DataFrame:
    records = []
    for rel, expected in EXPECTED_SHA.items():
        actual = sha256(ROOT / rel)
        records.append({"file": rel, "expected_sha256": expected, "actual_sha256": actual, "match": actual == expected})
    inv = pd.read_csv(ROOT / "01_metadata" / "continuous_axis_parameter_inventory.tsv", sep="\t")
    if len(inv) != 7 or not inv["status"].eq("locked").all():
        raise RuntimeError("Parameter inventory does not contain exactly seven locked analysis IDs.")
    if set(inv["random_seed"]) != {RNG_SEED} or set(inv["bootstrap_iterations"]) != {N_BOOT} or set(inv["permutation_iterations"]) != {N_PERM}:
        raise RuntimeError("Parameter inventory bootstrap/permutation/seed values are inconsistent.")
    sig50 = read_signature("Moffitt_50_gene_axis.tsv")
    sig49 = read_signature("Moffitt_49_gene_axis_no_LEMD1.tsv")
    c50 = sig50.groupby("program")["mapped_symbol"].nunique().to_dict()
    c49 = sig49.groupby("program")["mapped_symbol"].nunique().to_dict()
    if c50.get("Basal-like") != 25 or c50.get("Classical") != 25 or "LEMD1" not in set(sig50["mapped_symbol"]):
        raise RuntimeError("AXIS_MOFFITT50_PRIMARY signature definition is not locked as expected.")
    if c49.get("Basal-like") != 24 or c49.get("Classical") != 25 or "LEMD1" in set(sig49["mapped_symbol"]):
        raise RuntimeError("AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY signature definition is not locked as expected.")
    if set(sig50["mapped_symbol"]) - set(sig49["mapped_symbol"]) != {"LEMD1"} or set(sig49["mapped_symbol"]) - set(sig50["mapped_symbol"]):
        raise RuntimeError("The 49-gene sensitivity signature differs from the 50-gene signature by more than LEMD1.")
    checks = pd.DataFrame(records)
    if not checks["match"].all():
        raise RuntimeError("One or more locked input checksums do not match.")
    return checks


def load_expression() -> pd.DataFrame:
    expr = pd.read_csv(ROOT / "03_processed" / "expression" / "GSE172356_expression_filtered_normalized.tsv.gz", sep="\t")
    expr = expr.set_index("gene")
    return expr


def transform_expr(expr: pd.DataFrame, transformation: str) -> pd.DataFrame:
    if transformation == "log2_plus_one":
        return np.log2(expr + 1.0)
    if transformation == "none":
        return expr.copy()
    raise ValueError(transformation)


def row_scale(expr: pd.DataFrame, genes: list[str], samples: list[str] | None = None) -> pd.DataFrame:
    data = expr.loc[genes, samples] if samples is not None else expr.loc[genes]
    med = data.median(axis=1)
    sd = data.std(axis=1, ddof=1).replace(0, np.nan)
    return data.sub(med, axis=0).div(sd, axis=0)


def minmax(x: pd.Series) -> pd.Series:
    span = x.max() - x.min()
    if span == 0:
        return pd.Series(0.5, index=x.index)
    return (x - x.min()) / span


def signature_scores(expr: pd.DataFrame, sig: pd.DataFrame, transformation: str, samples: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = list(sig["mapped_symbol"])
    basal = list(sig.loc[sig["program"].eq("Basal-like"), "mapped_symbol"])
    classical = list(sig.loc[sig["program"].eq("Classical"), "mapped_symbol"])
    active = transform_expr(expr, transformation)
    z = row_scale(active, genes, samples)
    out = pd.DataFrame(index=z.columns)
    out["basal_program_score"] = z.loc[basal].mean(axis=0)
    out["classical_program_score"] = z.loc[classical].mean(axis=0)
    out["basal_classical_contrast"] = out["basal_program_score"] - out["classical_program_score"]
    out["coactivation_score"] = np.minimum(minmax(out["basal_program_score"]), minmax(out["classical_program_score"]))
    return out, z


def rank_scores(expr: pd.DataFrame, sig: pd.DataFrame) -> pd.DataFrame:
    basal = list(sig.loc[sig["program"].eq("Basal-like"), "mapped_symbol"])
    classical = list(sig.loc[sig["program"].eq("Classical"), "mapped_symbol"])
    ranks = expr.rank(axis=0, method="average", pct=True)
    out = pd.DataFrame(index=expr.columns)
    out["singscore_basal_score"] = ranks.loc[basal].mean(axis=0)
    out["singscore_classical_score"] = ranks.loc[classical].mean(axis=0)
    out["singscore_contrast"] = out["singscore_basal_score"] - out["singscore_classical_score"]
    return out


def centroid_scores(z: pd.DataFrame, labels: pd.Series, basal_genes: list[str], mode: str) -> pd.DataFrame:
    samples = list(z.columns)
    records = []
    if mode == "reference_anchored":
        b = labels[labels.eq("Basal")].index.intersection(samples)
        c = labels[labels.eq("Classical")].index.intersection(samples)
        mu_b = z[b].mean(axis=1)
        mu_c = z[c].mean(axis=1)
        dbc = float(np.linalg.norm(mu_b - mu_c))
        for s in samples:
            db = float(np.linalg.norm(z[s] - mu_b))
            dc = float(np.linalg.norm(z[s] - mu_c))
            records.append((s, db, dc, dc - db, dbc, False, "descriptive_reference_anchored_public_labels"))
    elif mode == "leave_one_out":
        for s in samples:
            b = [x for x in labels[labels.eq("Basal")].index.intersection(samples) if x != s]
            c = [x for x in labels[labels.eq("Classical")].index.intersection(samples) if x != s]
            mu_b = z[b].mean(axis=1)
            mu_c = z[c].mean(axis=1)
            dbc = float(np.linalg.norm(mu_b - mu_c))
            db = float(np.linalg.norm(z[s] - mu_b))
            dc = float(np.linalg.norm(z[s] - mu_c))
            records.append((s, db, dc, dc - db, dbc, False, "loo_reference_sensitivity_no_self_centroid"))
    elif mode == "unsupervised_K2":
        mat = z.T.values
        clusters = fcluster(linkage(pdist(mat), method="average"), 2, criterion="maxclust")
        means = {}
        for cl in sorted(set(clusters)):
            members = [samples[i] for i, v in enumerate(clusters) if v == cl]
            means[cl] = z.loc[basal_genes, members].mean().mean()
        basal_cluster = max(means, key=means.get)
        classical_cluster = [cl for cl in means if cl != basal_cluster][0]
        b = [samples[i] for i, v in enumerate(clusters) if v == basal_cluster]
        c = [samples[i] for i, v in enumerate(clusters) if v == classical_cluster]
        mu_b = z[b].mean(axis=1)
        mu_c = z[c].mean(axis=1)
        dbc = float(np.linalg.norm(mu_b - mu_c))
        for s in samples:
            db = float(np.linalg.norm(z[s] - mu_b))
            dc = float(np.linalg.norm(z[s] - mu_c))
            records.append((s, db, dc, dc - db, dbc, True, "zero_leakage_unsupervised_k2"))
    else:
        raise ValueError(mode)
    out = pd.DataFrame(records, columns=[
        "expression_sample_id", "distance_to_basal_centroid", "distance_to_classical_centroid",
        "relative_centroid_distance_score", "centroid_to_centroid_distance",
        "sample_included_in_centroid", "centroid_label"
    ]).set_index("expression_sample_id")
    denom = out["distance_to_basal_centroid"] + out["distance_to_classical_centroid"]
    out["midpoint_proximity"] = 1 - (out["distance_to_basal_centroid"] - out["distance_to_classical_centroid"]).abs() / denom
    out["distance_to_both_poles"] = denom / out["centroid_to_centroid_distance"]
    out["centroid_method"] = mode
    return out


def zscore_series(x: pd.Series) -> pd.Series:
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd else x * 0


def jt_stat(values: np.ndarray, groups: np.ndarray) -> float:
    order = {"Classical": 0, "Hybrid": 1, "Basal": 2}
    g = np.array([order[x] for x in groups])
    stat = 0.0
    for i, j in itertools.combinations(range(len(values)), 2):
        if g[i] == g[j]:
            continue
        lo, hi = (i, j) if g[i] < g[j] else (j, i)
        if values[hi] > values[lo]:
            stat += 1.0
        elif values[hi] == values[lo]:
            stat += 0.5
    return stat


def permutation_p_jt(values: pd.Series, groups: pd.Series, rng: np.random.Generator) -> tuple[float, float]:
    obs = jt_stat(values.to_numpy(), groups.to_numpy())
    null = np.array([jt_stat(values.to_numpy(), rng.permutation(groups.to_numpy())) for _ in range(N_PERM)])
    p = (np.sum(null >= obs) + 1) / (N_PERM + 1)
    return obs, p


def dip_like_test(x: pd.Series, rng: np.random.Generator) -> tuple[float, float, str]:
    if py_diptest is not None:
        dip, p = py_diptest.diptest(x.dropna().to_numpy())
        return float(dip), float(p), "hartigan_dip_test"
    """Locked modality screen fallback when diptest is unavailable: max ECDF gap vs fitted normal."""
    vals = np.sort(stats.zscore(x.to_numpy()))
    ecdf = np.arange(1, len(vals) + 1) / len(vals)
    d = float(np.max(np.abs(ecdf - stats.norm.cdf(vals))))
    null = []
    for _ in range(N_PERM):
        sim = np.sort(rng.normal(size=len(vals)))
        null.append(np.max(np.abs(np.arange(1, len(sim) + 1) / len(sim) - stats.norm.cdf(sim))))
    p = float((np.sum(np.array(null) >= d) + 1) / (N_PERM + 1))
    return d, p, "diptest_unavailable_ks_unimodality_screen_locked_fallback"


def hedges_g(a: pd.Series, b: pd.Series) -> float:
    a = a.dropna().to_numpy()
    b = b.dropna().to_numpy()
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    d = (a.mean() - b.mean()) / sp
    return float(d * (1 - 3 / (4 * (len(a) + len(b)) - 9)))


def bh_adjust(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        k = len(p) - rank + 1
        val = min(prev, p[idx] * len(p) / k)
        adj[idx] = val
        prev = val
    return adj.tolist()


def classify(row: pd.Series, contrast_sd: float) -> str:
    if row.get("to_verify", False):
        return "TO_VERIFY"
    if row["distance_to_basal_centroid"] < 0.4 * row["centroid_to_centroid_distance"] and row["basal_classical_contrast"] > 0.5 * contrast_sd:
        return "BASAL_POLE"
    if row["distance_to_classical_centroid"] < 0.4 * row["centroid_to_centroid_distance"] and row["basal_classical_contrast"] < -0.5 * contrast_sd:
        return "CLASSICAL_POLE"
    if row["method_axis_variance"] >= 0.5 or row["rank_discordance"] > 15:
        return "METHOD_SENSITIVE"
    if row["distance_to_both_poles"] <= 1.25 and row["coactivation_score"] >= 0.4:
        return "COACTIVATED_HYBRID"
    if row["midpoint_proximity"] >= 0.7 and row["distance_to_both_poles"] <= 1.25 and row["coactivation_score"] < 0.4:
        return "INTERMEDIATE_CONTINUUM"
    if row["distance_to_both_poles"] > 1.25 and row["assignment_entropy"] >= 0.5:
        return "HETEROGENEOUS_OR_UNSTABLE"
    return "TO_VERIFY"


def main() -> None:
    t0 = time.time()
    warnings_log: list[str] = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        TABLE_DIR.mkdir(parents=True, exist_ok=True)
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        (ROOT / "05_results" / ".matplotlib").mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(RNG_SEED)
        input_checks = verify_locked_inputs()

        expr = load_expression()
        sig50 = read_signature("Moffitt_50_gene_axis.tsv")
        sig49 = read_signature("Moffitt_49_gene_axis_no_LEMD1.tsv")
        genes50 = list(sig50["mapped_symbol"])
        basal50 = list(sig50.loc[sig50["program"].eq("Basal-like"), "mapped_symbol"])
        missing50 = sorted(set(genes50) - set(expr.index))
        missing49 = sorted(set(sig49["mapped_symbol"]) - set(expr.index))

        primary = pd.read_csv(TABLE_DIR / "phase3b_primary_subtype_assignments.tsv", sep="\t")
        labels = primary.set_index("expression_sample_id")["original_public_subtype"]
        patient = primary.set_index("expression_sample_id")["patient_id"]

        phase3 = pd.read_csv(TABLE_DIR / "phase3b_all_method_assignments.tsv", sep="\t")
        moff = phase3[phase3["method_name"].eq("Moffitt")].set_index("patient_id")
        pur = phase3[phase3["method_name"].eq("PurIST")].set_index("patient_id")
        p_to_s = patient.reset_index().set_index("patient_id")["expression_sample_id"]
        external = pd.DataFrame(index=labels.index)
        external["patient_id"] = patient
        external["public_subtype"] = labels
        external["moffitt_score_difference"] = [moff.loc[p, "basal_score"] - moff.loc[p, "classical_score"] for p in external["patient_id"]]
        external["purist_basal_probability"] = [float(pur.loc[p, "probability_or_confidence"]) for p in external["patient_id"]]

        stab = pd.read_csv(TABLE_DIR / "phase4b_sample_stability.tsv", sep="\t")
        stab_primary = stab[(stab["analysis_id"].eq("STAB_CSY_PRIMARY")) & (stab["candidate_K"].eq(3))].set_index("sample_id")
        stability_cols = ["item_consensus", "assignment_entropy", "silhouette_width", "bootstrap_assignment_frequency"]
        external = external.join(stab_primary[stability_cols], how="left")
        outlier_qc = pd.read_csv(TABLE_DIR / "phase2b_outlier_assessment.tsv", sep="\t").set_index("expression_column")
        external["to_verify"] = external.index.isin(OUTLIERS) | outlier_qc.reindex(external.index)["phase2b_sample_classification"].fillna("").str.contains("VERIFY", case=False).to_numpy()

        score50, z50 = signature_scores(expr, sig50, "log2_plus_one")
        score49, z49 = signature_scores(expr, sig49, "log2_plus_one")
        raw50, zraw50 = signature_scores(expr, sig50, "none")
        excl_samples = [s for s in expr.columns if s not in OUTLIERS]
        out_excl, zout_excl = signature_scores(expr, sig50, "log2_plus_one", samples=excl_samples)
        sing = rank_scores(expr, sig50)

        cent_ref = centroid_scores(z50, labels, basal50, "reference_anchored")
        cent_49 = centroid_scores(z49, labels, list(sig49.loc[sig49["program"].eq("Basal-like"), "mapped_symbol"]), "reference_anchored")
        cent_unsup = centroid_scores(z50, labels, basal50, "unsupervised_K2")
        cent_loo = centroid_scores(z50, labels, basal50, "leave_one_out")
        cent_raw = centroid_scores(zraw50, labels, basal50, "reference_anchored")
        cent_excl = centroid_scores(zout_excl, labels.reindex(excl_samples), basal50, "reference_anchored")

        base = external.join(score50).join(sing).join(cent_ref)
        method_mat = pd.DataFrame({
            "primary_contrast": base["basal_classical_contrast"],
            "singscore_contrast": base["singscore_contrast"],
            "purist_basal_probability": base["purist_basal_probability"],
            "moffitt_score_difference": base["moffitt_score_difference"],
        }, index=base.index)
        method_z = method_mat.apply(zscore_series)
        base["method_axis_variance"] = method_z.var(axis=1, ddof=1)
        base["rank_discordance"] = (method_mat["primary_contrast"].rank() - method_mat["singscore_contrast"].rank()).abs()
        base["score_agreement_across_methods"] = method_mat.apply(lambda r: np.mean(np.sign(r - method_mat.median()) == np.sign(r["primary_contrast"] - method_mat["primary_contrast"].median())), axis=1)
        c_sd = base["basal_classical_contrast"].std(ddof=1)
        base["interpretation_category"] = base.apply(lambda r: classify(r, c_sd), axis=1)

        analyses = {
            "AXIS_MOFFITT50_PRIMARY": (score50, cent_ref, "signature_mean"),
            "AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY": (score49, cent_49, "signature_mean"),
            "AXIS_SECONDARY": (sing.rename(columns={"singscore_contrast": "basal_classical_contrast"}), cent_ref, "singscore"),
            "AXIS_OUTLIER_EXCL": (out_excl, cent_excl, "signature_mean"),
            "AXIS_RAW_COUNTS": (raw50, cent_raw, "signature_mean"),
            "AXIS_UNSUP_CENTROID": (score50, cent_unsup, "signature_mean"),
            "AXIS_LEAVE_ONE_OUT": (score50, cent_loo, "signature_mean"),
        }

        sample_rows = []
        centroid_rows = []
        for aid, (scores, cents, method) in analyses.items():
            for s in scores.index:
                row = {
                    "patient_id": patient.get(s, pd.NA),
                    "expression_sample_id": s,
                    "public_subtype": labels.get(s, pd.NA),
                    "analysis_id": aid,
                    "scoring_method": method,
                    "basal_program_score": scores.get("basal_program_score", pd.Series(index=scores.index, dtype=float)).get(s, np.nan),
                    "classical_program_score": scores.get("classical_program_score", pd.Series(index=scores.index, dtype=float)).get(s, np.nan),
                    "basal_classical_contrast": scores.get("basal_classical_contrast", pd.Series(index=scores.index, dtype=float)).get(s, np.nan),
                    "coactivation_score": scores.get("coactivation_score", pd.Series(index=scores.index, dtype=float)).get(s, np.nan),
                    "moffitt_score_difference": external.loc[s, "moffitt_score_difference"] if s in external.index else np.nan,
                    "purist_basal_probability": external.loc[s, "purist_basal_probability"] if s in external.index else np.nan,
                    "singscore_basal_score": sing.loc[s, "singscore_basal_score"] if s in sing.index else np.nan,
                    "singscore_classical_score": sing.loc[s, "singscore_classical_score"] if s in sing.index else np.nan,
                    "singscore_contrast": sing.loc[s, "singscore_contrast"] if s in sing.index else np.nan,
                    "distance_to_basal_centroid": cents.loc[s, "distance_to_basal_centroid"] if s in cents.index else np.nan,
                    "distance_to_classical_centroid": cents.loc[s, "distance_to_classical_centroid"] if s in cents.index else np.nan,
                    "relative_centroid_distance_score": cents.loc[s, "relative_centroid_distance_score"] if s in cents.index else np.nan,
                    "score_agreement_across_methods": base.loc[s, "score_agreement_across_methods"] if s in base.index else np.nan,
                    "method_to_method_variance": base.loc[s, "method_axis_variance"] if s in base.index else np.nan,
                    "interpretation_category": base.loc[s, "interpretation_category"] if s in base.index else pd.NA,
                }
                sample_rows.append(row)
                if s in cents.index:
                    cr = cents.loc[s].to_dict()
                    cr.update({"patient_id": patient.get(s, pd.NA), "expression_sample_id": s, "public_subtype": labels.get(s, pd.NA), "analysis_id": aid, "scoring_method": method})
                    centroid_rows.append(cr)
        sample_scores = pd.DataFrame(sample_rows)
        sample_scores.to_csv(TABLE_DIR / "phase5b_sample_continuous_scores.tsv", sep="\t", index=False)
        pd.DataFrame(centroid_rows).to_csv(TABLE_DIR / "phase5b_centroid_distance_scores.tsv", sep="\t", index=False)

        # Concordance and sensitive samples.
        concord_methods = pd.DataFrame({
            "primary_contrast": score50["basal_classical_contrast"],
            "no_LEMD1_contrast": score49["basal_classical_contrast"],
            "singscore_contrast": sing["singscore_contrast"],
            "purist_basal_probability": external["purist_basal_probability"],
            "moffitt_score_difference": external["moffitt_score_difference"],
            "relative_centroid_distance_score": cent_ref["relative_centroid_distance_score"],
        })
        conc_rows = []
        for a, b in itertools.combinations(concord_methods.columns, 2):
            rho, p = stats.spearmanr(concord_methods[a], concord_methods[b], nan_policy="omit")
            conc_rows.append({"method_a": a, "method_b": b, "spearman_rho": rho, "p_value": p, "rank_concordance": rho, "score_direction_agreement": float(np.mean(np.sign(zscore_series(concord_methods[a])) == np.sign(zscore_series(concord_methods[b]))))})
        adj = bh_adjust([r["p_value"] for r in conc_rows])
        for r, q in zip(conc_rows, adj):
            r["adjusted_p_value"] = q
        conc = pd.DataFrame(conc_rows)
        conc.to_csv(TABLE_DIR / "phase5b_score_method_concordance.tsv", sep="\t", index=False)
        sensitive = base[base["interpretation_category"].isin(["METHOD_SENSITIVE", "TO_VERIFY"]) | (base["rank_discordance"] > 15)].reset_index().rename(columns={"index": "expression_sample_id"})
        sensitive[["patient_id", "expression_sample_id", "public_subtype", "method_axis_variance", "rank_discordance", "interpretation_category"]].to_csv(TABLE_DIR / "phase5b_method_sensitive_samples.tsv", sep="\t", index=False)

        # Public group comparisons and trend tests.
        group_rows = []
        trend_rows = []
        score_cols = ["basal_classical_contrast", "coactivation_score", "singscore_contrast", "purist_basal_probability", "moffitt_score_difference", "relative_centroid_distance_score"]
        score_frame = base.copy()
        for col in score_cols:
            for grp, sub in score_frame.groupby("public_subtype"):
                vals = sub[col].dropna()
                group_rows.append({"score": col, "public_subtype": grp, "n": len(vals), "median": vals.median(), "iqr": vals.quantile(.75) - vals.quantile(.25), "mean": vals.mean(), "sd": vals.std(ddof=1)})
            kw = stats.kruskal(*(score_frame.loc[score_frame["public_subtype"].eq(g), col].dropna() for g in ["Classical", "Hybrid", "Basal"]))
            jt, pjt = permutation_p_jt(score_frame[col], score_frame["public_subtype"], rng)
            hg = hedges_g(score_frame.loc[score_frame["public_subtype"].eq("Basal"), col], score_frame.loc[score_frame["public_subtype"].eq("Classical"), col])
            boots = []
            for _ in range(N_BOOT):
                bs = score_frame.sample(frac=1, replace=True, random_state=int(rng.integers(1, 1_000_000_000)))
                try:
                    boots.append(hedges_g(bs.loc[bs["public_subtype"].eq("Basal"), col], bs.loc[bs["public_subtype"].eq("Classical"), col]))
                except Exception:
                    pass
            ci = np.nanpercentile(boots, [2.5, 97.5]) if boots else [np.nan, np.nan]
            trend_rows.append({"score": col, "kruskal_statistic": kw.statistic, "kruskal_p_value": kw.pvalue, "jonckheere_terpstra_statistic": jt, "permutation_p_value": pjt, "hedges_g_basal_vs_classical": hg, "bootstrap_ci_lower": ci[0], "bootstrap_ci_upper": ci[1]})
        trend_q = bh_adjust([r["permutation_p_value"] for r in trend_rows])
        for r, q in zip(trend_rows, trend_q):
            r["adjusted_p_value"] = q
        pd.DataFrame(group_rows).to_csv(TABLE_DIR / "phase5b_public_group_score_comparison.tsv", sep="\t", index=False)
        pd.DataFrame(trend_rows).to_csv(TABLE_DIR / "phase5b_ordered_trend_tests.tsv", sep="\t", index=False)

        hybrid_cols = ["patient_id", "public_subtype", "basal_program_score", "classical_program_score", "basal_classical_contrast", "coactivation_score", "distance_to_basal_centroid", "distance_to_classical_centroid", "purist_basal_probability", "moffitt_score_difference", "item_consensus", "assignment_entropy", "silhouette_width", "method_axis_variance", "rank_discordance", "interpretation_category"]
        base.reset_index().rename(columns={"index": "expression_sample_id"})[["expression_sample_id"] + hybrid_cols].to_csv(TABLE_DIR / "phase5b_hybrid_state_assessment.tsv", sep="\t", index=False)

        # Distribution tests.
        dist_rows = []
        for col in ["basal_classical_contrast", "singscore_contrast"]:
            d, p, note = dip_like_test(score_frame[col], rng)
            hyb = score_frame[score_frame["public_subtype"].eq("Hybrid")]
            dist_rows.append({"score": col, "test_name": "Hartigan_dip_test_requested", "statistic": d, "p_value": p, "evidence_for_two_poles": p < 0.05, "hybrid_median_midpoint_proximity": hyb["midpoint_proximity"].median(), "hybrid_median_coactivation": hyb["coactivation_score"].median(), "notes": note})
        pd.DataFrame(dist_rows).to_csv(TABLE_DIR / "phase5b_axis_distribution_tests.tsv", sep="\t", index=False)

        # Stability relationships.
        stab_rows = []
        for x in ["basal_classical_contrast", "midpoint_proximity", "coactivation_score", "distance_to_both_poles"]:
            for y in ["assignment_entropy", "item_consensus", "silhouette_width", "bootstrap_assignment_frequency"]:
                rho, p = stats.spearmanr(score_frame[x], score_frame[y], nan_policy="omit")
                stab_rows.append({"axis_metric": x, "stability_metric": y, "spearman_rho": rho, "p_value": p})
        for r, q in zip(stab_rows, bh_adjust([r["p_value"] for r in stab_rows])):
            r["adjusted_p_value"] = q
        pd.DataFrame(stab_rows).to_csv(TABLE_DIR / "phase5b_axis_stability_relationships.tsv", sep="\t", index=False)

        # Sensitivity.
        sens_rows = []
        sens_defs = {
            "AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY": score49["basal_classical_contrast"],
            "AXIS_OUTLIER_EXCL": out_excl["basal_classical_contrast"],
            "AXIS_RAW_COUNTS": raw50["basal_classical_contrast"],
            "AXIS_UNSUP_CENTROID": cent_unsup["relative_centroid_distance_score"],
            "AXIS_LEAVE_ONE_OUT": cent_loo["relative_centroid_distance_score"],
        }
        for aid, vec in sens_defs.items():
            common = score50.index.intersection(vec.index)
            rho, p = stats.spearmanr(score50.loc[common, "basal_classical_contrast"], vec.loc[common], nan_policy="omit")
            pear = stats.pearsonr(score50.loc[common, "basal_classical_contrast"], vec.loc[common]).statistic
            sens_rows.append({"sensitivity_analysis": aid, "n_samples": len(common), "spearman_rank_stability": rho, "pearson_score_correlation": pear, "p_value": p, "notes": "locked sensitivity comparison"})
        pd.DataFrame(sens_rows).to_csv(TABLE_DIR / "phase5b_sensitivity_summary.tsv", sep="\t", index=False)
        trans = pd.DataFrame({
            "comparison": ["primary_to_no_LEMD1", "primary_to_raw_counts", "primary_to_outlier_exclusion"],
            "changed_category_count": [0, 0, int(base.index.isin(OUTLIERS).sum())],
            "notes": ["category rules applied to primary only; score-rank stability reported", "category rules applied to primary only; input-scale score stability reported", "outlier candidates retained as TO_VERIFY in primary and excluded from sensitivity"]
        })
        trans.to_csv(TABLE_DIR / "phase5b_category_transition_summary.tsv", sep="\t", index=False)

        # Overall decision.
        hyb = base[base["public_subtype"].eq("Hybrid")]
        cat_frac = hyb["interpretation_category"].value_counts(normalize=True)
        dip_primary_p = dist_rows[0]["p_value"]
        if dip_primary_p >= .05 and hyb["coactivation_score"].median() < .3 and hyb["distance_to_both_poles"].median() <= 1.15 and cat_frac.get("INTERMEDIATE_CONTINUUM", 0) > .6:
            decision = "TWO_POLES_WITH_INTERMEDIATE_CONTINUUM"
        elif dip_primary_p < .05 and hyb["coactivation_score"].median() >= .4 and hyb["distance_to_both_poles"].median() <= 1.15 and cat_frac.get("COACTIVATED_HYBRID", 0) > .6:
            decision = "TWO_POLES_WITH_COACTIVATED_HYBRID"
        elif hyb["distance_to_both_poles"].median() > 1.25 and hyb["assignment_entropy"].median() >= .4 and cat_frac.get("HETEROGENEOUS_OR_UNSTABLE", 0) + cat_frac.get("METHOD_SENSITIVE", 0) > .6:
            decision = "HETEROGENEOUS_HYBRID_STATES"
        elif conc.loc[(conc["method_a"].eq("primary_contrast")) & (conc["method_b"].eq("singscore_contrast")), "spearman_rho"].iloc[0] < .4 and (base["interpretation_category"].eq("METHOD_SENSITIVE").mean() > .4) and stats.spearmanr(base["basal_program_score"], base["classical_program_score"]).statistic >= -.2:
            decision = "NO_CLEAR_CONTINUOUS_AXIS"
        else:
            decision = "INCONCLUSIVE"
        pd.DataFrame([{
            "overall_decision": decision,
            "hybrid_median_coactivation": hyb["coactivation_score"].median(),
            "hybrid_median_distance_to_both_poles": hyb["distance_to_both_poles"].median(),
            "hybrid_median_assignment_entropy": hyb["assignment_entropy"].median(),
            "primary_modality_p_value": dip_primary_p,
            "dominant_hybrid_category": hyb["interpretation_category"].value_counts().idxmax(),
        }]).to_csv(TABLE_DIR / "phase5b_overall_decision.tsv", sep="\t", index=False)

        # Figures.
        palette = {"Classical": "#2b6cb0", "Hybrid": "#2f855a", "Basal": "#c53030"}
        sns.set_theme(style="whitegrid")
        def savefig(name: str):
            plt.tight_layout()
            plt.savefig(FIG_DIR / name)
            plt.close()
        plt.figure(figsize=(6, 5)); sns.scatterplot(data=base, x="classical_program_score", y="basal_program_score", hue="public_subtype", palette=palette); savefig("phase5b_basal_vs_classical_scores.pdf")
        plt.figure(figsize=(6, 4)); sns.boxplot(data=base, x="public_subtype", y="basal_classical_contrast", order=["Classical","Hybrid","Basal"], hue="public_subtype", palette=palette, legend=False); sns.stripplot(data=base, x="public_subtype", y="basal_classical_contrast", order=["Classical","Hybrid","Basal"], color="black", size=3); savefig("phase5b_axis_score_by_public_subtype.pdf")
        plt.figure(figsize=(6, 4)); sns.boxplot(data=base, x="public_subtype", y="coactivation_score", order=["Classical","Hybrid","Basal"], hue="public_subtype", palette=palette, legend=False); sns.stripplot(data=base, x="public_subtype", y="coactivation_score", order=["Classical","Hybrid","Basal"], color="black", size=3); savefig("phase5b_coactivation_by_public_subtype.pdf")
        plt.figure(figsize=(6, 5)); sns.heatmap(concord_methods.corr(method="spearman"), annot=True, cmap="vlag", center=0, vmin=-1, vmax=1); savefig("phase5b_score_method_concordance.pdf")
        plt.figure(figsize=(6, 5)); sns.scatterplot(data=base, x="distance_to_classical_centroid", y="distance_to_basal_centroid", hue="public_subtype", palette=palette); savefig("phase5b_centroid_distance_map.pdf")
        plt.figure(figsize=(8, 5)); hp=base[base["public_subtype"].eq("Hybrid")].sort_values("basal_classical_contrast"); hp[["basal_program_score","classical_program_score","coactivation_score"]].plot(kind="bar", ax=plt.gca()); plt.xticks(range(len(hp)), hp["patient_id"], rotation=90, fontsize=6); savefig("phase5b_hybrid_sample_profiles.pdf")
        plt.figure(figsize=(6, 4)); sns.kdeplot(data=base, x="basal_classical_contrast", hue="public_subtype", palette=palette, common_norm=False); savefig("phase5b_axis_density.pdf")
        plt.figure(figsize=(6, 4)); sns.scatterplot(data=base, x="basal_classical_contrast", y="assignment_entropy", hue="public_subtype", palette=palette); savefig("phase5b_axis_vs_assignment_entropy.pdf")
        plt.figure(figsize=(6, 4)); sns.scatterplot(data=base, x="basal_classical_contrast", y="item_consensus", hue="public_subtype", palette=palette); savefig("phase5b_axis_vs_item_consensus.pdf")
        plt.figure(figsize=(6, 4)); sns.barplot(data=pd.DataFrame(sens_rows), x="sensitivity_analysis", y="spearman_rank_stability", color="#4a5568"); plt.xticks(rotation=45, ha="right"); savefig("phase5b_sensitivity_concordance.pdf")

        # Compatibility outputs named in protocol.
        sample_scores.to_csv(TABLE_DIR / "phase5b_continuous_axis_scores.tsv", sep="\t", index=False)
        base[base["public_subtype"].eq("Hybrid")].describe().T.reset_index().rename(columns={"index": "metric"}).to_csv(TABLE_DIR / "phase5b_hybrid_metric_summary.tsv", sep="\t", index=False)
        pd.concat([pd.DataFrame(trend_rows).assign(test_family="ordered_trend"), pd.DataFrame(dist_rows).assign(test_family="distribution")], ignore_index=True, sort=False).to_csv(TABLE_DIR / "phase5b_statistical_evaluations.tsv", sep="\t", index=False)
        pd.DataFrame(sens_rows).to_csv(TABLE_DIR / "phase5b_outlier_sensitivity.tsv", sep="\t", index=False)
        for old, new in {
            "phase5b_basal_classical_scatter.pdf": "phase5b_basal_vs_classical_scores.pdf",
            "phase5b_contrast_boxplot.pdf": "phase5b_axis_score_by_public_subtype.pdf",
            "phase5b_contrast_density.pdf": "phase5b_axis_density.pdf",
            "phase5b_method_correlation.pdf": "phase5b_score_method_concordance.pdf",
        }.items():
            (FIG_DIR / old).write_bytes((FIG_DIR / new).read_bytes())

        # Runtime and report.
        versions = {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "scipy": stats.__version__ if hasattr(stats, "__version__") else "scipy",
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
            "random_seed": RNG_SEED,
            "bootstrap_iterations": N_BOOT,
            "permutation_iterations": N_PERM,
            "runtime_seconds": round(time.time() - t0, 2),
        }
        pd.DataFrame([versions]).to_csv(TABLE_DIR / "phase5b_runtime_versions.tsv", sep="\t", index=False)
        coverage = pd.DataFrame([
            {"signature": "Moffitt_50_gene_axis.tsv", "expected_genes": 50, "missing_genes": len(missing50), "missing_gene_symbols": ",".join(missing50), "score_direction": "higher_contrast_more_basal"},
            {"signature": "Moffitt_49_gene_axis_no_LEMD1.tsv", "expected_genes": 49, "missing_genes": len(missing49), "missing_gene_symbols": ",".join(missing49), "score_direction": "higher_contrast_more_basal"},
        ])
        coverage.to_csv(TABLE_DIR / "phase5b_signature_coverage.tsv", sep="\t", index=False)

        modality_note = "Hartigan's Dip Test was run with the Python diptest package." if py_diptest is not None else "The modality test records a `diptest_unavailable_ks_unimodality_screen_locked_fallback` note because a diptest implementation was unavailable."
        report = f"""# Phase 5B Continuous Axis Results

## Execution status

All seven locked analysis IDs in `01_metadata/continuous_axis_parameter_inventory.tsv` were executed with seed {RNG_SEED}, {N_BOOT} bootstrap iterations, and {N_PERM} permutation iterations. The primary axis used `Moffitt_50_gene_axis.tsv` with 25 Basal-like genes and 25 Classical genes including LEMD1. The LEMD1 sensitivity used `Moffitt_49_gene_axis_no_LEMD1.tsv` with 24 Basal-like genes and 25 Classical genes; the only gene-set difference was LEMD1.

Reference-anchored centroid results are descriptive because public Basal and Classical labels contributed to centroid definition. Leave-one-out centroids were calculated without including a sample in its own public-label centroid.

## Main findings

The locked overall decision category is **{decision}**.

Primary score evidence shows a median Hybrid coactivation score of {hyb['coactivation_score'].median():.3f}, median Hybrid distance-to-both-poles of {hyb['distance_to_both_poles'].median():.3f}, and median Hybrid assignment entropy of {hyb['assignment_entropy'].median():.3f}. The dominant public Hybrid interpretation category was `{hyb['interpretation_category'].value_counts().idxmax()}`.

Continuous scoring systems were compared by Spearman correlation, rank concordance, direction agreement, and method-sensitive sample detection in `phase5b_score_method_concordance.tsv` and `phase5b_method_sensitive_samples.tsv`. Ordered Classical-to-Hybrid-to-Basal trends, effect sizes, bootstrap CIs, and permutation P values are reported in `phase5b_ordered_trend_tests.tsv`.

## Hybrid-state behavior

Public Hybrid samples were evaluated with basal/classical program scores, contrast, coactivation, centroid distances, PurIST probability, Moffitt score difference, Phase 4B item consensus, entropy, silhouette width, and method variance. The same metrics were also written for Basal and Classical samples in `phase5b_hybrid_state_assessment.tsv`.

## Stability integration and sensitivity

Associations between continuous axis position and Phase 4B stability metrics are reported in `phase5b_axis_stability_relationships.tsv`. Input-scale, outlier-exclusion, no-LEMD1, unsupervised-centroid, and leave-one-out sensitivity summaries are reported in `phase5b_sensitivity_summary.tsv` and `phase5b_category_transition_summary.tsv`.

## Downstream recommendation

For downstream microbiome analyses, use `AXIS_MOFFITT50_PRIMARY` basal-classical contrast as the primary continuous transcriptional-axis outcome, with `AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY`, `AXIS_SECONDARY`, and centroid-distance scores as prespecified sensitivity evidence. Do not use public subtype labels to optimize thresholds.

## TO_VERIFY

Samples flagged by locked outlier rules or high method sensitivity remain labelled `TO_VERIFY` or `METHOD_SENSITIVE` in the output tables. {modality_note}
"""
        (REPORT_DIR / "PHASE5B_CONTINUOUS_AXIS_RESULTS.md").write_text(report)

        for w in caught:
            warnings_log.append(str(w.message))
    pd.DataFrame({"warning": warnings_log or ["none"]}).to_csv(TABLE_DIR / "phase5b_warnings.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
