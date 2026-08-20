#!/usr/bin/env python3
"""Execute locked Phase 6C PRJNA719915 tumor microbiome preprocessing."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MPL_DIR = ROOT / ".cache/matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(MPL_DIR)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr


AUDITED = ROOT / "03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz"
CROSSWALK = ROOT / "01_metadata/microbiome_sample_crosswalk.tsv"
PARAMS = ROOT / "01_metadata/microbiome_preprocessing_parameter_inventory.tsv"
PHASE6B_LOCK = ROOT / "04_analysis/04_microbiome_qc/PHASE6B_MICROBIOME_METHOD_LOCK.md"
PROTOCOL = ROOT / "09_docs/methods/PDAC_microbiome_preprocessing_protocol.md"
PHASE6A_SAMPLE_QC = ROOT / "05_results/tables/phase6a_microbiome_sample_qc.tsv"
PHASE6A_CONTAM = ROOT / "05_results/tables/phase6a_potential_contaminant_flags.tsv"
OUT_DIR = ROOT / "03_processed/microbiome"
SENS_DIR = OUT_DIR / "sensitivity"
TABLE_DIR = ROOT / "05_results/tables"
FIGURE_DIR = ROOT / "05_results/figures"
WARNINGS_PATH = TABLE_DIR / "phase6c_warnings.tsv"
SUMMARY_SCRIPT = ROOT / "06_scripts/python/09_summarize_phase6c_microbiome.py"

PRIMARY_FILTERED = OUT_DIR / "PRJNA719915_genus_primary_filtered.tsv.gz"
PRIMARY_CLR = OUT_DIR / "PRJNA719915_genus_primary_CLR.tsv.gz"
PRIMARY_DIST = OUT_DIR / "PRJNA719915_primary_aitchison_distance.tsv.gz"

PRIMARY_PSEUDOCOUNT = 0.889651
MIN_NONZERO_LOCKED = 1.77930272379619
RANDOM_SEED = 2026
EXPECTED_PRIMARY_FEATURES = 122
EXTREME_SAMPLES = ["Basal-like1", "Hybrid18", "Hybrid23"]

HIGH_RISK = {
    "Elizabethkingia",
    "Delftia",
    "Brevundimonas",
    "Comamonas",
    "Caulobacter",
    "Ralstonia",
}
MODERATE_RISK = {
    "Paraburkholderia",
    "Mesorhizobium",
    "Novosphingobium",
    "Dechloromonas",
    "Sphingopyxis",
    "Herbaspirillum",
}
BIO_PLAUSIBLE = {
    "Pseudomonas",
    "Acinetobacter",
    "Burkholderia",
    "Stenotrophomonas",
    "Sphingomonas",
    "Rhizobium",
    "Cupriavidus",
    "Methylobacterium",
    "Bradyrhizobium",
}


@dataclass
class Representation:
    analysis_id: str
    role: str
    filtered: pd.DataFrame
    transformed: pd.DataFrame
    distance: pd.DataFrame
    transform_name: str
    distance_metric: str
    pseudocount: str
    prevalence_threshold: float
    abundance_threshold: float
    contaminant_policy: str
    sample_exclusion_policy: str
    output_files: list[str]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", compression="gzip", float_format="%.12g")


def locked_file_consistency() -> pd.DataFrame:
    for path in [PHASE6B_LOCK, PROTOCOL, PARAMS]:
        require(path.exists(), f"Missing locked input: {path}")
    lock_text = PHASE6B_LOCK.read_text()
    protocol_text = PROTOCOL.read_text()
    param = pd.read_csv(PARAMS, sep="\t")

    checks = [
        ("method_lock_primary_prevalence", "at least 20% of samples" in lock_text),
        ("method_lock_primary_detection", "abundance $> 0.0$" in lock_text),
        ("method_lock_primary_pseudocount", "`0.889651`" in lock_text),
        ("method_lock_primary_clr", "Centered Log-Ratio" in lock_text),
        ("method_lock_primary_aitchison", "Aitchison distance" in lock_text),
        ("protocol_primary_prevalence", "less than 20%" in protocol_text and ("\\ge 13" in protocol_text or ">= 13" in protocol_text)),
        ("protocol_primary_detection", "strictly greater than 0.0" in protocol_text),
        ("protocol_primary_pseudocount", "0.889651" in protocol_text),
        ("protocol_primary_clr", "Centered Log-Ratio" in protocol_text),
        ("protocol_primary_aitchison", "Aitchison distance" in protocol_text),
    ]
    primary = param.loc[param["analysis_id"] == "MICRO_PRIMARY"].iloc[0]
    checks.extend(
        [
            ("inventory_primary_prevalence", math.isclose(float(primary["prevalence_threshold"]), 0.2)),
            ("inventory_primary_abundance", math.isclose(float(primary["abundance_threshold"]), 0.0)),
            ("inventory_primary_pseudocount", primary["pseudocount_method"] == "fixed_0.889651"),
            ("inventory_primary_transform", primary["transformation"] == "centered_log_ratio"),
            ("inventory_primary_distance", primary["distance_metric"] == "aitchison"),
        ]
    )
    failed = [name for name, ok in checks if not ok]
    require(not failed, "Locked filtering, pseudocount, or transformation parameters are inconsistent: " + ", ".join(failed))
    require(set(param["status"]) == {"locked"}, "Parameter inventory contains unlocked rows")
    return param


def prevalence_filter(df: pd.DataFrame, prevalence_threshold: float, abundance_threshold: float) -> pd.DataFrame:
    min_samples = math.ceil(prevalence_threshold * df.shape[1] - 1e-12)
    detected = df > abundance_threshold
    keep = detected.sum(axis=1) >= min_samples
    return df.loc[keep].copy()


def clr(df: pd.DataFrame, pseudocount: float) -> pd.DataFrame:
    shifted = df + pseudocount
    log_df = np.log(shifted)
    return log_df.sub(log_df.mean(axis=0), axis=1)


def robust_clr(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(0.0, index=df.index, columns=df.columns)
    for sample in df.columns:
        values = df[sample].astype(float)
        nonzero = values > 0
        if nonzero.sum() == 0:
            continue
        logs = np.log(values.loc[nonzero])
        centered = logs - logs.mean()
        out.loc[nonzero, sample] = centered
    return out


def euclidean_distance(transformed: pd.DataFrame) -> pd.DataFrame:
    distances = squareform(pdist(transformed.T, metric="euclidean"))
    return pd.DataFrame(distances, index=transformed.columns, columns=transformed.columns)


def jaccard_distance(binary: pd.DataFrame) -> pd.DataFrame:
    distances = squareform(pdist(binary.T.astype(bool), metric="jaccard"))
    return pd.DataFrame(distances, index=binary.columns, columns=binary.columns)


def pcoa(distance: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    d = distance.to_numpy(dtype=float)
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    pos = np.maximum(vals[:n_components], 0)
    coords = vecs[:, :n_components] * np.sqrt(pos)
    cols = [f"PCoA{i + 1}" for i in range(n_components)]
    return pd.DataFrame(coords, index=distance.index, columns=cols)


def validate_distance(distance: pd.DataFrame, expected_n: int | None = None) -> None:
    d = distance.to_numpy(dtype=float)
    if expected_n is not None:
        require(distance.shape == (expected_n, expected_n), f"Distance matrix shape {distance.shape}, expected {expected_n} x {expected_n}")
    require(np.isfinite(d).all(), "Distance matrix contains missing or infinite values")
    require((d >= -1e-12).all(), "Distance matrix contains negative values")
    require(np.allclose(d, d.T, atol=1e-10), "Distance matrix is not symmetric")
    require(np.allclose(np.diag(d), 0.0, atol=1e-10), "Distance matrix diagonal is not zero")


def save_representation(rep: Representation) -> None:
    base = SENS_DIR / rep.analysis_id
    if rep.analysis_id == "MICRO_PRIMARY":
        write_tsv_gz(rep.filtered, PRIMARY_FILTERED)
        write_tsv_gz(rep.transformed, PRIMARY_CLR)
        write_tsv_gz(rep.distance, PRIMARY_DIST)
        rep.output_files.extend([str(PRIMARY_FILTERED.relative_to(ROOT)), str(PRIMARY_CLR.relative_to(ROOT)), str(PRIMARY_DIST.relative_to(ROOT))])
    else:
        filtered_path = base.with_name(base.name + "_filtered.tsv.gz")
        transform_path = base.with_name(base.name + f"_{rep.transform_name}.tsv.gz")
        dist_path = base.with_name(base.name + f"_{rep.distance_metric}_distance.tsv.gz")
        write_tsv_gz(rep.filtered, filtered_path)
        write_tsv_gz(rep.transformed, transform_path)
        write_tsv_gz(rep.distance, dist_path)
        rep.output_files.extend([str(filtered_path.relative_to(ROOT)), str(transform_path.relative_to(ROOT)), str(dist_path.relative_to(ROOT))])


def make_rep(
    analysis_id: str,
    source: pd.DataFrame,
    prevalence_threshold: float,
    abundance_threshold: float,
    pseudocount: float | None,
    transform_name: str,
    distance_metric: str,
    contaminant_policy: str,
    sample_exclusion_policy: str,
    role: str = "sensitivity",
) -> Representation:
    filtered = prevalence_filter(source, prevalence_threshold, abundance_threshold)
    if transform_name == "centered_log_ratio":
        require(pseudocount is not None, f"{analysis_id} requires a pseudocount")
        transformed = clr(filtered, pseudocount)
    elif transform_name == "robust_clr":
        transformed = robust_clr(filtered)
    elif transform_name == "presence_absence":
        transformed = (filtered > abundance_threshold).astype(int)
    else:
        raise ValueError(f"Unsupported transform: {transform_name}")

    if distance_metric == "aitchison":
        distance = euclidean_distance(transformed)
    elif distance_metric == "jaccard":
        distance = jaccard_distance(transformed)
    else:
        raise ValueError(f"Unsupported distance metric: {distance_metric}")

    if transform_name == "centered_log_ratio":
        require(np.isfinite(transformed.to_numpy()).all(), f"{analysis_id} CLR contains invalid values")
        require(np.allclose(transformed.sum(axis=0), 0.0, atol=1e-8), f"{analysis_id} CLR sample sums are not zero")
    validate_distance(distance, expected_n=source.shape[1])
    return Representation(
        analysis_id=analysis_id,
        role=role,
        filtered=filtered,
        transformed=transformed,
        distance=distance,
        transform_name=transform_name,
        distance_metric=distance_metric,
        pseudocount="none" if pseudocount is None else f"{pseudocount:g}",
        prevalence_threshold=prevalence_threshold,
        abundance_threshold=abundance_threshold,
        contaminant_policy=contaminant_policy,
        sample_exclusion_policy=sample_exclusion_policy,
        output_files=[],
    )


def matrix_inventory(reps: list[Representation], package_versions: dict[str, str], original_sha: str) -> pd.DataFrame:
    rows = []
    for rep in reps:
        values = rep.filtered.to_numpy(dtype=float)
        transformed = rep.transformed.to_numpy(dtype=float)
        clr_col_sum = np.nan
        if rep.transform_name == "centered_log_ratio":
            clr_col_sum = float(np.abs(rep.transformed.sum(axis=0)).max())
        rows.append(
            {
                "analysis_id": rep.analysis_id,
                "analysis_role": rep.role,
                "n_features": rep.filtered.shape[0],
                "n_samples": rep.filtered.shape[1],
                "prevalence_threshold": rep.prevalence_threshold,
                "abundance_threshold": rep.abundance_threshold,
                "pseudocount": rep.pseudocount,
                "transformation": rep.transform_name,
                "distance_metric": rep.distance_metric,
                "contaminant_policy": rep.contaminant_policy,
                "sample_exclusion_policy": rep.sample_exclusion_policy,
                "zero_fraction_before_replacement": float((values == 0).mean()),
                "clr_min": float(np.nanmin(transformed)),
                "clr_max": float(np.nanmax(transformed)),
                "max_abs_clr_column_sum": clr_col_sum,
                "random_seed": RANDOM_SEED,
                "original_audited_sha256": original_sha,
                "python_version": sys.version.split()[0],
                "platform": platform.platform(),
                "package_versions_json": json.dumps(package_versions, sort_keys=True),
                "output_files": ";".join(rep.output_files),
            }
        )
    return pd.DataFrame(rows)


def contamination_category(genus: str) -> tuple[str, str, str]:
    if genus in HIGH_RISK:
        return ("HIGH_RISK_POTENTIAL_CONTAMINANT", "Phase 6B locked environmental/reagent flag list", "Potential contaminant flag only; not confirmed contamination.")
    if genus in MODERATE_RISK:
        return ("MODERATE_RISK_ENVIRONMENTAL", "Phase 6B locked environmental/reagent flag list", "Environmental risk flag only; not confirmed contamination.")
    if genus in BIO_PLAUSIBLE:
        return ("BIOLOGICALLY_PLAUSIBLE_BUT_CONTAMINATION_SENSITIVE", "Phase 6B locked contamination-sensitivity list", "Retained in primary; interpret with sensitivity analyses.")
    return ("LOW_CURRENT_CONCERN", "No Phase 6B contamination flag", "No current contamination-specific flag.")


def write_qc_tables(audited: pd.DataFrame, primary: Representation, reps: list[Representation], sample_qc6a: pd.DataFrame) -> None:
    prevalence = (audited > 0).sum(axis=1) / audited.shape[1]
    taxa_rows = []
    no_high = set(reps_by_id(reps)["MICRO_SENS_NO_HIGH_RISK"].filtered.index)
    no_high_mod = set(reps_by_id(reps)["MICRO_SENS_NO_CONTAMINANTS"].filtered.index)
    for genus in audited.index:
        cat, source, note = contamination_category(genus)
        primary_retained = genus in primary.filtered.index
        sens_members = []
        if genus in no_high:
            sens_members.append("MICRO_SENS_NO_HIGH_RISK")
        if genus in no_high_mod:
            sens_members.append("MICRO_SENS_NO_CONTAMINANTS")
        taxa_rows.append(
            {
                "genus": genus,
                "primary_retained": primary_retained,
                "prevalence": prevalence.loc[genus],
                "mean_abundance": audited.loc[genus].mean(),
                "median_abundance": audited.loc[genus].median(),
                "zero_fraction": float((audited.loc[genus] == 0).mean()),
                "contamination_risk_category": cat,
                "evidence_source": source,
                "included_in_primary": primary_retained,
                "included_in_sensitivity": ";".join(sens_members) if sens_members else "not_in_contaminant_removal_sensitivity_or_filtered",
                "notes": note,
            }
        )
    taxa_flags = pd.DataFrame(taxa_rows)
    taxa_flags.loc[taxa_flags["primary_retained"]].to_csv(TABLE_DIR / "phase6c_retained_taxa_with_contamination_flags.tsv", sep="\t", index=False)

    sample_rows = []
    for sample in primary.filtered.columns:
        vals = primary.filtered[sample]
        sample_rows.append(
            {
                "sample": sample,
                "patient_id": sample_qc6a.loc[sample, "patient_id"] if sample in sample_qc6a.index else "",
                "taxa_detected_per_sample": int((vals > 0).sum()),
                "sample_total_abundance_proxy": float(audited[sample].sum()),
                "filtered_total_abundance_proxy": float(vals.sum()),
                "zero_fraction_before_replacement": float((vals == 0).mean()),
                "pseudocount_used": PRIMARY_PSEUDOCOUNT,
                "clr_min": float(primary.transformed[sample].min()),
                "clr_max": float(primary.transformed[sample].max()),
                "technical_outlier_status": "RETAIN_WITH_SENSITIVITY_ANALYSIS" if sample in EXTREME_SAMPLES else "RETAIN",
                "degenerate_after_filtering": bool((vals > 0).sum() == 0 or vals.sum() == 0),
            }
        )
    pd.DataFrame(sample_rows).to_csv(TABLE_DIR / "phase6c_processed_sample_qc.tsv", sep="\t", index=False)

    taxon_rows = []
    for genus in primary.filtered.index:
        vals = primary.filtered.loc[genus]
        taxon_rows.append(
            {
                "genus": genus,
                "prevalence_count": int((vals > 0).sum()),
                "prevalence_fraction": float((vals > 0).mean()),
                "mean_abundance": float(vals.mean()),
                "median_abundance": float(vals.median()),
                "zero_fraction_before_replacement": float((vals == 0).mean()),
                "pseudocount_used": PRIMARY_PSEUDOCOUNT,
                "clr_min": float(primary.transformed.loc[genus].min()),
                "clr_max": float(primary.transformed.loc[genus].max()),
                "contamination_risk_category": contamination_category(genus)[0],
            }
        )
    pd.DataFrame(taxon_rows).to_csv(TABLE_DIR / "phase6c_processed_taxon_qc.tsv", sep="\t", index=False)


def reps_by_id(reps: list[Representation]) -> dict[str, Representation]:
    return {rep.analysis_id: rep for rep in reps}


def upper_triangle(distance: pd.DataFrame, samples: list[str]) -> np.ndarray:
    sub = distance.loc[samples, samples].to_numpy(dtype=float)
    return sub[np.triu_indices_from(sub, k=1)]


def concordance_table(reps: list[Representation], primary: Representation) -> pd.DataFrame:
    rows = []
    primary_mean = primary.filtered.mean(axis=1).sort_values(ascending=False)
    primary_coords = pcoa(primary.distance)
    for rep in reps:
        if rep.analysis_id == "MICRO_PRIMARY":
            continue
        shared_samples = [s for s in primary.distance.index if s in rep.distance.index]
        shared_taxa = [g for g in primary.filtered.index if g in rep.filtered.index]
        x = upper_triangle(primary.distance, shared_samples)
        y = upper_triangle(rep.distance, shared_samples)
        mantel = spearmanr(x, y).statistic if len(x) > 2 else np.nan

        procrustes_corr = np.nan
        sample_order = np.nan
        if len(shared_samples) >= 3:
            pc1 = primary_coords.loc[shared_samples]
            pc2 = pcoa(rep.distance.loc[shared_samples, shared_samples]).loc[shared_samples]
            _, _, disparity = procrustes(pc1.to_numpy(), pc2.to_numpy())
            procrustes_corr = 1.0 - disparity
            sample_order = abs(spearmanr(pc1["PCoA1"], pc2["PCoA1"]).statistic)

        taxon_rank = np.nan
        if len(shared_taxa) >= 3:
            taxon_rank = spearmanr(primary_mean.loc[shared_taxa], rep.filtered.loc[shared_taxa].mean(axis=1)).statistic

        clr_corr = np.nan
        if rep.transform_name in {"centered_log_ratio", "robust_clr"} and len(shared_taxa) > 0:
            shared = primary.transformed.loc[shared_taxa, shared_samples].to_numpy().ravel()
            other = rep.transformed.loc[shared_taxa, shared_samples].to_numpy().ravel()
            clr_corr = spearmanr(shared, other).statistic

        ord_shift = np.nan
        if rep.analysis_id in {"MICRO_SENS_NO_CONTAMINANTS", "MICRO_SENS_NO_HIGH_RISK"} and len(shared_samples) >= 3:
            pc1 = primary_coords.loc[shared_samples]
            pc2 = pcoa(rep.distance.loc[shared_samples, shared_samples]).loc[shared_samples]
            aligned1, aligned2, _ = procrustes(pc1.to_numpy(), pc2.to_numpy())
            ord_shift = float(np.linalg.norm(aligned1 - aligned2, axis=1).mean())

        rows.append(
            {
                "comparison": f"MICRO_PRIMARY_vs_{rep.analysis_id}",
                "analysis_id": rep.analysis_id,
                "shared_samples": len(shared_samples),
                "shared_taxa": len(shared_taxa),
                "mantel_spearman_r": mantel,
                "procrustes_concordance_1_minus_disparity": procrustes_corr,
                "sample_order_stability_abs_spearman_pc1": sample_order,
                "taxon_rank_stability_spearman_mean_abundance": taxon_rank,
                "clr_correlation_spearman": clr_corr,
                "mean_ordination_shift_after_contaminant_removal": ord_shift,
                "biological_outcomes_used": False,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "phase6c_preprocessing_sensitivity_concordance.tsv", sep="\t", index=False)
    return out


def make_figures(audited: pd.DataFrame, primary: Representation, reps: list[Representation], concordance: pd.DataFrame, sample_qc: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")

    plt.figure(figsize=(10, 8))
    order = primary.transformed.var(axis=1).sort_values(ascending=False).head(60).index
    sns.heatmap(primary.transformed.loc[order], cmap="vlag", center=0, xticklabels=False, yticklabels=True, cbar_kws={"label": "CLR"})
    plt.xlabel("Tumor samples")
    plt.ylabel("Top variable retained genera")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase6c_primary_CLR_heatmap.pdf")
    plt.close()

    coords = pcoa(primary.distance)
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = ["#b33f62" if s in EXTREME_SAMPLES else "#35618f" for s in coords.index]
    ax.scatter(coords["PCoA1"], coords["PCoA2"], c=colors, s=36, edgecolor="black", linewidth=0.4)
    ax.set_xlabel("PCoA1")
    ax.set_ylabel("PCoA2")
    ax.set_title("Primary Aitchison PCoA")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase6c_primary_Aitchison_PCoA.pdf")
    plt.close()

    inv_rows = []
    for rep in reps:
        inv_rows.append({"analysis_id": rep.analysis_id, "features": rep.filtered.shape[0], "zero_fraction": (rep.filtered.to_numpy() == 0).mean()})
    inv = pd.DataFrame(inv_rows)
    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax1.bar(inv["analysis_id"], inv["features"], color="#477c6e")
    ax1.set_ylabel("Genera retained")
    ax1.tick_params(axis="x", rotation=60)
    ax2 = ax1.twinx()
    ax2.plot(inv["analysis_id"], inv["zero_fraction"], color="#8c3b3b", marker="o")
    ax2.set_ylabel("Zero fraction")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase6c_filtering_sensitivity_summary.pdf")
    plt.close()

    pseudo = concordance[concordance["analysis_id"].isin(["MICRO_SENS_PSEUDO_0.1", "MICRO_SENS_PSEUDO_1.0", "MICRO_SENS_ROBUST_CLR"])]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(pseudo["analysis_id"], pseudo["clr_correlation_spearman"], color="#5f6f9c")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Spearman correlation vs primary CLR")
    ax.tick_params(axis="x", rotation=30)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase6c_pseudocount_sensitivity.pdf")
    plt.close()

    cats = pd.Series([contamination_category(g)[0] for g in primary.filtered.index]).value_counts().rename_axis("category").reset_index(name="n")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(cats["category"], cats["n"], color="#7a5c8d")
    ax.set_xlabel("Primary retained genera")
    ax.set_ylabel("Flag category")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase6c_contamination_flag_summary.pdf")
    plt.close()

    sq = sample_qc.reset_index().rename(columns={"index": "sample"})
    fig, ax = plt.subplots(figsize=(8, 4))
    bar_colors = ["#b33f62" if s in EXTREME_SAMPLES else "#3f7b7b" for s in sq["microbiome_matrix_sample"]]
    ax.bar(np.arange(len(sq)), sq["total_abundance"], color=bar_colors)
    ax.set_yscale("log")
    ax.set_xlabel("Tumor samples in locked order")
    ax.set_ylabel("Matrix total-abundance proxy")
    ax.set_xticks([])
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase6c_sample_depth_proxy.pdf")
    plt.close()


def package_versions() -> dict[str, str]:
    versions = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "seaborn": sns.__version__,
    }
    try:
        import scipy

        versions["scipy"] = scipy.__version__
    except Exception as exc:  # pragma: no cover
        versions["scipy"] = f"unavailable: {exc}"
    return versions


def main() -> None:
    np.random.seed(RANDOM_SEED)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SENS_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    warning_records: list[dict[str, str]] = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        params = locked_file_consistency()
        original_sha = sha256(AUDITED)
        audited = pd.read_csv(AUDITED, sep="\t", index_col=0)
        crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
        sample_qc = pd.read_csv(PHASE6A_SAMPLE_QC, sep="\t", dtype={"tumor_number": str}).set_index("microbiome_matrix_sample")

        require(audited.shape == (365, 62), f"Unexpected audited matrix shape: {audited.shape}")
        require(audited.index.is_unique, "Duplicated genera in audited matrix")
        require(audited.columns.is_unique, "Duplicated samples in audited matrix")
        require(list(crosswalk["microbiome_matrix_sample"]) == list(audited.columns), "Unexpected sample ordering")
        require(np.isfinite(audited.to_numpy(dtype=float)).all(), "Audited matrix has missing or infinite values")
        require((audited.to_numpy(dtype=float) >= 0).all(), "Audited matrix has negative values")
        min_nonzero = float(audited.to_numpy()[audited.to_numpy() > 0].min())
        require(math.isclose(min_nonzero, MIN_NONZERO_LOCKED, rel_tol=1e-12), f"Minimum non-zero value {min_nonzero} differs from lock")

        primary = make_rep("MICRO_PRIMARY", audited, 0.2, 0.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "keep_all_flagged", "retain_all", role="primary")
        if primary.filtered.shape[0] != EXPECTED_PRIMARY_FEATURES:
            raise AssertionError(f"Primary retained {primary.filtered.shape[0]} genera; Phase 6B expected approximately {EXPECTED_PRIMARY_FEATURES}. Stopping without forcing count.")

        no_high_source = audited.drop(index=[g for g in HIGH_RISK if g in audited.index])
        no_high_mod_source = audited.drop(index=[g for g in (HIGH_RISK | MODERATE_RISK) if g in audited.index])
        no_extreme_source = audited.drop(columns=EXTREME_SAMPLES)

        reps = [
            primary,
            make_rep("MICRO_SENS_PREV_10", audited, 0.1, 0.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "keep_all_flagged", "retain_all"),
            make_rep("MICRO_SENS_PREV_30", audited, 0.3, 0.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "keep_all_flagged", "retain_all"),
            make_rep("MICRO_SENS_DET_10_P20", audited, 0.2, 10.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "keep_all_flagged", "retain_all"),
            make_rep("MICRO_SENS_PSEUDO_0.1", audited, 0.2, 0.0, 0.1, "centered_log_ratio", "aitchison", "keep_all_flagged", "retain_all"),
            make_rep("MICRO_SENS_PSEUDO_1.0", audited, 0.2, 0.0, 1.0, "centered_log_ratio", "aitchison", "keep_all_flagged", "retain_all"),
            make_rep("MICRO_SENS_ROBUST_CLR", audited, 0.2, 0.0, None, "robust_clr", "aitchison", "keep_all_flagged", "retain_all"),
            make_rep("MICRO_SENS_NO_HIGH_RISK", no_high_source, 0.2, 0.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "exclude_high_risk", "retain_all"),
            make_rep("MICRO_SENS_NO_CONTAMINANTS", no_high_mod_source, 0.2, 0.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "exclude_high_and_moderate_risk", "retain_all"),
            make_rep("MICRO_SENS_EXCLUDE_EXTREME", no_extreme_source, 0.2, 0.0, PRIMARY_PSEUDOCOUNT, "centered_log_ratio", "aitchison", "keep_all_flagged", "exclude_extreme"),
            make_rep("MICRO_SENS_PRESENCE_ABSENCE", audited, 0.2, 0.0, None, "presence_absence", "jaccard", "keep_all_flagged", "retain_all"),
        ]
        for rep in reps:
            save_representation(rep)

        write_qc_tables(audited, primary, reps, sample_qc)
        concordance = concordance_table(reps, primary)
        inventory = matrix_inventory(reps, package_versions(), original_sha)
        inventory.to_csv(TABLE_DIR / "phase6c_matrix_inventory.tsv", sep="\t", index=False)
        make_figures(audited, primary, reps, concordance, sample_qc)

        require(sha256(AUDITED) == original_sha, "Original audited matrix was modified")
        subprocess.run([sys.executable, str(SUMMARY_SCRIPT)], cwd=ROOT, check=True)

        for warning in caught:
            warning_records.append(
                {
                    "warning_category": warning.category.__name__,
                    "warning_message": str(warning.message),
                    "source_file": str(warning.filename),
                    "line_number": warning.lineno,
                }
            )

    if not warning_records:
        warning_records.append({"warning_category": "NONE", "warning_message": "No warnings captured", "source_file": "", "line_number": ""})
    pd.DataFrame(warning_records).to_csv(WARNINGS_PATH, sep="\t", index=False)

    print("Phase 6C preprocessing complete")
    print(f"Primary retained genera: {primary.filtered.shape[0]}")
    print(f"Primary samples: {primary.filtered.shape[1]}")
    print(f"Primary CLR max abs sample sum: {np.abs(primary.transformed.sum(axis=0)).max():.3e}")
    print(f"Primary Aitchison distance max: {primary.distance.to_numpy().max():.6f}")


if __name__ == "__main__":
    main()
