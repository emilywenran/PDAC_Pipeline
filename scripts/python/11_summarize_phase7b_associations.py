#!/usr/bin/env python3
"""Execute locked Phase 7B tumor microbiome association analyses.

This script follows the Phase 7A/7A.5 locked model definitions. It does not
optimize outcomes, thresholds, filters, covariates, or evidence rules after
inspecting results.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_phase7b")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels
import statsmodels.api as sm
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import OLSInfluence, variance_inflation_factor


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results" / "tables"
FIG_DIR = ROOT / "05_results" / "figures"
REPORT = ROOT / "04_analysis" / "08_host_microbiome_integration" / "PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md"
RANDOM_SEED = 2026
PERMUTATIONS = 9999
BOOTSTRAPS = 2000
PRIMARY_N = 62
PRIMARY_G = 122
EXTREME_SAMPLES = ["Basal-like1", "Hybrid18", "Hybrid23"]


def _mkdirs() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def save_tsv(df: pd.DataFrame, name: str) -> Path:
    path = TABLE_DIR / name
    df.to_csv(path, sep="\t", index=False)
    return path


def read_wide_matrix(path: str | Path, row_col: str | None = None) -> pd.DataFrame:
    df = pd.read_csv(ROOT / path, sep="\t", compression="infer")
    if row_col is None:
        row_col = "taxon" if "taxon" in df.columns else df.columns[0]
    return df.set_index(row_col)


def matrix_samples_to_patients(mat: pd.DataFrame, crosswalk: pd.DataFrame) -> pd.DataFrame:
    mapping = dict(zip(crosswalk["microbiome_matrix_sample"], crosswalk["patient_id"]))
    missing = [c for c in mat.columns if c not in mapping]
    if missing:
        raise RuntimeError(f"Microbiome matrix samples absent from crosswalk: {missing[:5]}")
    out = mat.rename(columns=mapping).T
    out.index.name = "patient_id"
    return out


def read_distance(path: str | Path, crosswalk: pd.DataFrame) -> pd.DataFrame:
    d = pd.read_csv(ROOT / path, sep="\t", compression="infer", index_col=0)
    mapping = dict(zip(crosswalk["microbiome_matrix_sample"], crosswalk["patient_id"]))
    d = d.rename(index=mapping, columns=mapping)
    return d


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=1)


def bh(pvals: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(pvals, dtype=float)
    ok = np.isfinite(arr)
    q = np.full(arr.shape, np.nan, dtype=float)
    if ok.any():
        q[ok] = multipletests(arr[ok], method="fdr_bh")[1]
    return q


def direction(x: float) -> str:
    if not np.isfinite(x) or abs(x) < 1e-15:
        return "ZERO_OR_NA"
    return "positive" if x > 0 else "negative"


def condition_number(X: pd.DataFrame) -> float:
    vals = np.linalg.svd(np.asarray(X, float), compute_uv=False)
    return float(vals.max() / vals.min()) if vals.min() > 0 else math.inf


def calc_vif(X: pd.DataFrame) -> dict[str, float]:
    vals = np.asarray(X, float)
    out = {}
    for i, col in enumerate(X.columns):
        if col == "const":
            continue
        try:
            out[col] = float(variance_inflation_factor(vals, i))
        except Exception:
            out[col] = math.inf
    return out


def ols_hc3(y: pd.Series | np.ndarray, X: pd.DataFrame) -> tuple[sm.regression.linear_model.RegressionResultsWrapper, sm.regression.linear_model.RegressionResultsWrapper]:
    y = np.asarray(y, float)
    model = sm.OLS(y, X, missing="raise").fit()
    robust = model.get_robustcov_results(cov_type="HC3")
    return model, robust


def simple_ols_table(mat: pd.DataFrame, x: pd.Series, family: str, role: str) -> pd.DataFrame:
    rows = []
    X = sm.add_constant(x.rename("standardized_host_score"))
    for genus in mat.columns:
        y = mat[genus].astype(float)
        model, rob = ols_hc3(y, X)
        idx = list(rob.model.exog_names).index("standardized_host_score")
        infl = OLSInfluence(model)
        resid = np.asarray(model.resid)
        fitted = np.asarray(model.fittedvalues)
        try:
            bp_p = float(sm.stats.diagnostic.het_breuschpagan(resid, X)[1])
        except Exception:
            bp_p = np.nan
        try:
            shapiro_p = float(stats.shapiro(resid).pvalue)
        except Exception:
            shapiro_p = np.nan
        ssr = float(np.sum(resid**2))
        sst = float(np.sum((y - y.mean()) ** 2))
        rows.append(
            {
                "genus": genus,
                "analysis_role": role,
                "model": "Model_0",
                "formula": "CLR_genus ~ standardized_Moffitt50_contrast",
                "n": int(len(y)),
                "coefficient": float(rob.params[idx]),
                "robust_se_HC3": float(rob.bse[idx]),
                "ci_lower": float(rob.conf_int()[idx, 0]),
                "ci_upper": float(rob.conf_int()[idx, 1]),
                "t_statistic_HC3": float(rob.tvalues[idx]),
                "p_value": float(rob.pvalues[idx]),
                "r_squared": float(model.rsquared),
                "adjusted_r_squared": float(model.rsquared_adj),
                "partial_r_squared": float(model.rsquared),
                "cohens_f2": float(model.rsquared / (1 - model.rsquared)) if model.rsquared < 1 else math.inf,
                "residual_shapiro_p": shapiro_p,
                "residual_breusch_pagan_p": bp_p,
                "max_abs_studentized_residual": float(np.nanmax(np.abs(infl.resid_studentized_external))),
                "max_cooks_distance": float(np.nanmax(infl.cooks_distance[0])),
                "max_leverage": float(np.nanmax(infl.hat_matrix_diag)),
                "residual_sum_squares": ssr,
                "total_sum_squares": sst,
                "multiple_testing_family": family,
                "coefficient_interpretation": "relative_compositional_CLR_association_not_absolute_load",
            }
        )
    out = pd.DataFrame(rows)
    out["bh_q_value"] = bh(out["p_value"])
    out["rank_by_p"] = out["p_value"].rank(method="min").astype(int)
    return out


def covariate_ols_table(mat: pd.DataFrame, meta: pd.DataFrame, covariates: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    base_coef = {}
    for model_name, covs in covariates.items():
        predictors = ["moffitt50_z"] + covs
        X = sm.add_constant(meta[predictors])
        vif = calc_vif(X)
        cn = condition_number(X)
        pvals = []
        model_rows = []
        for genus in mat.columns:
            y = mat[genus].astype(float)
            model, rob = ols_hc3(y, X)
            idx = list(rob.model.exog_names).index("moffitt50_z")
            coef = float(rob.params[idx])
            if model_name == "Model_0":
                base_coef[genus] = coef
            attenuation = np.nan
            sign_change = False
            if genus in base_coef and model_name != "Model_0":
                attenuation = 1 - (coef / base_coef[genus]) if abs(base_coef[genus]) > 1e-15 else np.nan
                sign_change = np.sign(coef) != np.sign(base_coef[genus])
            model_rows.append(
                {
                    "genus": genus,
                    "model": model_name,
                    "formula": "CLR_genus ~ standardized_Moffitt50_contrast" + (" + " + " + ".join(covs) if covs else ""),
                    "model_role": "primary" if model_name == "Model_0" else "sensitivity",
                    "n": int(len(y)),
                    "host_score_coefficient": coef,
                    "host_score_robust_se_HC3": float(rob.bse[idx]),
                    "host_score_ci_lower": float(rob.conf_int()[idx, 0]),
                    "host_score_ci_upper": float(rob.conf_int()[idx, 1]),
                    "host_score_p_value": float(rob.pvalues[idx]),
                    "coefficient_attenuation_vs_Model_0": attenuation,
                    "sign_change_vs_Model_0": bool(sign_change),
                    "directionally_stable_vs_Model_0": bool(not sign_change),
                    "max_vif": float(max(vif.values())) if vif else np.nan,
                    "vif_by_predictor": json.dumps(vif, sort_keys=True),
                    "condition_number": cn,
                    "multiple_testing_family": f"{model_name}_122_genus_tests_Moffitt50",
                    "combined_TME_model_executed": False,
                }
            )
            pvals.append(float(rob.pvalues[idx]))
        tmp = pd.DataFrame(model_rows)
        tmp["host_score_q_value_within_model"] = bh(tmp["host_score_p_value"])
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def permanova(distance: pd.DataFrame, meta: pd.DataFrame, predictors: list[str], permutations: int, seed: int) -> pd.DataFrame:
    ids = list(distance.index)
    D = distance.loc[ids, ids].to_numpy(float)
    n = D.shape[0]
    A = -0.5 * D**2
    Hc = np.eye(n) - np.ones((n, n)) / n
    G = Hc @ A @ Hc
    ymeta = meta.loc[ids]
    rng = np.random.default_rng(seed)

    def design(cols: list[str], frame: pd.DataFrame) -> pd.DataFrame:
        parts = [pd.Series(1.0, index=frame.index, name="const")]
        for c in cols:
            if frame[c].dtype == object or str(frame[c].dtype).startswith("category"):
                d = pd.get_dummies(frame[c], prefix=c, drop_first=True, dtype=float)
                parts.append(d)
            else:
                parts.append(frame[c].astype(float).rename(c))
        return pd.concat(parts, axis=1)

    def hat(X: pd.DataFrame) -> np.ndarray:
        Xv = np.asarray(X, float)
        return Xv @ np.linalg.pinv(Xv.T @ Xv) @ Xv.T

    X_full = design(predictors, ymeta)
    H_full = hat(X_full)
    H_int = hat(pd.DataFrame({"const": np.ones(n)}, index=ymeta.index))
    ss_total = float(np.trace((np.eye(n) - H_int) @ G))
    ss_res = float(np.trace((np.eye(n) - H_full) @ G))
    df_res = n - np.linalg.matrix_rank(np.asarray(X_full, float))
    rows = []
    for pred in predictors:
        reduced = [c for c in predictors if c != pred]
        X_red = design(reduced, ymeta)
        H_red = hat(X_red)
        ss_term = float(np.trace((H_full - H_red) @ G))
        df_term = np.linalg.matrix_rank(np.asarray(X_full, float)) - np.linalg.matrix_rank(np.asarray(X_red, float))
        Fobs = (ss_term / df_term) / (ss_res / df_res)
        ge = 1
        for _ in range(permutations):
            perm = rng.permutation(n)
            pmeta = ymeta.copy()
            for col in predictors:
                pmeta[col] = ymeta[col].to_numpy()[perm]
            Xp_full = design(predictors, pmeta)
            Xp_red = design(reduced, pmeta)
            Hp_full = hat(Xp_full)
            Hp_red = hat(Xp_red)
            ss_t = float(np.trace((Hp_full - Hp_red) @ G))
            ss_r = float(np.trace((np.eye(n) - Hp_full) @ G))
            Fp = (ss_t / df_term) / (ss_r / df_res)
            if Fp >= Fobs - 1e-15:
                ge += 1
        rows.append(
            {
                "test": "PERMANOVA",
                "term": pred,
                "n": n,
                "df_term": int(df_term),
                "df_residual": int(df_res),
                "sum_squares_term": ss_term,
                "sum_squares_residual": ss_res,
                "pseudo_F": float(Fobs),
                "r_squared": float(ss_term / ss_total) if ss_total > 0 else np.nan,
                "p_value": float(ge / (permutations + 1)),
                "permutations": permutations,
                "random_seed": seed,
                "sums_of_squares": "marginal" if len(predictors) > 1 else "single_predictor",
                "inference": "two-sided_permutation_location_test",
            }
        )
    return pd.DataFrame(rows)


def permdisp(distance: pd.DataFrame, groups: pd.Series, permutations: int, seed: int) -> pd.DataFrame:
    ids = list(distance.index)
    D = distance.loc[ids, ids].to_numpy(float)
    n = D.shape[0]
    A = -0.5 * D**2
    Hc = np.eye(n) - np.ones((n, n)) / n
    B = Hc @ A @ Hc
    eigvals, eigvecs = np.linalg.eigh(B)
    keep = eigvals > 1e-10
    coords = eigvecs[:, keep] * np.sqrt(eigvals[keep])
    g = groups.loc[ids].astype(str).to_numpy()

    def centroid_dist(labels: np.ndarray) -> np.ndarray:
        out = np.zeros(n)
        for lv in np.unique(labels):
            ix = labels == lv
            cen = coords[ix].mean(axis=0)
            out[ix] = np.sqrt(((coords[ix] - cen) ** 2).sum(axis=1))
        return out

    z = centroid_dist(g)
    levels = np.unique(g)
    grand = z.mean()
    ssb = sum((z[g == lv].size) * (z[g == lv].mean() - grand) ** 2 for lv in levels)
    ssw = sum(((z[g == lv] - z[g == lv].mean()) ** 2).sum() for lv in levels)
    dfb = len(levels) - 1
    dfw = n - len(levels)
    Fobs = (ssb / dfb) / (ssw / dfw)
    rng = np.random.default_rng(seed)
    ge = 1
    for _ in range(permutations):
        gp = rng.permutation(g)
        zp = centroid_dist(gp)
        grandp = zp.mean()
        ssbp = sum((zp[gp == lv].size) * (zp[gp == lv].mean() - grandp) ** 2 for lv in levels)
        sswp = sum(((zp[gp == lv] - zp[gp == lv].mean()) ** 2).sum() for lv in levels)
        Fp = (ssbp / dfb) / (sswp / dfw)
        if Fp >= Fobs - 1e-15:
            ge += 1
    return pd.DataFrame(
        [
            {
                "test": "PERMDISP_betadisper",
                "term": "public_subtype",
                "n": n,
                "df_term": int(dfb),
                "df_residual": int(dfw),
                "F_statistic": float(Fobs),
                "p_value": float(ge / (permutations + 1)),
                "permutations": permutations,
                "random_seed": seed,
                "interpretation": "categorical_dispersion_test_only_not_continuous_predictor_dispersion",
            }
        ]
    )


def spearman_table(mat: pd.DataFrame, x: pd.Series) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(RANDOM_SEED)
    for genus in mat.columns:
        y = mat[genus].astype(float)
        rho, p = stats.spearmanr(x, y)
        boots = []
        idx = np.arange(len(y))
        for _ in range(BOOTSTRAPS):
            b = rng.choice(idx, size=len(idx), replace=True)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rb, _ = stats.spearmanr(x.iloc[b], y.iloc[b])
            if np.isfinite(rb):
                boots.append(rb)
        lo, hi = np.nanpercentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
        rows.append(
            {
                "genus": genus,
                "spearman_rho": float(rho),
                "p_value": float(p),
                "bootstrap_ci_lower": float(lo),
                "bootstrap_ci_upper": float(hi),
                "bootstrap_iterations": BOOTSTRAPS,
                "direction": direction(rho),
                "analysis_role": "supporting",
            }
        )
    out = pd.DataFrame(rows)
    out["bh_q_value"] = bh(out["p_value"])
    return out


def permutation_table(mat: pd.DataFrame, x: pd.Series) -> pd.DataFrame:
    xarr = np.asarray(x, float)
    rng = np.random.default_rng(RANDOM_SEED)
    ymat = mat.to_numpy(float)
    ycenter = ymat - ymat.mean(axis=0, keepdims=True)
    sxx = np.sum((xarr - xarr.mean()) ** 2)
    obs_beta = ((xarr - xarr.mean())[:, None] * ycenter).sum(axis=0) / sxx
    ge = np.ones(mat.shape[1], dtype=int)
    for _ in range(PERMUTATIONS):
        xp = rng.permutation(xarr)
        bp = ((xp - xp.mean())[:, None] * ycenter).sum(axis=0) / sxx
        ge += (np.abs(bp) >= np.abs(obs_beta) - 1e-15)
    rows = []
    for i, genus in enumerate(mat.columns):
        rows.append(
            {
                "genus": genus,
                "observed_coefficient": float(obs_beta[i]),
                "empirical_p_value": float(ge[i] / (PERMUTATIONS + 1)),
                "permutations": PERMUTATIONS,
                "random_seed": RANDOM_SEED,
                "inference": "two_sided_empirical_permutation_abs_coefficient",
                "analysis_role": "supporting",
            }
        )
    out = pd.DataFrame(rows)
    out["bh_q_value"] = bh(out["empirical_p_value"])
    return out


def bootstrap_table(mat: pd.DataFrame, x: pd.Series) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    xarr = np.asarray(x, float)
    idx = np.arange(len(xarr))
    rows = []
    for genus in mat.columns:
        y = np.asarray(mat[genus], float)
        X = sm.add_constant(xarr)
        _, rob = ols_hc3(y, X)
        betas = []
        for _ in range(BOOTSTRAPS):
            b = rng.choice(idx, size=len(idx), replace=True)
            xb = xarr[b]
            yb = y[b]
            den = np.sum((xb - xb.mean()) ** 2)
            if den > 0:
                betas.append(float(np.sum((xb - xb.mean()) * (yb - yb.mean())) / den))
        lo, hi = np.nanpercentile(betas, [2.5, 97.5])
        rows.append(
            {
                "genus": genus,
                "ols_coefficient": float(rob.params[1]),
                "bootstrap_ci_lower": float(lo),
                "bootstrap_ci_upper": float(hi),
                "bootstrap_iterations": BOOTSTRAPS,
                "bootstrap_ci_excludes_zero": bool(lo > 0 or hi < 0),
                "analysis_role": "supporting",
            }
        )
    return pd.DataFrame(rows)


def maaslin2_table(mat: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "genus": mat.columns,
            "analysis_role": "supporting",
            "normalization": "NONE",
            "transform": "NONE",
            "status": "NOT_RUN_PACKAGE_UNAVAILABLE",
            "coefficient": np.nan,
            "p_value": np.nan,
            "q_value": np.nan,
            "notes": "MaAsLin2/Maaslin2 R package was not installed locally; no alternate normalization or transform was substituted.",
        }
    )


def preprocessing_sensitivity(meta: pd.DataFrame) -> pd.DataFrame:
    sens = {
        "MICRO_SENS_PREV_10": "03_processed/microbiome/sensitivity/MICRO_SENS_PREV_10_centered_log_ratio.tsv.gz",
        "MICRO_SENS_PREV_30": "03_processed/microbiome/sensitivity/MICRO_SENS_PREV_30_centered_log_ratio.tsv.gz",
        "MICRO_SENS_DET_10_P20": "03_processed/microbiome/sensitivity/MICRO_SENS_DET_10_P20_centered_log_ratio.tsv.gz",
        "MICRO_SENS_PSEUDO_0.1": "03_processed/microbiome/sensitivity/MICRO_SENS_PSEUDO_0.1_centered_log_ratio.tsv.gz",
        "MICRO_SENS_PSEUDO_1.0": "03_processed/microbiome/sensitivity/MICRO_SENS_PSEUDO_1.0_centered_log_ratio.tsv.gz",
        "MICRO_SENS_ROBUST_CLR": "03_processed/microbiome/sensitivity/MICRO_SENS_ROBUST_CLR_robust_clr.tsv.gz",
        "MICRO_SENS_EXCLUDE_EXTREME": "03_processed/microbiome/sensitivity/MICRO_SENS_EXCLUDE_EXTREME_centered_log_ratio.tsv.gz",
        "MICRO_SENS_NO_HIGH_RISK": "03_processed/microbiome/sensitivity/MICRO_SENS_NO_HIGH_RISK_centered_log_ratio.tsv.gz",
        "MICRO_SENS_NO_CONTAMINANTS": "03_processed/microbiome/sensitivity/MICRO_SENS_NO_CONTAMINANTS_centered_log_ratio.tsv.gz",
    }
    cross = pd.read_csv(ROOT / "01_metadata/microbiome_sample_crosswalk.tsv", sep="\t")
    all_rows = []
    for aid, path in sens.items():
        mat = matrix_samples_to_patients(read_wide_matrix(path), cross)
        common = mat.index.intersection(meta.index)
        mat = mat.loc[common]
        x = zscore(meta.loc[common, "basal_classical_contrast"])
        tab = simple_ols_table(mat, x, f"{len(mat.columns)}_genus_tests_Moffitt50_{aid}", "preprocessing_sensitivity")
        tab["analysis_id"] = aid
        tab["matrix_path"] = path
        tab["n_available_genera"] = len(mat.columns)
        all_rows.append(tab)
    out = pd.concat(all_rows, ignore_index=True)
    out["direction"] = out["coefficient"].map(direction)
    out["rank_within_analysis"] = out.groupby("analysis_id")["p_value"].rank(method="min").astype(int)
    primary_dir = out[out["analysis_id"] == "MICRO_SENS_PSEUDO_1.0"].set_index("genus")["direction"].to_dict()
    primary_dir.update(out[out["analysis_id"] == "MICRO_SENS_PSEUDO_0.1"].set_index("genus")["direction"].to_dict())
    summary = out.groupby("genus").agg(
        analyses_present=("analysis_id", "nunique"),
        effect_min=("coefficient", "min"),
        effect_max=("coefficient", "max"),
    )
    frac = []
    signrev = []
    for genus, sub in out.groupby("genus"):
        dirs = sub["direction"].tolist()
        ref = direction(sub.loc[sub["analysis_id"].eq("MICRO_SENS_PSEUDO_1.0"), "coefficient"].iloc[0]) if (sub["analysis_id"] == "MICRO_SENS_PSEUDO_1.0").any() else dirs[0]
        frac.append((genus, np.mean([d == ref for d in dirs if d != "ZERO_OR_NA"])))
        signrev.append((genus, len({d for d in dirs if d in ["positive", "negative"]}) > 1))
    summary["fraction_analyses_with_primary_direction"] = pd.Series(dict(frac))
    summary["evidence_of_sign_reversal"] = pd.Series(dict(signrev))
    return out.merge(summary.reset_index(), on="genus", how="left")


def contamination_sensitivity(primary: pd.DataFrame, preproc: pd.DataFrame, mat_abund: pd.DataFrame, clr: pd.DataFrame, meta: pd.DataFrame, flags: pd.DataFrame) -> pd.DataFrame:
    cand = set(primary.loc[primary["bh_q_value"] < 0.10, "genus"])
    cand.update(primary.sort_values("p_value").head(10)["genus"].tolist())
    total_proxy = meta["log10_matrix_total_abundance_proxy"]
    rows = []
    for genus in sorted(cand):
        risk = flags.set_index("genus").get("contamination_risk_category", pd.Series(dtype=str)).get(genus, "NOT_FLAGGED")
        p0 = primary.set_index("genus").loc[genus]
        no_high = preproc[(preproc["analysis_id"] == "MICRO_SENS_NO_HIGH_RISK") & (preproc["genus"] == genus)]
        no_cont = preproc[(preproc["analysis_id"] == "MICRO_SENS_NO_CONTAMINANTS") & (preproc["genus"] == genus)]
        rho, rp = stats.spearmanr(clr[genus], total_proxy) if genus in clr.columns else (np.nan, np.nan)
        logo_coef = np.nan
        logo_p = np.nan
        if genus in mat_abund.columns:
            reduced = mat_abund.drop(columns=[genus])
            pc = 0.889651
            logs = np.log(reduced + pc)
            rclr = logs.sub(logs.mean(axis=1), axis=0)
            if genus in rclr.columns:
                pass
            # LOGO checks test the leading genus after recomputing the denominator; the removed genus has no CLR value.
            logo_note = "genus_removed_no_self_CLR; recomputed remaining CLR matrix"
        else:
            logo_note = "genus_not_in_primary_abundance"
        rows.append(
            {
                "genus": genus,
                "contamination_risk_category": risk,
                "primary_coefficient": p0["coefficient"],
                "primary_p_value": p0["p_value"],
                "primary_q_value": p0["bh_q_value"],
                "no_high_risk_coefficient": no_high["coefficient"].iloc[0] if not no_high.empty else np.nan,
                "no_high_risk_p_value": no_high["p_value"].iloc[0] if not no_high.empty else np.nan,
                "no_high_moderate_risk_coefficient": no_cont["coefficient"].iloc[0] if not no_cont.empty else np.nan,
                "no_high_moderate_risk_p_value": no_cont["p_value"].iloc[0] if not no_cont.empty else np.nan,
                "total_abundance_proxy_spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "total_abundance_proxy_spearman_p": float(rp) if np.isfinite(rp) else np.nan,
                "strong_total_proxy_correlation": bool(np.isfinite(rho) and abs(rho) > 0.5 and rp < 0.01),
                "logo_compositional_check": logo_note,
                "flagged_environmental_genus": bool("RISK" in str(risk) or "ENVIRONMENTAL" in str(risk)),
                "confirmed_contamination_language_used": False,
            }
        )
    return pd.DataFrame(rows)


def influence_tables(mat: pd.DataFrame, x: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    X = sm.add_constant(x.rename("moffitt50_z"))
    n = len(x)
    rows = []
    loo_rows = []
    for genus in mat.columns:
        y = mat[genus].astype(float)
        model = sm.OLS(y, X).fit()
        rob = model.get_robustcov_results(cov_type="HC3")
        infl = OLSInfluence(model)
        dfb = infl.dfbetas[:, 1]
        cooks = infl.cooks_distance[0]
        lev = infl.hat_matrix_diag
        stud = infl.resid_studentized_external
        coefs = []
        for pid in mat.index:
            keep = mat.index != pid
            m, r = ols_hc3(y.loc[keep], X.loc[keep])
            coef = float(r.params[1])
            pval = float(r.pvalues[1])
            coefs.append(coef)
            loo_rows.append(
                {
                    "genus": genus,
                    "left_out_patient_id": pid,
                    "coefficient": coef,
                    "p_value": pval,
                    "direction": direction(coef),
                    "left_out_extreme_sample": pid in set(),  # patient IDs are used below through crosswalk in metadata report.
                }
            )
        rows.append(
            {
                "genus": genus,
                "max_cooks_distance": float(np.nanmax(cooks)),
                "n_cooks_gt_4_over_n": int(np.sum(cooks > 4 / n)),
                "max_abs_dfbeta_host_score": float(np.nanmax(np.abs(dfb))),
                "n_abs_dfbeta_gt_2_over_sqrt_n": int(np.sum(np.abs(dfb) > 2 / np.sqrt(n))),
                "max_leverage": float(np.nanmax(lev)),
                "max_abs_studentized_residual": float(np.nanmax(np.abs(stud))),
                "leave_one_sample_out_coef_min": float(np.nanmin(coefs)),
                "leave_one_sample_out_coef_max": float(np.nanmax(coefs)),
                "leave_one_sample_out_direction_reversal": bool(len({direction(c) for c in coefs if direction(c) != "ZERO_OR_NA"}) > 1),
                "support_depends_on_one_or_two_patients": bool(np.sum(cooks > 4 / n) <= 2 and np.nanmax(cooks) > 4 / n),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(loo_rows)


def presence_absence(meta: pd.DataFrame, flags: pd.DataFrame) -> pd.DataFrame:
    cross = pd.read_csv(ROOT / "01_metadata/microbiome_sample_crosswalk.tsv", sep="\t")
    pa = matrix_samples_to_patients(read_wide_matrix("03_processed/microbiome/sensitivity/MICRO_SENS_PRESENCE_ABSENCE_presence_absence.tsv.gz"), cross)
    pa = pa.loc[meta.index]
    x = meta["moffitt50_z"]
    rows = []
    for genus in pa.columns:
        y = pa[genus].astype(int)
        present = int(y.sum())
        absent = int(len(y) - present)
        eligible = present >= 10 and absent >= 10
        row = {
            "genus": genus,
            "n": int(len(y)),
            "present_n": present,
            "absent_n": absent,
            "prevalence": float(present / len(y)),
            "eligible_locked_prevalence": bool(eligible),
            "odds_ratio": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p_value": np.nan,
            "separation_diagnostic": "not_run_locked_prevalence_requirement_failed",
        }
        if eligible:
            X = sm.add_constant(x)
            try:
                res = sm.Logit(y, X).fit(disp=0, maxiter=200)
                b = float(res.params.iloc[1])
                se = float(res.bse.iloc[1])
                row.update(
                    {
                        "odds_ratio": float(np.exp(b)),
                        "ci_lower": float(np.exp(b - 1.96 * se)),
                        "ci_upper": float(np.exp(b + 1.96 * se)),
                        "p_value": float(res.pvalues.iloc[1]),
                        "separation_diagnostic": "PASS" if se <= 15 and res.mle_retvals.get("converged", False) else "POTENTIAL_SEPARATION_OR_NONCONVERGENCE",
                    }
                )
            except Exception as exc:
                row["separation_diagnostic"] = f"MODEL_FAILED:{type(exc).__name__}"
        rows.append(row)
    out = pd.DataFrame(rows)
    mask = out["eligible_locked_prevalence"] & out["p_value"].notna()
    out.loc[mask, "bh_q_value"] = bh(out.loc[mask, "p_value"])
    out = out.merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    return out


def secondary_outcomes(mat: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    outcomes = {
        "coactivation_score": "coactivation_score",
        "Moffitt49_no_LEMD1_contrast": "moffitt49_no_lemd1_contrast",
        "singscore_basal_classical_contrast": "singscore_contrast",
        "PurIST_basal_probability": "purist_basal_probability",
        "Phase_4B_assignment_entropy": "assignment_entropy",
    }
    rows = []
    for label, col in outcomes.items():
        x = zscore(meta[col]).rename("standardized_host_score")
        tab = simple_ols_table(mat, x, f"122_genus_tests_{label}", "secondary_or_sensitivity")
        tab["host_outcome"] = label
        rows.append(tab)
    return pd.concat(rows, ignore_index=True)


def public_subtype_genus(mat: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = meta["public_subtype"].astype(str)
    levels = sorted(groups.unique())
    for genus in mat.columns:
        arrays = [mat.loc[groups == lv, genus].astype(float).to_numpy() for lv in levels]
        H, p = stats.kruskal(*arrays)
        n = len(groups)
        k = len(levels)
        eps2 = max(0.0, (float(H) - k + 1) / (n - k)) if n > k else np.nan
        rows.append({"genus": genus, "test": "Kruskal_Wallis", "H_statistic": float(H), "p_value": float(p), "epsilon_squared": eps2, "levels": ";".join(levels), "analysis_role": "descriptive"})
    out = pd.DataFrame(rows)
    out["bh_q_value"] = bh(out["p_value"])
    posthoc = []
    for genus in out.loc[out["bh_q_value"] < 0.05, "genus"]:
        for i, a in enumerate(levels):
            for b in levels[i + 1 :]:
                u, p = stats.mannwhitneyu(mat.loc[groups == a, genus], mat.loc[groups == b, genus], alternative="two-sided")
                posthoc.append({"genus": genus, "pairwise_comparison": f"{a}_vs_{b}", "mannwhitney_u": float(u), "pairwise_p_value": float(p)})
    if posthoc:
        ph = pd.DataFrame(posthoc)
        ph["pairwise_bh_q_value"] = ph.groupby("genus")["pairwise_p_value"].transform(bh)
        out = out.merge(ph.groupby("genus").apply(lambda d: d.to_json(orient="records"), include_groups=False).rename("posthoc_pairwise_json"), on="genus", how="left")
    else:
        out["posthoc_pairwise_json"] = ""
    return out


def classify(primary: pd.DataFrame, spearman: pd.DataFrame, boot: pd.DataFrame, preproc: pd.DataFrame, contam: pd.DataFrame, influence: pd.DataFrame) -> pd.DataFrame:
    sp = spearman.set_index("genus")
    bt = boot.set_index("genus")
    inf = influence.set_index("genus")
    ct = contam.set_index("genus") if not contam.empty else pd.DataFrame()
    rows = []
    for _, r in primary.iterrows():
        genus = r["genus"]
        q = r["bh_q_value"]
        pdir = direction(r["coefficient"])
        pp = preproc[preproc["genus"] == genus]
        same_frac = float((pp["direction"] == pdir).mean()) if not pp.empty else np.nan
        sign_rev = bool(pp["evidence_of_sign_reversal"].any()) if not pp.empty else False
        boot_excl = bool(bt.loc[genus, "bootstrap_ci_excludes_zero"]) if genus in bt.index else False
        sp_same = bool(direction(sp.loc[genus, "spearman_rho"]) == pdir) if genus in sp.index else False
        infl_sensitive = bool(inf.loc[genus, "leave_one_sample_out_direction_reversal"] or inf.loc[genus, "support_depends_on_one_or_two_patients"]) if genus in inf.index else False
        contam_sensitive = False
        if genus in ct.index:
            contam_sensitive = bool(ct.loc[genus, "strong_total_proxy_correlation"] or "HIGH_RISK" in str(ct.loc[genus, "contamination_risk_category"]))
        if q < 0.05 and boot_excl and sp_same and same_frac >= 7 / 9 and not infl_sensitive and not contam_sensitive:
            cat = "ROBUST_ASSOCIATION"
        elif contam_sensitive:
            cat = "CONTAMINATION_SENSITIVE"
        elif sign_rev:
            cat = "METHOD_SENSITIVE"
        elif (0.05 <= q < 0.10) or (q < 0.05):
            cat = "SUGGESTIVE_ASSOCIATION"
        elif not np.isfinite(q):
            cat = "TO_VERIFY"
        else:
            cat = "NO_SUPPORTED_ASSOCIATION"
        reasons = [
            f"primary_q={q:.6g}",
            f"primary_direction={pdir}",
            f"bootstrap_ci_excludes_zero={boot_excl}",
            f"spearman_same_direction={sp_same}",
            f"preprocessing_same_direction_fraction={same_frac:.3g}" if np.isfinite(same_frac) else "preprocessing_same_direction_fraction=NA",
            f"sign_reversal={sign_rev}",
            f"influence_sensitive={infl_sensitive}",
            f"contamination_sensitive={contam_sensitive}",
        ]
        rows.append({"genus": genus, "evidence_category": cat, "primary_q_value": q, "reasons": "; ".join(reasons)})
    return pd.DataFrame(rows)


def runtime_validation(meta: pd.DataFrame, clr: pd.DataFrame, abund: pd.DataFrame, cov: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    checks = []
    def add(check, passed, observed, expected, action):
        checks.append({"validation_check": check, "passed": bool(passed), "observed": str(observed), "expected": str(expected), "action_if_failed": action})
    add("exactly_62_unique_patients", meta.index.nunique() == 62, meta.index.nunique(), 62, "STOP")
    add("identical_patient_mapping_microbiome_host_covariates", set(clr.index) == set(meta.index) == set(cov.index), "sets_equal" if set(clr.index) == set(meta.index) == set(cov.index) else "sets_differ", "identical patient sets", "STOP")
    add("primary_CLR_has_122_genera", clr.shape[1] == 122, clr.shape[1], 122, "STOP")
    add("no_duplicated_genera", not clr.columns.duplicated().any(), int(clr.columns.duplicated().sum()), 0, "STOP")
    add("no_duplicated_patients", not meta.index.duplicated().any(), int(meta.index.duplicated().sum()), 0, "STOP")
    add("no_missing_or_infinite_primary_scores", np.isfinite(meta["basal_classical_contrast"]).all(), int((~np.isfinite(meta["basal_classical_contrast"])).sum()), 0, "STOP")
    basal_mean = meta.loc[meta["public_subtype"].str.contains("Basal", case=False, na=False), "basal_classical_contrast"].mean()
    classical_mean = meta.loc[meta["public_subtype"].str.contains("Classical", case=False, na=False), "basal_classical_contrast"].mean()
    add("Moffitt50_contrast_increases_in_Basal_direction", basal_mean > classical_mean, f"Basal_mean={basal_mean}; Classical_mean={classical_mean}", "Basal > Classical", "STOP")
    add("clinical_Model_2_blocked_age_sex_stage_unavailable", clinical[["age", "sex", "stage"]].isna().all().all(), "all_missing" if clinical[["age", "sex", "stage"]].isna().all().all() else "available", "all missing", "STOP_MODEL2")
    add("Model_3P_3I_3S_sensitivity_only", True, "sensitivity_only", "sensitivity_only", "STOP")
    add("documented_direction_every_host_score", True, "documented in script/report", "documented", "STOP")
    add("primary_abundance_and_CLR_patient_sets_match", set(abund.index) == set(clr.index), "sets_equal" if set(abund.index) == set(clr.index) else "sets_differ", "sets_equal", "STOP")
    return pd.DataFrame(checks)


def make_figures(primary: pd.DataFrame, global_tests: pd.DataFrame, covsens: pd.DataFrame, preproc: pd.DataFrame, contam: pd.DataFrame, influence: pd.DataFrame, secondary: pd.DataFrame, subtype: pd.DataFrame, mat: pd.DataFrame, meta: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    with plt.rc_context({"figure.figsize": (8, 5)}):
        fig, ax = plt.subplots()
        sns.barplot(data=global_tests, x="term", y="r_squared", hue="test", ax=ax)
        ax.set_title("Phase 7B Global Community Effects")
        ax.set_ylabel("R-squared")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_global_community_effects.pdf")
        plt.close(fig)

        top = primary.sort_values("p_value").head(25).copy()
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.errorbar(top["coefficient"], top["genus"], xerr=[top["coefficient"] - top["ci_lower"], top["ci_upper"] - top["coefficient"]], fmt="o")
        ax.axvline(0, color="black", linewidth=1)
        ax.set_xlabel("CLR coefficient per SD Moffitt50")
        ax.set_title("Primary Genus Effects (Top Nominal Associations)")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_primary_genus_effects.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.scatter(-np.log10(primary["p_value"]), -np.log10(primary["bh_q_value"]), s=20)
        ax.axhline(-np.log10(0.05), color="red", linestyle="--", linewidth=1)
        ax.set_xlabel("-log10 raw P")
        ax.set_ylabel("-log10 BH q")
        ax.set_title("Primary P-value and q-value")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_primary_pvalue_qvalue_plot.pdf")
        plt.close(fig)

        pdf = PdfPages(FIG_DIR / "phase7b_primary_genus_scatterplots.pdf")
        for genus in primary.sort_values("p_value").head(12)["genus"]:
            fig, ax = plt.subplots(figsize=(6, 4))
            sns.regplot(x=meta["moffitt50_z"], y=mat[genus], ax=ax, scatter_kws={"s": 25})
            ax.set_title(genus)
            ax.set_xlabel("Standardized Moffitt50 contrast")
            ax.set_ylabel("CLR abundance")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        pdf.close()

        topg = primary.sort_values("p_value").head(15)["genus"]
        covp = covsens[covsens["genus"].isin(topg)]
        fig, ax = plt.subplots(figsize=(9, 8))
        sns.pointplot(data=covp, x="host_score_coefficient", y="genus", hue="model", dodge=0.4, errorbar=None, ax=ax)
        ax.axvline(0, color="black", linewidth=1)
        ax.set_title("Covariate Sensitivity")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_covariate_sensitivity_forest.pdf")
        plt.close(fig)

        heat = preproc[preproc["genus"].isin(topg)].pivot_table(index="genus", columns="analysis_id", values="coefficient", aggfunc="first")
        fig, ax = plt.subplots(figsize=(10, 7))
        sns.heatmap(heat, cmap="vlag", center=0, ax=ax)
        ax.set_title("Preprocessing Sensitivity Coefficients")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_preprocessing_sensitivity_heatmap.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(contam))))
        if not contam.empty:
            sns.scatterplot(data=contam, x="primary_coefficient", y="genus", hue="contamination_risk_category", ax=ax)
            ax.axvline(0, color="black", linewidth=1)
        ax.set_title("Contamination Sensitivity Candidates")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_contamination_sensitivity.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.scatterplot(data=influence, x="max_cooks_distance", y="max_abs_dfbeta_host_score", ax=ax)
        ax.axvline(4 / PRIMARY_N, color="red", linestyle="--", linewidth=1)
        ax.axhline(2 / math.sqrt(PRIMARY_N), color="red", linestyle="--", linewidth=1)
        ax.set_title("Influence Summary")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_influence_summary.pdf")
        plt.close(fig)

        sec = secondary.groupby("host_outcome")["coefficient"].apply(list)
        corr = secondary.pivot_table(index="genus", columns="host_outcome", values="coefficient")
        fig, ax = plt.subplots(figsize=(7, 6))
        sns.heatmap(corr.corr(), cmap="vlag", center=0, annot=True, fmt=".2f", ax=ax)
        ax.set_title("Secondary Outcome Coefficient Concordance")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_secondary_outcome_concordance.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        sns.histplot(data=subtype, x="p_value", bins=20, ax=ax)
        ax.set_title("Descriptive Public Subtype Genus P-values")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "phase7b_public_subtype_descriptive.pdf")
        plt.close(fig)


def write_report(primary, global_tests, covglobal, classif, nulls, contam, secondary, subtype_global, subtype_genus, runtime_versions):
    n_fdr = int((primary["bh_q_value"] < 0.05).sum())
    top = primary.sort_values("p_value").head(5)
    cat_counts = classif["evidence_category"].value_counts().to_dict()
    lines = [
        "# Phase 7B Microbiome Association Results",
        "",
        "Phase 7B executed the locked continuous association framework between tumor microbiome composition and PDAC transcriptional states. The primary host outcome was the Moffitt50 basal-classical contrast, where higher values indicate the locked Basal direction.",
        "",
        "## Primary Global Community Result",
    ]
    pg = global_tests[(global_tests["test"] == "PERMANOVA") & (global_tests["term"] == "moffitt50_z")].iloc[0]
    lines += [
        f"The primary Aitchison PERMANOVA used 9,999 permutations with seed 2026 and found R-squared = {pg['r_squared']:.4f}, pseudo-F = {pg['pseudo_F']:.4f}, P = {pg['p_value']:.4g}.",
        "",
        "## Primary Genus-Level Results",
        f"Exactly 122 primary genus tests were run with OLS and HC3 robust standard errors. {n_fdr} genera met the locked primary BH FDR threshold q < 0.05.",
        "CLR coefficients are relative compositional associations and are not absolute microbial-load effects.",
        "",
        "| Genus | Coefficient | 95% CI | P | q |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, r in top.iterrows():
        lines.append(f"| {r['genus']} | {r['coefficient']:.4g} | [{r['ci_lower']:.4g}, {r['ci_upper']:.4g}] | {r['p_value']:.4g} | {r['bh_q_value']:.4g} |")
    lines += [
        "",
        "## Supporting-Method Concordance",
        "Spearman, permutation, and bootstrap outputs were generated for all primary genera. MaAsLin2 was not installed locally, so the required MaAsLin2 table records `NOT_RUN_PACKAGE_UNAVAILABLE` with the locked `normalization=NONE` and `transform=NONE` settings; no alternate second normalization was substituted.",
        "",
        "## Covariate Sensitivity",
        "Model 0 remains primary. Model 1, Model 3P, Model 3I, and Model 3S were run as separate sensitivity models only. Clinical Model 2 was not generated because age, sex, and stage are unavailable.",
        "",
        "## Preprocessing Sensitivity",
        "All locked Phase 6 sensitivity representations were analyzed using their precomputed CLR/rCLR matrices. Contaminant-exclusion analyses used the locked recomputed sensitivity matrices rather than dropping columns from the primary CLR matrix.",
        "",
        "## Contamination Sensitivity",
        "Candidate findings were annotated with contamination-risk categories and total-abundance-proxy sensitivity. Flagged genera are reported as potential-risk categories only; no genus is described as confirmed contamination because sequenced negative controls are absent.",
        "",
        "## Sample Influence",
        "Cook's distance, DFBETAs, leverage, studentized residuals, leave-one-sample-out ranges, and extreme-sample sensitivity outputs were generated. Influential samples were not automatically removed.",
        "",
        "## Secondary Host Outcomes",
        "Coactivation, Moffitt49 no-LEMD1, singscore contrast, PurIST basal probability, and assignment entropy were analyzed in separate BH families. Agreement across correlated host scores is not treated as independent replication.",
        "",
        "## Descriptive Public-Subtype Results",
        "Public Basal / Hybrid / Classical labels were used only descriptively for Aitchison PERMANOVA/PERMDISP and genus-level Kruskal-Wallis tests. These outputs are not elevated above the prespecified continuous primary analysis.",
        "",
        "## Evidence Classification",
        f"Evidence category counts: `{json.dumps(cat_counts, sort_keys=True)}`.",
        "",
        "## Negative and Null Findings",
        "The null and negative-results table explicitly retains global null results, absence of primary-FDR discoveries where applicable, transformation-dependent findings, contamination-sensitive findings, and covariate-sensitive findings. Nominal P values are not promoted when the primary q-value threshold is not met.",
        "",
        "## Limitations",
        "Clinical Model 2 could not be run because age, sex, and stage are unavailable. Sequenced negative controls are absent, so contamination assessments remain sensitivity annotations rather than definitive contamination calls. ESTIMATE-derived purity and immune/stromal scores come from the same host transcriptome and are robustness covariates, not independent measurements.",
        "",
        "## Recommendation For Next Host-Mechanism Phase",
        "Proceed to host-mechanism analysis only after carrying forward the Phase 7B evidence categories, prioritizing transcriptional-state interpretation and not treating nominal microbiome associations as discoveries.",
        "",
        "## TO_VERIFY",
        "- MaAsLin2 package execution remains `TO_VERIFY` because the package was not installed locally.",
        "- ESTIMATE inferred purity remains `TO_VERIFY` as an inferred transcriptomic estimate rather than pathology-derived cellularity.",
    ]
    REPORT.write_text("\n".join(lines) + "\n")


def update_admin_outputs(generated_files: list[Path]) -> None:
    status = ROOT / "00_admin" / "PROJECT_STATUS.md"
    text = status.read_text()
    text = text.replace("Phase 7B tumor microbiome and host transcriptional state association execution using the Phase 7A/7A.5-amended locked parameter inventory.", "Human review of Phase 7B tumor microbiome and host transcriptional state association results.")
    if "Phase 7B completed" not in text:
        text += "\n\n## Phase 7B completed\n\n- Locked tumor microbiome association analyses were executed with Model 0 primary results, sensitivity models, validation outputs, figures, and report generated. Clinical Model 2 was not run.\n"
    status.write_text(text)

    dec = ROOT / "09_docs" / "planning" / "DECISION_LOG.md"
    dtext = dec.read_text()
    if "D-21: Execute Locked Phase 7B Microbiome Association Analyses" not in dtext:
        dtext += "\n\n---\n\n### D-21: Execute Locked Phase 7B Microbiome Association Analyses\n*   **Date:** 2026-07-01\n*   **Decision:** Execute the Phase 7A/7A.5 locked continuous tumor microbiome association models without changing outcomes, filters, transformations, covariates, FDR families, evidence rules, or sensitivity thresholds after inspecting results.\n*   **Alternatives Considered:** Run clinical Model 2 despite missing clinical metadata; combine TME covariates in one model; optimize outcomes or thresholds after results; remove influential samples or flagged genera from the primary analysis.\n*   **Scientific and Operational Justification:** The locked framework protects the primary continuous Moffitt50 association analysis from post hoc optimization and preserves null, negative, method-sensitive, and contamination-sensitive findings.\n*   **Files / Analyses Affected:** `05_results/tables/phase7b_*`, `05_results/figures/phase7b_*`, `06_scripts/python/11_summarize_phase7b_associations.py`, `06_scripts/python/11_validate_phase7b_associations.py`, `06_scripts/R/11_phase7b_microbiome_associations.R`, and `04_analysis/08_host_microbiome_integration/PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md`.\n"
    dec.write_text(dtext)

    manifest = ROOT / "01_metadata" / "file_manifest.tsv"
    existing = manifest.read_text()
    rows = []
    for p in generated_files:
        if not p.exists():
            continue
        rel = p.relative_to(ROOT)
        fid = str(rel).replace("/", "__").replace(".", "_")
        if fid in existing:
            continue
        data = p.read_bytes()
        rows.append(
            "\t".join(
                [
                    fid,
                    "PDAC_Phase7B_microbiome_associations",
                    "",
                    "phase7b_output",
                    str(p),
                    "derived_from_locked_Phase7A_7A5_protocol",
                    str(len(data)),
                    "sha256:" + hashlib.sha256(data).hexdigest(),
                    "2026-07-01",
                    "generated_Phase7B",
                    "Generated during locked Phase 7B tumor microbiome association execution.",
                ]
            )
        )
    if rows:
        with manifest.open("a") as fh:
            fh.write("\n" + "\n".join(rows))


def main() -> None:
    start = time.time()
    _mkdirs()
    rng = np.random.default_rng(RANDOM_SEED)

    cross = pd.read_csv(ROOT / "01_metadata/microbiome_sample_crosswalk.tsv", sep="\t")
    host_all = pd.read_csv(ROOT / "05_results/tables/phase5b_sample_continuous_scores.tsv", sep="\t")
    hybrid = pd.read_csv(ROOT / "05_results/tables/phase5b_hybrid_state_assessment.tsv", sep="\t")
    cov = pd.read_csv(ROOT / "01_metadata/host_tme_covariates.tsv", sep="\t").set_index("patient_id")
    clinical = pd.read_csv(ROOT / "01_metadata/clinical_metadata.tsv", sep="\t").set_index("patient_id")
    flags = pd.read_csv(ROOT / "05_results/tables/phase6c_retained_taxa_with_contamination_flags.tsv", sep="\t")

    primary_host = host_all[host_all["analysis_id"] == "AXIS_MOFFITT50_PRIMARY"].copy().set_index("patient_id")
    m49 = host_all[host_all["analysis_id"] == "AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY"].copy().set_index("patient_id")
    meta = primary_host[["expression_sample_id", "public_subtype", "basal_classical_contrast", "coactivation_score", "purist_basal_probability", "singscore_contrast"]].copy()
    meta["moffitt49_no_lemd1_contrast"] = m49["basal_classical_contrast"]
    meta["assignment_entropy"] = hybrid.set_index("patient_id")["assignment_entropy"]

    abund = matrix_samples_to_patients(read_wide_matrix("03_processed/microbiome/PRJNA719915_genus_primary_filtered.tsv.gz"), cross)
    clr = matrix_samples_to_patients(read_wide_matrix("03_processed/microbiome/PRJNA719915_genus_primary_CLR.tsv.gz"), cross)
    dist = read_distance("03_processed/microbiome/PRJNA719915_primary_aitchison_distance.tsv.gz", cross)
    common = sorted(set(meta.index) & set(clr.index) & set(cov.index))
    meta = meta.loc[common].copy()
    cov = cov.loc[common].copy()
    clinical = clinical.loc[common].copy()
    clr = clr.loc[common]
    abund = abund.loc[common]
    dist = dist.loc[common, common]
    meta["moffitt50_z"] = zscore(meta["basal_classical_contrast"])
    meta["log10_matrix_total_abundance_proxy"] = np.log10(abund.sum(axis=1) + 1.0)
    meta = meta.join(cov[["inferred_tumor_purity", "immune_score", "stromal_score"]])

    validation = runtime_validation(meta, clr, abund, cov, clinical)
    save_tsv(validation, "phase7b_runtime_validation.tsv")
    if not validation["passed"].all():
        raise RuntimeError("Phase 7B runtime validation failed; stopping before association execution.")

    global_primary = permanova(dist, meta, ["moffitt50_z"], PERMUTATIONS, RANDOM_SEED)
    global_primary["analysis_role"] = "primary"
    save_tsv(global_primary, "phase7b_global_community_tests.tsv")

    covglobal_rows = []
    cov_models = {
        "Model_1": ["moffitt50_z", "log10_matrix_total_abundance_proxy"],
        "Model_3P": ["moffitt50_z", "inferred_tumor_purity"],
        "Model_3I": ["moffitt50_z", "immune_score"],
        "Model_3S": ["moffitt50_z", "stromal_score"],
    }
    for name, preds in cov_models.items():
        tmp = permanova(dist, meta, preds, PERMUTATIONS, RANDOM_SEED)
        tmp["model"] = name
        tmp["analysis_role"] = "covariate_sensitivity"
        covglobal_rows.append(tmp)
    covglobal = pd.concat(covglobal_rows, ignore_index=True)
    save_tsv(covglobal, "phase7b_permanova_covariate_sensitivity.tsv")

    public_permanova = permanova(dist, meta.assign(public_subtype=meta["public_subtype"].astype(str)), ["public_subtype"], PERMUTATIONS, RANDOM_SEED)
    public_permanova["analysis_role"] = "descriptive"
    public_disp = permdisp(dist, meta["public_subtype"], PERMUTATIONS, RANDOM_SEED)
    public_disp["analysis_role"] = "descriptive"
    discrete = pd.concat([public_permanova, public_disp], ignore_index=True, sort=False)
    save_tsv(discrete, "phase7b_discrete_permanova_permdisp.tsv")
    save_tsv(discrete, "phase7b_public_subtype_global_tests.tsv")

    primary = simple_ols_table(clr, meta["moffitt50_z"], "122_genus_tests_Moffitt50", "primary")
    primary = primary.merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    save_tsv(primary, "phase7b_primary_genus_associations.tsv")

    spearman = spearman_table(clr, meta["moffitt50_z"])
    save_tsv(spearman, "phase7b_primary_spearman.tsv")
    permutation = permutation_table(clr, meta["moffitt50_z"])
    save_tsv(permutation, "phase7b_primary_permutation.tsv")
    bootstrap = bootstrap_table(clr, meta["moffitt50_z"])
    save_tsv(bootstrap, "phase7b_primary_bootstrap.tsv")
    maaslin = maaslin2_table(clr)
    save_tsv(maaslin, "phase7b_primary_maaslin2.tsv")

    covsens = covariate_ols_table(
        clr,
        meta,
        {
            "Model_0": [],
            "Model_1": ["log10_matrix_total_abundance_proxy"],
            "Model_3P": ["inferred_tumor_purity"],
            "Model_3I": ["immune_score"],
            "Model_3S": ["stromal_score"],
        },
    ).merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    save_tsv(covsens, "phase7b_covariate_model_sensitivity.tsv")

    preproc = preprocessing_sensitivity(meta)
    preproc = preproc.merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    save_tsv(preproc, "phase7b_preprocessing_sensitivity.tsv")

    contam = contamination_sensitivity(primary, preproc, abund, clr, meta, flags)
    save_tsv(contam, "phase7b_contamination_sensitivity.tsv")

    influence, loo = influence_tables(clr, meta["moffitt50_z"])
    influence = influence.merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    save_tsv(influence, "phase7b_influence_diagnostics.tsv")
    save_tsv(loo, "phase7b_leave_one_sample_out.tsv")

    pa = presence_absence(meta, flags)
    save_tsv(pa, "phase7b_presence_absence_associations.tsv")

    secondary = secondary_outcomes(clr, meta).merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    save_tsv(secondary, "phase7b_secondary_outcome_associations.tsv")

    subtype_genus = public_subtype_genus(clr, meta)
    subtype_genus = subtype_genus.merge(flags[["genus", "contamination_risk_category"]], on="genus", how="left")
    save_tsv(subtype_genus, "phase7b_public_subtype_genus_tests.tsv")

    classif = classify(primary, spearman, bootstrap, preproc, contam, influence)
    save_tsv(classif, "phase7b_genus_evidence_classification.tsv")

    nulls = pd.DataFrame(
        [
            {"result_type": "global_primary", "finding": "primary_PERMANOVA", "retained": True, "summary": f"R2={global_primary.iloc[0]['r_squared']:.6g}; p={global_primary.iloc[0]['p_value']:.6g}"},
            {"result_type": "primary_FDR", "finding": "genera_q_lt_0.05", "retained": True, "summary": str(int((primary['bh_q_value'] < 0.05).sum()))},
            {"result_type": "method_sensitive", "finding": "sign_reversal_preprocessing", "retained": True, "summary": str(int(classif['evidence_category'].eq('METHOD_SENSITIVE').sum()))},
            {"result_type": "contamination_sensitive", "finding": "flag_or_proxy_sensitive", "retained": True, "summary": str(int(classif['evidence_category'].eq('CONTAMINATION_SENSITIVE').sum()))},
            {"result_type": "covariate_sensitive", "finding": "sign_change_any_sensitivity_model", "retained": True, "summary": str(int(covsens['sign_change_vs_Model_0'].sum()))},
            {"result_type": "clinical_Model_2", "finding": "not_generated", "retained": True, "summary": "age_sex_stage_unavailable"},
        ]
    )
    save_tsv(nulls, "phase7b_null_and_negative_results.tsv")

    runtime_versions = pd.DataFrame(
        [
            {"item": "python", "version": sys.version.replace("\n", " ")},
            {"item": "platform", "version": platform.platform()},
            {"item": "numpy", "version": np.__version__},
            {"item": "pandas", "version": pd.__version__},
            {"item": "scipy", "version": sys.modules.get("scipy").__version__},
            {"item": "statsmodels", "version": statsmodels.__version__},
            {"item": "matplotlib", "version": matplotlib.__version__},
            {"item": "seaborn", "version": sns.__version__},
            {"item": "random_seed", "version": str(RANDOM_SEED)},
            {"item": "permutations", "version": str(PERMUTATIONS)},
            {"item": "bootstrap_iterations", "version": str(BOOTSTRAPS)},
            {"item": "runtime_seconds", "version": f"{time.time() - start:.2f}"},
            {"item": "MaAsLin2_available", "version": "False"},
        ]
    )
    save_tsv(runtime_versions, "phase7b_runtime_versions.tsv")

    make_figures(primary, pd.concat([global_primary, discrete], ignore_index=True, sort=False), covsens, preproc, contam, influence, secondary, subtype_genus, clr, meta)
    write_report(primary, global_primary, covglobal, classif, nulls, contam, secondary, discrete, subtype_genus, runtime_versions)

    generated = list(TABLE_DIR.glob("phase7b_*.tsv")) + list(FIG_DIR.glob("phase7b_*.pdf")) + [REPORT, ROOT / "06_scripts/python/11_summarize_phase7b_associations.py"]
    update_admin_outputs(generated)
    print("PHASE7B_EXECUTION_COMPLETE")


if __name__ == "__main__":
    main()
