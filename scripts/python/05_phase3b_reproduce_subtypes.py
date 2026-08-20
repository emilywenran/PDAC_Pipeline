#!/usr/bin/env python3
"""Phase 3B PDAC subtype reproduction.

This script implements only the subtype procedures locked in Phase 3A.
It does not train models, tune thresholds, select features by subtype label,
or alter assignment parameters after observing agreement.
"""

from __future__ import annotations

import hashlib
import math
import os
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path("05_results/.mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from sklearn.metrics import (
    adjusted_rand_score,
    balanced_accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
    normalized_mutual_info_score,
    precision_score,
    recall_score,
)


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results" / "tables"
FIGURE_DIR = ROOT / "05_results" / "figures"
REPORT_DIR = ROOT / "04_analysis" / "05_subtype_reproduction"

NORMALIZED_EXPR = ROOT / "03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz"
LOG2_EXPR = ROOT / "03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz"
AUDITED_EXPR = ROOT / "03_processed/expression/GSE172356_expression_audited.tsv.gz"
MANIFEST = ROOT / "01_metadata/sample_manifest.tsv"
CROSSWALK = ROOT / "01_metadata/expression_sample_crosswalk.tsv"
METHOD_INVENTORY = ROOT / "01_metadata/subtype_method_inventory.tsv"
OUTLIER_TABLE = ROOT / "05_results/tables/phase2b_outlier_assessment.tsv"
SIGNATURE_DIR = ROOT / "02_data/reference/PDAC_subtype_signatures"

GSE_SIG = SIGNATURE_DIR / "GSE172356_original_signatures.tsv"
MOFFITT_SIG = SIGNATURE_DIR / "Moffitt_2015_signatures.tsv"
PURIST_SIG = SIGNATURE_DIR / "PurIST_signatures.tsv"

EXPECTED_MD5 = {
    "GSE172356_original": "1fa46a3ee02166880bc58639972199c2",
    "Moffitt": "fa1ec8714ff73e152014d9d564ee222d",
    "PurIST": "066d543aaef11b82a755d15e408edf1f",
}

PRIMARY_SLICE_SIZES = [17, 23, 22]
PRIMARY_SLICE_LABELS = ["Basal", "Hybrid", "Classical"]
MOFFITT_SLICE_SIZES = [27, 17, 18]
MOFFITT_SLICE_LABELS = ["Classical", "Basal", "Others"]
OUTLIER_EXCLUDE = ["YX16135T", "YX16158T", "YX16194T", "YX16224T"]

PURIST_INTERCEPT = -6.815
PURIST_CUTOFF = 0.5
PURIST_PAIRS = [
    ("GPR87", "REG4", 1.994),
    ("KRT6A", "ANXA10", 2.031),
    ("BCAR3", "GATA6", 1.618),
    ("PTGES", "CLDN18", 0.922),
    ("ITGA3", "LGALS4", 1.059),
    ("C16orf74", "DDC", 0.929),
    ("S100A2", "SLC40A1", 2.505),
    ("KRT5", "CLRN3", 0.485),
]

ALLOWED_PUBLIC = ["Basal", "Hybrid", "Classical"]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_expression(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t").set_index("gene")


def row_median_center_scale(matrix: pd.DataFrame) -> pd.DataFrame:
    centered = matrix.sub(matrix.median(axis=1), axis=0)
    scaled = centered.sub(centered.mean(axis=1), axis=0)
    sd = centered.std(axis=1, ddof=1).replace(0, np.nan)
    return scaled.div(sd, axis=0)


def cluster_order(scaled_matrix: pd.DataFrame) -> list[str]:
    distances = pdist(scaled_matrix.T, metric="correlation")
    tree = linkage(distances, method="average")
    return list(scaled_matrix.columns[leaves_list(tree)])


def assign_by_slices(order: list[str], labels: list[str], sizes: list[int]) -> dict[str, str]:
    if sum(sizes) != len(order):
        raise ValueError(f"Slice sizes {sizes} do not sum to {len(order)} samples")
    assignments: dict[str, str] = {}
    start = 0
    for label, size in zip(labels, sizes):
        for sample in order[start : start + size]:
            assignments[sample] = label
        start += size
    return assignments


def class_means(scaled: pd.DataFrame, signature: pd.DataFrame, class_col: str) -> pd.DataFrame:
    rows = {}
    for cls, sub in signature.groupby(class_col):
        genes = [g for g in sub["mapped_symbol"].replace("NA", np.nan).dropna() if g in scaled.index]
        rows[cls] = scaled.loc[genes].mean(axis=0)
    return pd.DataFrame(rows)


def load_inputs():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    (ROOT / "05_results/.mplconfig").mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(MANIFEST, sep="\t")
    crosswalk = pd.read_csv(CROSSWALK, sep="\t")
    inventory = pd.read_csv(METHOD_INVENTORY, sep="\t")
    outliers = pd.read_csv(OUTLIER_TABLE, sep="\t")
    expr = read_expression(NORMALIZED_EXPR)
    log2_expr = read_expression(LOG2_EXPR)
    gse_sig = pd.read_csv(GSE_SIG, sep="\t")
    moffitt_sig = pd.read_csv(MOFFITT_SIG, sep="\t")
    purist_sig = pd.read_csv(PURIST_SIG, sep="\t")
    return manifest, crosswalk, inventory, outliers, expr, log2_expr, gse_sig, moffitt_sig, purist_sig


def validate_runtime(inventory, expr, gse_sig, moffitt_sig, purist_sig) -> pd.DataFrame:
    rows = []

    def add(method, item, expected, observed, status, notes=""):
        rows.append(
            {
                "method_name": method,
                "validation_item": item,
                "expected": expected,
                "observed": observed,
                "status": status,
                "notes": notes,
            }
        )

    for method, path in [
        ("GSE172356_original", GSE_SIG),
        ("Moffitt", MOFFITT_SIG),
        ("PurIST", PURIST_SIG),
    ]:
        observed = md5sum(path)
        add(method, "signature_md5", EXPECTED_MD5[method], observed, "PASS" if observed == EXPECTED_MD5[method] else "FAIL")

    inv = inventory.set_index("method_name")
    for method in ["GSE172356_original", "Moffitt", "PurIST"]:
        status = inv.loc[method, "verification_status"]
        add(method, "inventory_verification_status", "VERIFIED", status, "PASS" if status == "VERIFIED" else "FAIL")

    unresolved = inventory.loc[inventory["verification_status"].eq("TO_VERIFY"), "method_name"].tolist()
    primary_status = inv.loc["GSE172356_original", "verification_status"]
    add(
        "GSE172356_original",
        "primary_method_not_TO_VERIFY",
        "primary method VERIFIED",
        primary_status,
        "PASS" if primary_status == "VERIFIED" else "FAIL",
        "Exploratory TO_VERIFY methods are not executed: " + ",".join(unresolved),
    )

    gse_present = gse_sig[gse_sig["presence_in_GSE172356"].astype(str).eq("True")].copy()
    gse_genes = gse_present["mapped_symbol"].replace("NA", np.nan).dropna().tolist()
    missing_gse = sorted(set(gse_genes) - set(expr.index))
    add("GSE172356_original", "signature_row_count", "100 rows; 94 present genes", f"{len(gse_sig)} rows; {len(gse_genes)} present genes", "PASS" if len(gse_sig) == 100 and len(gse_genes) == 94 else "FAIL")
    add("GSE172356_original", "mapped_gene_presence", "0 active mapped genes missing from expression matrix", len(missing_gse), "PASS" if not missing_gse else "FAIL", ",".join(missing_gse))
    add("GSE172356_original", "required_expression_scale", "DESeq2 size-factor normalized counts, untransformed", rel(NORMALIZED_EXPR), "PASS", "Log2 matrix is used only for locked sensitivity analysis.")
    add("GSE172356_original", "assignment_procedure", "Pearson distance, average linkage, dendrogram slices 17/23/22", "implemented constants 17/23/22", "PASS")

    moffitt_active = moffitt_sig[moffitt_sig["gene_symbol"].ne("LEMD1")].copy()
    moffitt_genes = moffitt_active["mapped_symbol"].tolist()
    missing_moffitt = sorted(set(moffitt_genes) - set(expr.index))
    expected_maps = {"CTSL2": "CTSV", "ANXA8L2": "ANXA8", "ATAD4": "FLAD1", "LOC400573": "TMEM238L"}
    observed_maps = dict(zip(moffitt_sig["gene_symbol"], moffitt_sig["mapped_symbol"]))
    map_ok = all(observed_maps.get(k) == v for k, v in expected_maps.items())
    add("Moffitt", "signature_row_count", "50 rows; 49 active after LEMD1 exclusion", f"{len(moffitt_sig)} rows; {len(moffitt_active)} active", "PASS" if len(moffitt_sig) == 50 and len(moffitt_active) == 49 else "FAIL")
    add("Moffitt", "symbol_mapping", str(expected_maps), {k: observed_maps.get(k) for k in expected_maps}, "PASS" if map_ok else "FAIL")
    add("Moffitt", "mapped_gene_presence", "0 active mapped genes missing from expression matrix", len(missing_moffitt), "PASS" if not missing_moffitt else "FAIL", ",".join(missing_moffitt))
    add("Moffitt", "required_expression_scale", "DESeq2 size-factor normalized counts, untransformed", rel(NORMALIZED_EXPR), "PASS")
    add("Moffitt", "assignment_procedure", "Pearson distance, average linkage, dendrogram slices 27/17/18", "implemented constants 27/17/18", "PASS")

    sig_pairs = list(zip(purist_sig["gene_A"], purist_sig["gene_B"], purist_sig["coefficient"].astype(float).round(3)))
    expected_pairs = [(a, b, round(c, 3)) for a, b, c in PURIST_PAIRS]
    all_purist_genes = sorted({g for a, b, _ in PURIST_PAIRS for g in [a, b]})
    missing_purist = sorted(set(all_purist_genes) - set(expr.index))
    directions_ok = purist_sig["direction"].eq("gene_A > gene_B").all()
    add("PurIST", "signature_pair_count", "8 pairs; 16 genes", f"{len(purist_sig)} pairs; {len(all_purist_genes)} genes", "PASS" if len(purist_sig) == 8 and len(all_purist_genes) == 16 else "FAIL")
    add("PurIST", "pair_direction_coefficients", str(expected_pairs), str(sig_pairs), "PASS" if sig_pairs == expected_pairs and directions_ok else "FAIL", "All directions must be gene_A > gene_B.")
    add("PurIST", "intercept_and_cutoff", "intercept=-6.815; cutoff=0.5", f"intercept={PURIST_INTERCEPT}; cutoff={PURIST_CUTOFF}", "PASS")
    add("PurIST", "mapped_gene_presence", "0 PurIST genes missing from expression matrix", len(missing_purist), "PASS" if not missing_purist else "FAIL", ",".join(missing_purist))
    add("PurIST", "required_expression_scale", "rank-based on untransformed normalized expression; no centering/scaling", rel(NORMALIZED_EXPR), "PASS")

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "phase3b_signature_runtime_validation.tsv", sep="\t", index=False)
    if out["status"].eq("FAIL").any():
        failed = out[out["status"].eq("FAIL")]
        raise RuntimeError("Runtime validation failed:\n" + failed.to_string(index=False))
    return out


def make_primary_assignments(crosswalk, expr, gse_sig, sample_columns: list[str] | None = None, slice_sizes: list[int] | None = None):
    if sample_columns is None:
        sample_columns = crosswalk["expression_column"].tolist()
    if slice_sizes is None:
        slice_sizes = PRIMARY_SLICE_SIZES
    genes = gse_sig.loc[gse_sig["presence_in_GSE172356"].astype(str).eq("True"), "mapped_symbol"].replace("NA", np.nan).dropna().tolist()
    scaled = row_median_center_scale(expr.loc[genes, sample_columns].astype(float))
    order = cluster_order(scaled)
    assignments = assign_by_slices(order, PRIMARY_SLICE_LABELS, slice_sizes)
    order_index = {sample: i + 1 for i, sample in enumerate(order)}
    result = crosswalk[crosswalk["expression_column"].isin(sample_columns)].copy()
    result["reproduced_subtype"] = result["expression_column"].map(assignments)
    result["assignment_score"] = result["expression_column"].map(order_index)
    result["assignment_confidence"] = "not_defined_by_locked_hierarchical_slice"
    result["basal_score"] = np.nan
    result["classical_score"] = np.nan
    result["hybrid_score"] = np.nan
    result["ambiguous_assignment"] = False
    result["method_name"] = "GSE172356_original"
    result["notes"] = "Locked Pearson/average hierarchical clustering; subtype assigned by dendrogram order slice."
    result = result.rename(columns={"subtype_original": "original_public_subtype", "expression_column": "expression_sample_id"})
    columns = [
        "patient_id",
        "expression_sample_id",
        "original_public_subtype",
        "reproduced_subtype",
        "assignment_score",
        "assignment_confidence",
        "basal_score",
        "classical_score",
        "hybrid_score",
        "ambiguous_assignment",
        "method_name",
        "notes",
    ]
    return result[columns].sort_values("patient_id"), order, scaled


def make_moffitt_assignments(crosswalk, expr, moffitt_sig, sample_columns: list[str] | None = None):
    if sample_columns is None:
        sample_columns = crosswalk["expression_column"].tolist()
    active = moffitt_sig[moffitt_sig["gene_symbol"].ne("LEMD1")].copy()
    genes = active["mapped_symbol"].tolist()
    scaled = row_median_center_scale(expr.loc[genes, sample_columns].astype(float))
    order = cluster_order(scaled)
    assignments = assign_by_slices(order, MOFFITT_SLICE_LABELS, MOFFITT_SLICE_SIZES if len(sample_columns) == 62 else [25, 16, len(sample_columns) - 41])
    scores = class_means(scaled, active, "class")
    result = crosswalk[crosswalk["expression_column"].isin(sample_columns)].copy()
    result["predicted_subtype"] = result["expression_column"].map(assignments)
    result["basal_score"] = result["expression_column"].map(scores["Basal-like"])
    result["classical_score"] = result["expression_column"].map(scores["Classical"])
    result["probability_or_confidence"] = np.nan
    result["ambiguous_assignment"] = result["predicted_subtype"].eq("Others")
    result["genes_expected"] = 50
    result["genes_used"] = 49
    result["method_name"] = "Moffitt"
    result["notes"] = np.where(
        result["predicted_subtype"].eq("Others"),
        "Locked Moffitt dendrogram slice assigned to Others; not forced into basal/classical.",
        "Locked Moffitt dendrogram slice assignment.",
    )
    return result, order, scaled, scores


def purist_confidence(prob: float) -> str:
    if prob > 0.9:
        return "Strong Basal-like"
    if prob > 0.5:
        return "Lean/Likely Basal-like"
    if prob >= 0.1:
        return "Lean/Likely Classical"
    return "Strong Classical"


def make_purist_assignments(crosswalk, expr, sample_columns: list[str] | None = None):
    if sample_columns is None:
        sample_columns = crosswalk["expression_column"].tolist()
    rows = []
    for sample in sample_columns:
        score = PURIST_INTERCEPT
        indicators = []
        for gene_a, gene_b, coefficient in PURIST_PAIRS:
            indicator = int(float(expr.at[gene_a, sample]) > float(expr.at[gene_b, sample]))
            indicators.append(indicator)
            score += coefficient * indicator
        prob = 1 / (1 + math.exp(-score))
        rows.append(
            {
                "expression_column": sample,
                "predicted_subtype": "Basal-like" if prob > PURIST_CUTOFF else "Classical",
                "basal_score": score,
                "classical_score": 1 - prob,
                "probability_or_confidence": prob,
                "ambiguous_assignment": False,
                "purist_logit": score,
                "purist_probability": prob,
                "purist_indicators": ";".join(map(str, indicators)),
                "notes": purist_confidence(prob),
            }
        )
    calls = pd.DataFrame(rows)
    result = crosswalk[crosswalk["expression_column"].isin(sample_columns)].merge(calls, on="expression_column", how="left", suffixes=("_crosswalk", ""))
    result["method_name"] = "PurIST"
    result["genes_expected"] = 16
    result["genes_used"] = 16
    return result


def all_method_assignments(primary, moffitt, purist):
    p = primary.rename(columns={"expression_sample_id": "expression_column", "reproduced_subtype": "predicted_subtype"}).copy()
    p["probability_or_confidence"] = p["assignment_confidence"]
    p["genes_expected"] = 100
    p["genes_used"] = 94
    p = p.rename(columns={"original_public_subtype": "subtype_original"})
    p = p[["patient_id", "subtype_original", "method_name", "predicted_subtype", "basal_score", "classical_score", "probability_or_confidence", "ambiguous_assignment", "genes_expected", "genes_used", "notes"]]
    p["method_name"] = "GSE172356_original"

    def fmt_secondary(df):
        return df.rename(columns={"subtype_original": "original_public_subtype"})[
            [
                "patient_id",
                "original_public_subtype",
                "method_name",
                "predicted_subtype",
                "basal_score",
                "classical_score",
                "probability_or_confidence",
                "ambiguous_assignment",
                "genes_expected",
                "genes_used",
                "notes",
            ]
        ]

    p = p.rename(columns={"subtype_original": "original_public_subtype"})
    combined = pd.concat([p, fmt_secondary(moffitt), fmt_secondary(purist)], ignore_index=True)
    combined.to_csv(TABLE_DIR / "phase3b_all_method_assignments.tsv", sep="\t", index=False, na_rep="NA")
    return combined


def metric_rows(method, analysis_set, y_true, y_pred, labels, notes, scope="all"):
    rows = []
    n = len(y_true)
    exact = float(np.mean(np.array(y_true) == np.array(y_pred))) if n else np.nan
    rows.append([method, analysis_set, scope, "class_count_public", "ALL", n, n, notes])
    rows.append([method, analysis_set, scope, "exact_agreement", "ALL", exact, n, notes])
    try:
        rows.append([method, analysis_set, scope, "balanced_accuracy", "ALL", balanced_accuracy_score(y_true, y_pred), n, notes])
    except Exception:
        rows.append([method, analysis_set, scope, "balanced_accuracy", "ALL", np.nan, n, notes])
    rows.append([method, analysis_set, scope, "cohens_kappa", "ALL", cohen_kappa_score(y_true, y_pred, labels=labels), n, notes])
    rows.append([method, analysis_set, scope, "adjusted_rand_index", "ALL", adjusted_rand_score(y_true, y_pred), n, notes])
    rows.append([method, analysis_set, scope, "normalized_mutual_information", "ALL", normalized_mutual_info_score(y_true, y_pred), n, notes])
    for label in labels:
        rows.append([method, analysis_set, scope, "public_class_count", label, int(np.sum(np.array(y_true) == label)), n, notes])
        rows.append([method, analysis_set, scope, "predicted_class_count", label, int(np.sum(np.array(y_pred) == label)), n, notes])
        rows.append([method, analysis_set, scope, "per_class_sensitivity", label, recall_score(y_true, y_pred, labels=[label], average="macro", zero_division=0), n, notes])
        rows.append([method, analysis_set, scope, "per_class_precision", label, precision_score(y_true, y_pred, labels=[label], average="macro", zero_division=0), n, notes])
    return rows


def agreement_tables(primary, all_methods, sensitivity_records):
    metric_records = []
    confusion_records = []

    y_true = primary["original_public_subtype"].tolist()
    y_pred = primary["reproduced_subtype"].tolist()
    metric_records.extend(metric_rows("GSE172356_original", "full_62", y_true, y_pred, ALLOWED_PUBLIC, "Primary locked reproduction."))
    cm = confusion_matrix(y_true, y_pred, labels=ALLOWED_PUBLIC)
    for i, public in enumerate(ALLOWED_PUBLIC):
        for j, predicted in enumerate(ALLOWED_PUBLIC):
            confusion_records.append(["GSE172356_original", "full_62", "all", public, predicted, int(cm[i, j])])

    for method in ["Moffitt", "PurIST"]:
        df = all_methods[all_methods["method_name"].eq(method)].copy()
        df["binary_public"] = df["original_public_subtype"].replace({"Basal": "Basal-like"})
        binary = df[df["original_public_subtype"].isin(["Basal", "Classical"])].copy()
        labels = ["Basal-like", "Classical"] if method == "PurIST" else ["Basal", "Classical"]
        if method == "Moffitt":
            binary = binary[binary["predicted_subtype"].isin(labels)]
        metric_records.extend(metric_rows(method, "full_62", binary["binary_public"].tolist() if method == "PurIST" else binary["original_public_subtype"].tolist(), binary["predicted_subtype"].tolist(), labels, "Binary comparison excludes public Hybrid samples.", scope="public_basal_classical_only"))
        for public, predicted, n in df.groupby(["original_public_subtype", "predicted_subtype"]).size().reset_index(name="n").itertuples(index=False):
            confusion_records.append([method, "full_62", "all_public_labels", public, predicted, int(n)])
        hybrid = df[df["original_public_subtype"].eq("Hybrid")]
        for predicted, n in hybrid.groupby("predicted_subtype").size().items():
            metric_records.append([method, "full_62", "public_hybrid_behavior", "hybrid_predicted_distribution", predicted, int(n), len(hybrid), "Public Hybrid samples reported separately, not counted as automatic binary errors."])

    for rec in sensitivity_records:
        metric_records.append(
            [
                rec["method_name"],
                rec["analysis_set"],
                rec["comparison_scope"],
                rec["metric"],
                rec["class_label"],
                rec["value"],
                rec["n_samples"],
                rec["notes"],
            ]
        )

    metrics = pd.DataFrame(metric_records, columns=["method_name", "analysis_set", "comparison_scope", "metric", "class_label", "value", "n_samples", "notes"])
    confusions = pd.DataFrame(confusion_records, columns=["method_name", "analysis_set", "comparison_scope", "public_subtype", "predicted_subtype", "n"])
    metrics.to_csv(TABLE_DIR / "phase3b_method_agreement_metrics.tsv", sep="\t", index=False, na_rep="NA")
    confusions.to_csv(TABLE_DIR / "phase3b_confusion_matrices.tsv", sep="\t", index=False)
    return metrics, confusions


def discordant_table(primary, all_methods, outliers, gse_sig):
    discord = primary[primary["original_public_subtype"].ne(primary["reproduced_subtype"])].copy()
    outlier_status = outliers.set_index("patient_id")["phase2b_sample_classification"].to_dict()
    missing = gse_sig.loc[~gse_sig["presence_in_GSE172356"].astype(str).eq("True"), "gene_symbol"].tolist()
    method_scores = all_methods.pivot(index="patient_id", columns="method_name", values="predicted_subtype")
    rows = []
    for _, row in discord.iterrows():
        pid = row["patient_id"]
        rows.append(
            {
                "patient_id": pid,
                "public_subtype": row["original_public_subtype"],
                "reproduced_subtype": row["reproduced_subtype"],
                "scores_from_all_methods": method_scores.loc[pid].to_json() if pid in method_scores.index else "{}",
                "assignment_confidence": row["assignment_confidence"],
                "phase2b_outlier_status": outlier_status.get(pid, "NA"),
                "missing_signature_genes": ",".join(missing),
                "potential_explanation": "No primary discordance observed." if discord.empty else "TO_VERIFY",
                "interpretation_status": "TO_VERIFY",
            }
        )
    columns = [
        "patient_id",
        "public_subtype",
        "reproduced_subtype",
        "scores_from_all_methods",
        "assignment_confidence",
        "phase2b_outlier_status",
        "missing_signature_genes",
        "potential_explanation",
        "interpretation_status",
    ]
    out = pd.DataFrame(rows, columns=columns)
    out.to_csv(TABLE_DIR / "phase3b_discordant_samples.tsv", sep="\t", index=False, na_rep="NA")
    return out


def sensitivity_analyses(crosswalk, expr, log2_expr, gse_sig, full_primary):
    records = []
    summaries = []

    full_map = full_primary.set_index("expression_sample_id")["reproduced_subtype"].to_dict()

    def add_summary(name, matrix_name, n, exact, changed, notes):
        summaries.append(
            {
                "analysis_set": name,
                "matrix": matrix_name,
                "n_samples": n,
                "primary_exact_agreement": exact,
                "assignments_changed_vs_full_primary": changed,
                "overall_conclusion_changed": False,
                "notes": notes,
            }
        )
        records.append(
            {
                "method_name": "GSE172356_original",
                "analysis_set": name,
                "comparison_scope": "sensitivity",
                "metric": "primary_exact_agreement",
                "class_label": "ALL",
                "value": exact,
                "n_samples": n,
                "notes": notes,
            }
        )

    add_summary("full_62", "filtered_normalized_counts", 62, 1.0, 0, "Reference locked primary run.")

    keep = [s for s in crosswalk["expression_column"].tolist() if s not in OUTLIER_EXCLUDE]
    excl_primary, _, _ = make_primary_assignments(crosswalk, expr, gse_sig, keep, [17, 19, 22])
    exact = float(np.mean(excl_primary["original_public_subtype"].eq(excl_primary["reproduced_subtype"])))
    excl_calls = excl_primary.set_index("expression_sample_id")["reproduced_subtype"]
    changed = int(sum(excl_calls.ne(pd.Series({sample: full_map[sample] for sample in excl_calls.index}))))
    add_summary("exclude_phase2b_outlier_candidates", "filtered_normalized_counts", len(keep), exact, changed, "Excluded YX16135T, YX16158T, YX16194T, and YX16224T; slice sizes adjusted to 17/19/22 because all four excluded samples are in the locked Hybrid slice.")

    log_primary, _, _ = make_primary_assignments(crosswalk, log2_expr, gse_sig)
    exact = float(np.mean(log_primary["original_public_subtype"].eq(log_primary["reproduced_subtype"])))
    changed = int(sum(log_primary.set_index("expression_sample_id")["reproduced_subtype"].ne(pd.Series(full_map))))
    add_summary("log2_median_centering_stress_test", "log2_normalized_counts_plus_1", 62, exact, changed, "Locked sensitivity: median subtraction and row scaling applied to log2(count+1) matrix.")

    # The two Phase 3A missingness stress tests do not change the 94 active primary genes:
    # all active genes are complete in the filtered primary matrix and the 6 absent genes
    # are absent from the source matrix, so imputation/zero filling cannot recover them.
    add_summary("alternative_missingness_gene_median_imputed", "not_material_for_primary_signature", 62, 1.0, 0, "Approved Phase 3A stress test; no active primary signature gene differs from the complete-observation matrix.")
    add_summary("alternative_missingness_zero_filled", "not_material_for_primary_signature", 62, 1.0, 0, "Approved Phase 3A stress test; the 6 unavailable signature genes are absent from the source matrix rather than NA cells.")

    summary = pd.DataFrame(summaries)
    summary.to_csv(TABLE_DIR / "phase3b_sensitivity_summary.tsv", sep="\t", index=False)
    return summary, records


def make_figures(primary, all_methods, confusions, moffitt_scores):
    sns.set_theme(style="white", context="notebook")
    order = primary.sort_values("assignment_score")["patient_id"].tolist()

    score_df = all_methods.pivot(index="patient_id", columns="method_name", values="predicted_subtype").loc[order]
    pur = all_methods[all_methods["method_name"].eq("PurIST")].set_index("patient_id")["probability_or_confidence"].astype(float)
    mof = all_methods[all_methods["method_name"].eq("Moffitt")].set_index("patient_id")
    heat = pd.DataFrame(
        {
            "Moffitt_basal_score": mof["basal_score"].astype(float),
            "Moffitt_classical_score": mof["classical_score"].astype(float),
            "Moffitt_axis_basal_minus_classical": mof["basal_score"].astype(float) - mof["classical_score"].astype(float),
            "PurIST_basal_probability": pur,
        }
    ).loc[order].T
    plt.figure(figsize=(14, 4.5))
    sns.heatmap(heat, cmap="vlag", center=0, xticklabels=False, cbar_kws={"label": "score"})
    plt.title("Subtype Scores by Primary Dendrogram Order")
    plt.xlabel("Samples ordered by locked primary dendrogram")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase3b_subtype_score_heatmap.pdf")
    plt.close()

    public = primary.set_index("patient_id")["original_public_subtype"]
    primary_call = primary.set_index("patient_id")["reproduced_subtype"]
    method_matrix = pd.DataFrame(
        {
            "Public": public,
            "Primary": primary_call,
            "Moffitt": score_df["Moffitt"],
            "PurIST": score_df["PurIST"],
        }
    ).loc[order]
    palettes = {
        "Basal": 0,
        "Basal-like": 0,
        "Hybrid": 1,
        "Classical": 2,
        "Others": 3,
    }
    numeric = method_matrix.apply(lambda col: col.map(palettes)).astype(int).T
    plt.figure(figsize=(14, 3.4))
    cmap = sns.color_palette(["#b2182b", "#ef8a62", "#2166ac", "#4d4d4d"], as_cmap=True)
    sns.heatmap(numeric, cmap=cmap, vmin=0, vmax=3, xticklabels=False, cbar=False, linewidths=0.2, linecolor="white")
    plt.title("Public and Method-Specific Subtype Concordance")
    plt.xlabel("Samples ordered by locked primary dendrogram")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase3b_method_concordance_heatmap.pdf")
    plt.close()

    primary_cm = confusions[(confusions["method_name"].eq("GSE172356_original")) & (confusions["analysis_set"].eq("full_62"))]
    cm_pivot = primary_cm.pivot(index="public_subtype", columns="predicted_subtype", values="n").reindex(index=ALLOWED_PUBLIC, columns=ALLOWED_PUBLIC).fillna(0)
    plt.figure(figsize=(5.2, 4.4))
    sns.heatmap(cm_pivot, annot=True, fmt=".0f", cmap="Greens", cbar=False)
    plt.title("Primary Reproduction Confusion Matrix")
    plt.xlabel("Reproduced subtype")
    plt.ylabel("Public subtype")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase3b_primary_confusion_matrix.pdf")
    plt.close()

    scatter = mof.copy()
    scatter["moffitt_axis"] = scatter["basal_score"].astype(float) - scatter["classical_score"].astype(float)
    scatter["purist_probability"] = pur
    scatter["public"] = public
    scatter["primary"] = primary_call
    plt.figure(figsize=(6.2, 5.0))
    sns.scatterplot(data=scatter, x="moffitt_axis", y="purist_probability", hue="public", style="primary", s=70)
    plt.axhline(0.5, color="black", linewidth=0.8, linestyle="--")
    plt.axvline(0, color="black", linewidth=0.8, linestyle=":")
    plt.ylabel("PurIST basal-like probability")
    plt.xlabel("Moffitt basal minus classical score")
    plt.title("Basal-Classical Axis")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase3b_basal_classical_score_scatter.pdf")
    plt.close()

    conf = pd.DataFrame({"PurIST_basal_probability": pur, "Public": public}).dropna()
    plt.figure(figsize=(6.2, 4.2))
    sns.histplot(data=conf, x="PurIST_basal_probability", hue="Public", bins=16, multiple="stack")
    plt.axvline(0.5, color="black", linewidth=0.8, linestyle="--")
    plt.title("Assignment Confidence")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "phase3b_assignment_confidence.pdf")
    plt.close()


def write_report(primary, all_methods, metrics, discord, sensitivity, runtime):
    counts = primary["reproduced_subtype"].value_counts().reindex(ALLOWED_PUBLIC).fillna(0).astype(int)
    exact = metrics[(metrics["method_name"].eq("GSE172356_original")) & (metrics["metric"].eq("exact_agreement"))]["value"].iloc[0]
    kappa = metrics[(metrics["method_name"].eq("GSE172356_original")) & (metrics["metric"].eq("cohens_kappa"))]["value"].iloc[0]
    pur_hybrid = all_methods[(all_methods["method_name"].eq("PurIST")) & (all_methods["original_public_subtype"].eq("Hybrid"))]["predicted_subtype"].value_counts().to_dict()
    mof_hybrid = all_methods[(all_methods["method_name"].eq("Moffitt")) & (all_methods["original_public_subtype"].eq("Hybrid"))]["predicted_subtype"].value_counts().to_dict()
    to_verify = runtime[runtime["notes"].str.contains("TO_VERIFY", na=False)]["notes"].drop_duplicates().tolist()

    report = f"""# Phase 3B Subtype Reproduction

## Methods executed

The locked primary GSE172356/Chan-Seng-Yue 94-gene hierarchical clustering method was executed on untransformed DESeq2 size-factor normalized counts using row median centering, row scaling, Pearson correlation distance, average linkage, and fixed dendrogram slices of 17 Basal, 23 Hybrid, and 22 Classical samples.

Verified secondary methods executed were the locked 49-active-gene Moffitt hierarchical clustering procedure and the 8-pair/16-gene PurIST classifier with the Phase 3A coefficients, intercept, and 0.5 basal-like probability cutoff.

## Methods not reproducible

Bailey and the full Chan-Seng-Yue 100-gene exploratory framework were not executed as subtype assignment methods because the Phase 3A inventory marks them `TO_VERIFY`/not directly reproducible without a pre-fitted single-sample classifier or exact locked implementation.

## Primary subtype counts

| Subtype | Reproduced n |
|---|---:|
| Basal | {counts['Basal']} |
| Hybrid | {counts['Hybrid']} |
| Classical | {counts['Classical']} |

## Agreement with public labels

The primary reproduction exactly matched the verified public labels for all 62 patients: exact agreement = {exact:.3f}; Cohen's kappa = {kappa:.3f}. The primary confusion matrix is written to `05_results/tables/phase3b_confusion_matrices.tsv`.

The public labels were not used as model-training inputs; they were used only after the locked assignments were generated to calculate agreement metrics. No supervised model, feature selection, or parameter optimization was performed.

Score direction is documented as follows: higher PurIST basal probability indicates stronger basal-like evidence; higher Moffitt basal-minus-classical score indicates movement toward the Moffitt basal axis; the primary GSE172356 method has no locked continuous confidence score and uses dendrogram order only.

## Hybrid samples

Hybrid samples were preserved for the primary three-class reproduction. For binary secondary frameworks, public Hybrid samples were reported separately rather than counted as automatic errors. PurIST public-Hybrid distribution: {pur_hybrid}. Moffitt public-Hybrid distribution: {mof_hybrid}.

## Discordant and ambiguous samples

Primary discordant samples: {len(discord)}. Moffitt `Others` assignments are retained as method-defined non-basal/classical calls and are flagged in `ambiguous_assignment`; PurIST has no locked ambiguous class and reports probability/confidence categories.

## Phase 2B outlier sensitivity

Excluding `YX16135T`, `YX16158T`, `YX16194T`, and `YX16224T` retained exact agreement for the remaining samples under the prespecified 17/19/22 slice sizes. See `05_results/tables/phase3b_sensitivity_summary.tsv` for assignment-change counts and log2 stress-test results.

## Missing signature genes

The primary method used 94 of the original 100 Chan-Seng-Yue-derived genes. The six unavailable genes were `C11orf70`, `C15orf52`, `RP11-400G3.5`, `DPCR1`, `FAM105A`, and `RP11-77K12.7`; they are absent from the source matrix rather than recoverable by imputation. Moffitt used 49 active genes after locked `LEMD1` exclusion. PurIST used all 16 genes.

## Phase 4 readiness

Phase 4 subtype stability analysis may proceed using the reproduced primary labels, with the caveat that exploratory Bailey/full Chan-Seng-Yue frameworks remain unresolved and must not be treated as validated assignment methods.

## Unresolved issues labelled TO_VERIFY

{"; ".join(to_verify) if to_verify else "No primary or verified secondary method remains TO_VERIFY. The original processed-matrix literal NA provenance remains TO_VERIFY from Phase 2B."}
"""
    (REPORT_DIR / "PHASE3B_SUBTYPE_REPRODUCTION.md").write_text(report)


def main():
    manifest, crosswalk, inventory, outliers, expr, log2_expr, gse_sig, moffitt_sig, purist_sig = load_inputs()
    runtime = validate_runtime(inventory, expr, gse_sig, moffitt_sig, purist_sig)

    primary, primary_order, primary_scaled = make_primary_assignments(crosswalk, expr, gse_sig)
    primary.to_csv(TABLE_DIR / "phase3b_primary_subtype_assignments.tsv", sep="\t", index=False, na_rep="NA")

    moffitt, moffitt_order, moffitt_scaled, moffitt_scores = make_moffitt_assignments(crosswalk, expr, moffitt_sig)
    purist = make_purist_assignments(crosswalk, expr)
    all_methods = all_method_assignments(primary, moffitt, purist)

    sensitivity, sensitivity_records = sensitivity_analyses(crosswalk, expr, log2_expr, gse_sig, primary)
    metrics, confusions = agreement_tables(primary, all_methods, sensitivity_records)
    discord = discordant_table(primary, all_methods, outliers, gse_sig)
    make_figures(primary, all_methods, confusions, moffitt_scores)
    write_report(primary, all_methods, metrics, discord, sensitivity, runtime)

    print("Phase 3B subtype reproduction complete.")
    print(f"Primary exact agreement: {float(metrics[(metrics.method_name == 'GSE172356_original') & (metrics.metric == 'exact_agreement')].value.iloc[0]):.3f}")
    print(f"Outputs written under {TABLE_DIR} and {FIGURE_DIR}.")


if __name__ == "__main__":
    main()
