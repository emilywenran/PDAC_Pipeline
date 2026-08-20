#!/usr/bin/env python3
"""Validate Phase 4B subtype-stability outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "05_results" / "tables"
FIGURES = ROOT / "05_results" / "figures"

ANALYSES = {
    "STAB_CSY_PRIMARY",
    "STAB_CSY_LOG2",
    "STAB_UNSUP_HVG",
    "STAB_CSY_OUTLIER_EXCL",
    "STAB_HVG_OUTLIER_EXCL",
    "STAB_CSY_FEAT_RESAMP",
    "STAB_CSY_IMPUTED",
    "STAB_HVG_VAR_FILTER",
}
K_VALUES = {2, 3, 4, 5, 6}

REQUIRED_TABLES = {
    "phase4b_cluster_stability_metrics.tsv": {"analysis_id", "candidate_K", "PAC", "overall_mean_silhouette", "mean_Jaccard_stability", "prediction_strength"},
    "phase4b_cluster_size_summary.tsv": {"analysis_id", "candidate_K", "cluster_id", "cluster_size", "within_cluster_consensus", "bootstrap_jaccard_stability"},
    "phase4b_k_selection_summary.tsv": {"analysis_id", "preferred_K", "preferred_K_basis"},
    "phase4b_sample_stability.tsv": {"analysis_id", "candidate_K", "sample_id", "final_cluster_assignment", "item_consensus", "assignment_entropy", "silhouette_width"},
    "phase4b_sample_assignment_probabilities.tsv": {"analysis_id", "candidate_K", "sample_id", "cluster_id", "assignment_probability"},
    "phase4b_public_label_comparison.tsv": {"analysis_id", "candidate_K", "adjusted_rand_index", "normalized_mutual_information", "cohens_kappa", "confusion_matrix", "per_class_agreement"},
    "phase4b_cluster_label_crosswalk.tsv": {"analysis_id", "candidate_K", "cluster_id", "aligned_public_label", "overlap_n"},
    "phase4b_hybrid_stability_assessment.tsv": {"sample_id", "public_subtype", "analysis_id", "assignment_entropy", "purist_basal_probability", "distance_to_basal_centroid", "distance_to_classical_centroid", "interpretation_category"},
    "phase4b_sensitivity_comparison.tsv": {"comparison_id", "reference_analysis", "sensitivity_analysis", "candidate_K", "assignment_changes", "ARI_between_analyses", "PAC_change", "silhouette_change"},
    "phase4b_recurrently_unstable_samples.tsv": {"sample_id", "public_subtype", "unstable_analysis_count", "recurrently_unstable"},
}

REQUIRED_FIGURES = [
    "phase4b_consensus_matrix_primary_K2.pdf",
    "phase4b_consensus_matrix_primary_K3.pdf",
    "phase4b_consensus_matrix_primary_K4.pdf",
    "phase4b_consensus_cdf.pdf",
    "phase4b_pac_by_K.pdf",
    "phase4b_silhouette_by_K.pdf",
    "phase4b_jaccard_stability_by_K.pdf",
    "phase4b_sample_item_consensus.pdf",
    "phase4b_sample_assignment_entropy.pdf",
    "phase4b_analysis_concordance_heatmap.pdf",
    "phase4b_hybrid_stability_summary.pdf",
    "phase4b_basal_classical_axis_with_clusters.pdf",
]


def fail(message: str) -> None:
    raise SystemExit(f"VALIDATION_FAILED: {message}")


def check_file(path: Path) -> None:
    if not path.exists():
        fail(f"missing file {path.relative_to(ROOT)}")
    if path.stat().st_size == 0:
        fail(f"empty file {path.relative_to(ROOT)}")


def check_analysis_k(df: pd.DataFrame, name: str) -> None:
    if set(df["analysis_id"]) != ANALYSES:
        fail(f"{name} analysis IDs mismatch: {sorted(set(df['analysis_id']))}")
    for aid, sub in df.groupby("analysis_id"):
        if set(sub["candidate_K"].astype(int)) != K_VALUES:
            fail(f"{name} K coverage mismatch for {aid}: {sorted(set(sub['candidate_K']))}")


def main() -> None:
    loaded: dict[str, pd.DataFrame] = {}
    for filename, columns in REQUIRED_TABLES.items():
        path = TABLES / filename
        check_file(path)
        df = pd.read_csv(path, sep="\t")
        missing = columns - set(df.columns)
        if missing:
            fail(f"{filename} missing columns {sorted(missing)}")
        loaded[filename] = df

    for figure in REQUIRED_FIGURES:
        path = FIGURES / figure
        check_file(path)
        if path.stat().st_size < 1000:
            fail(f"figure appears too small: {path.relative_to(ROOT)}")

    check_analysis_k(loaded["phase4b_cluster_stability_metrics.tsv"], "cluster metrics")
    check_analysis_k(loaded["phase4b_public_label_comparison.tsv"], "public comparison")

    ksel = loaded["phase4b_k_selection_summary.tsv"]
    if set(ksel["analysis_id"]) != ANALYSES or not set(ksel["preferred_K"].astype(int)).issubset(K_VALUES):
        fail("K-selection summary has invalid analysis IDs or preferred K values")

    sample = loaded["phase4b_sample_stability.tsv"]
    check_analysis_k(sample[["analysis_id", "candidate_K"]].drop_duplicates(), "sample stability")
    if sample[["assignment_entropy", "silhouette_width", "assignment_confidence"]].isna().any().any():
        fail("sample stability contains missing entropy, silhouette, or confidence values")
    sizes = loaded["phase4b_cluster_size_summary.tsv"]
    sample_with_size = sample.merge(
        sizes[["analysis_id", "candidate_K", "cluster_id", "cluster_size"]],
        left_on=["analysis_id", "candidate_K", "final_cluster_assignment"],
        right_on=["analysis_id", "candidate_K", "cluster_id"],
        how="left",
    )
    bad_item = sample_with_size[sample_with_size["item_consensus"].isna() & (sample_with_size["cluster_size"] > 1)]
    if not bad_item.empty:
        fail("sample item consensus is missing for non-singleton clusters")
    if not sample["assignment_confidence"].between(0, 1).all():
        fail("assignment confidence outside [0, 1]")

    probs = loaded["phase4b_sample_assignment_probabilities.tsv"]
    sums = probs.groupby(["analysis_id", "candidate_K", "sample_id"])["assignment_probability"].sum().reset_index()
    bad = sums[(sums.assignment_probability < 0.99) | (sums.assignment_probability > 1.01)]
    if not bad.empty:
        fail("assignment probabilities do not sum to approximately 1 for every sample-analysis-K")

    metrics = loaded["phase4b_cluster_stability_metrics.tsv"]
    bounded = ["PAC", "mean_Jaccard_stability", "min_cluster_Jaccard_stability", "prediction_strength"]
    for col in bounded:
        if not metrics[col].between(0, 1).all():
            fail(f"{col} outside [0, 1]")

    report = ROOT / "04_analysis" / "06_subtype_stability" / "PHASE4B_SUBTYPE_STABILITY_RESULTS.md"
    check_file(report)
    text = report.read_text()
    for phrase in [
        "Prespecified overall decision category",
        "Proceed to continuous basal-classical axis",
        "TO_VERIFY",
        "Public subtype labels were not used during clustering or K selection",
    ]:
        if phrase not in text:
            fail(f"report missing required phrase: {phrase}")

    print("Phase 4B validation passed.")
    print(f"Validated {len(REQUIRED_TABLES)} tables and {len(REQUIRED_FIGURES)} figures.")


if __name__ == "__main__":
    main()
