#!/usr/bin/env python3
"""Validate Phase 8B host-mechanism outputs."""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "05_results" / "tables"
FIG = ROOT / "05_results" / "figures"
REPORT = ROOT / "04_analysis" / "08_host_microbiome_integration" / "PHASE8B_HOST_MECHANISM_RESULTS.md"

PRIMARY_TAXA = {
    "Azoarcus", "Candida", "Ensifer", "Cutibacterium", "Chryseobacterium",
    "Ochrobactrum", "Burkholderia", "Rhizobium", "Herbaspirillum",
}


def require(ok: bool, msg: str, failures: list[str]) -> None:
    if not ok:
        failures.append(msg)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def count_gzip_rows(path: Path) -> int:
    with gzip.open(path, "rt") as handle:
        return sum(1 for _ in handle)


def main() -> int:
    failures: list[str] = []

    runtime = read_tsv(TABLE / "phase8b_runtime_validation.tsv")
    pathway = read_tsv(TABLE / "phase8b_primary_pathway_associations.tsv")
    tf = read_tsv(TABLE / "phase8b_primary_tf_associations.tsv")
    cov = read_tsv(TABLE / "phase8b_host_covariate_sensitivity.tsv")
    trans = read_tsv(TABLE / "phase8b_transformation_sensitivity.tsv")
    moff = read_tsv(TABLE / "phase8b_moffitt_gene_exclusion_sensitivity.tsv")
    wgcna_soft = read_tsv(TABLE / "phase8b_wgcna_soft_threshold.tsv")
    wgcna_summary = read_tsv(TABLE / "phase8b_wgcna_module_summary.tsv")
    gene_summary = read_tsv(TABLE / "phase8b_host_gene_associations_summary.tsv")
    evidence = read_tsv(TABLE / "phase8b_host_mechanism_evidence.tsv")

    require(runtime["passed"].astype(str).str.upper().eq("TRUE").all(), "runtime validation contains failed checks", failures)
    require(set(pathway["taxon"].unique()) == PRIMARY_TAXA, "pathway table does not contain exactly the nine primary taxa", failures)
    require(set(tf["taxon"].unique()) == PRIMARY_TAXA, "TF table does not contain exactly the nine primary taxa", failures)
    require(pathway.groupby(["taxon", "host_feature_collection"])["bh_q_value"].apply(lambda s: s.notna().all()).all(), "pathway BH q values missing within taxon x collection families", failures)
    require(tf.groupby(["taxon", "host_feature_collection"])["bh_q_value"].apply(lambda s: s.notna().all()).all(), "TF BH q values missing within taxon x collection families", failures)
    require(not cov["covariate"].astype(str).str.contains("estimate_score", case=False, na=False).any(), "ESTIMATE score appears as a covariate", failures)
    require(cov["covariate"].isin(["inferred_tumor_purity", "immune_score", "stromal_score"]).all(), "unexpected or combined TME covariate model present", failures)
    require(trans["rCLR_direction"].notna().any(), "rCLR direction checks are missing", failures)
    require(len(moff) > 0 and {"score_correlation", "direction_consistency"}.issubset(moff.columns), "Moffitt50 exclusion sensitivity missing", failures)
    require((wgcna_soft["SFT.R.sq"].max() >= 0.85), "WGCNA soft-threshold search did not meet R2 >= 0.85", failures)
    require((wgcna_summary["selected_soft_power"] == 5).all(), "WGCNA selected soft power does not match locked rule result", failures)
    require((gene_summary["n_genes"] == 42654).all(), "genome-wide summaries do not use all 42,654 eligible genes", failures)
    require(set(gene_summary["taxon"].unique()) == PRIMARY_TAXA, "gene summary does not contain exactly the nine primary taxa", failures)
    require(set(gene_summary["model"].unique()) == {"primary_CLR", "purity_adjusted", "immune_adjusted", "stromal_adjusted", "rCLR", "exclude_extreme_samples"}, "gene summary model set is incomplete", failures)

    full_dir = TABLE / "phase8b_host_gene_full"
    for taxon in sorted(PRIMARY_TAXA):
        path = full_dir / f"phase8b_host_gene_full_{taxon}.tsv.gz"
        require(path.exists(), f"missing full gene table for {taxon}", failures)
        if path.exists():
            require(count_gzip_rows(path) == 42655, f"full gene table for {taxon} does not have 42,654 result rows plus header", failures)

    allowed = {
        "ROBUST_HOST_MECHANISM", "TRANSFORMATION_SENSITIVE_MECHANISM",
        "COMPOSITION_SENSITIVE_MECHANISM", "SAMPLE_SENSITIVE_MECHANISM",
        "EXPLORATORY_HOST_MECHANISM", "NO_SUPPORTED_MECHANISM", "TO_VERIFY",
    }
    require(set(evidence["evidence_category"].unique()).issubset(allowed), "unexpected evidence category present", failures)
    require(evidence["criteria"].notna().all(), "evidence criteria missing", failures)

    expected_figs = [
        "phase8b_taxon_hallmark_associations.pdf",
        "phase8b_taxon_progeny_associations.pdf",
        "phase8b_taxon_tf_associations.pdf",
        "phase8b_covariate_sensitivity.pdf",
        "phase8b_rclr_direction_sensitivity.pdf",
        "phase8b_wgcna_soft_threshold.pdf",
        "phase8b_wgcna_module_taxon_heatmap.pdf",
        "phase8b_supported_module_annotations.pdf",
        "phase8b_host_gene_enrichment.pdf",
        "phase8b_shared_mechanism_network.pdf",
        "phase8b_mechanism_evidence_summary.pdf",
    ]
    for fig in expected_figs:
        path = FIG / fig
        require(path.exists() and path.stat().st_size > 1000, f"missing or empty figure: {fig}", failures)

    if REPORT.exists():
        text = REPORT.read_text()
        for pat in [r"\bcauses?\b", r"\bcaused\b", r"\bdrives?\b", r"\bdriven\b", r"\bmediates?\b", r"\binduces?\b"]:
            require(re.search(pat, text, flags=re.IGNORECASE) is None, f"causal language found in report: {pat}", failures)

    if failures:
        print("VALIDATION_FAIL")
        for failure in failures:
            print(f"- {failure}")
        return 1
    print("VALIDATION_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
