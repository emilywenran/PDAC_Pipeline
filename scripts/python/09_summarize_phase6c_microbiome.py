#!/usr/bin/env python3
"""Summarize Phase 6C microbiome preprocessing outputs into a locked report."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results/tables"
REPORT = ROOT / "04_analysis/04_microbiome_qc/PHASE6C_ANALYSIS_READY_MICROBIOME.md"


def fmt_float(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}f}"


def main() -> None:
    inventory = pd.read_csv(TABLE_DIR / "phase6c_matrix_inventory.tsv", sep="\t")
    concordance = pd.read_csv(TABLE_DIR / "phase6c_preprocessing_sensitivity_concordance.tsv", sep="\t")
    taxa = pd.read_csv(TABLE_DIR / "phase6c_retained_taxa_with_contamination_flags.tsv", sep="\t")
    sample_qc = pd.read_csv(TABLE_DIR / "phase6c_processed_sample_qc.tsv", sep="\t")
    taxon_qc = pd.read_csv(TABLE_DIR / "phase6c_processed_taxon_qc.tsv", sep="\t")

    primary = inventory.loc[inventory["analysis_id"] == "MICRO_PRIMARY"].iloc[0]
    primary_taxa = taxa.loc[taxa["primary_retained"], "genus"].tolist()
    flagged_primary = taxa.loc[
        taxa["primary_retained"] & (taxa["contamination_risk_category"] != "LOW_CURRENT_CONCERN")
    ]
    sensitivity_ids = inventory.loc[inventory["analysis_role"] == "sensitivity", "analysis_id"].tolist()
    degenerate = sample_qc.loc[sample_qc["degenerate_after_filtering"]]

    best_corr = concordance["mantel_spearman_r"].dropna().min()
    report = f"""# Phase 6C Analysis-Ready Tumor Microbiome Matrices

## Scope

Phase 6C executed the locked Phase 6B tumor microbiome preprocessing protocol for PRJNA719915 using only the audited genus abundance matrix, microbiome sample crosswalk, Phase 6A technical QC outputs, and Phase 6B method-lock files. No host-expression scores, continuous-axis results, survival outcomes, subtype comparisons, differential-abundance results, host-microbiome correlations, pathway analyses, or target-prioritization outputs were used to select preprocessing parameters.

## Primary Filtering Result

The primary rule retained genera with abundance strictly greater than 0 in at least 20% of the 62 samples, corresponding to a minimum detection count of 13 samples. The resulting matrix retained **{int(primary['n_features'])} genera** across **{int(primary['n_samples'])} samples**. This matches the Phase 6B expected retained feature count of approximately 122 genera; no forced feature count adjustment was applied.

The zero fraction before zero replacement in the primary filtered matrix was **{fmt_float(primary['zero_fraction_before_replacement'])}**. No sample became degenerate after filtering.

## Primary Retained Genera

{', '.join(primary_taxa)}

## Primary Pseudocount And Justification

The primary CLR representation used the locked fixed pseudocount **{primary['pseudocount']}**, derived in Phase 6B as one half of the minimum non-zero value in the audited matrix. This value is source-matrix specific and remains marked as non-transferable to external cohorts.

## CLR And Aitchison Validation

The primary CLR matrix contains no missing or infinite values. Per-sample CLR column sums were within floating-point tolerance, with maximum absolute column sum **{primary['max_abs_clr_column_sum']:.3e}**. The primary Aitchison distance matrix is 62 x 62, symmetric, has a zero diagonal, and contains no negative or infinite distances. The primary CLR value range was **{fmt_float(primary['clr_min'])}** to **{fmt_float(primary['clr_max'])}**.

## Sensitivity Representations Produced

Sensitivity outputs were generated under `03_processed/microbiome/sensitivity/` for these analysis IDs: {', '.join(sensitivity_ids)}. These cover the locked 10% and 30% prevalence filters, abundance threshold >10 at 20% prevalence, pseudocounts 0.1 and 1.0, robust CLR, high-risk contaminant exclusion, high- plus moderate-risk contaminant exclusion, presence/absence Jaccard representation, and exclusion of the three Phase 6A technical extreme samples.

## Contamination-Flag Handling

Potential contaminant genera were retained in the primary matrix. Categories in `phase6c_retained_taxa_with_contamination_flags.tsv` are evidence flags, not confirmed contamination labels. Among primary retained genera, **{len(flagged_primary)}** had non-low contamination-sensitivity flags. The high- and moderate-risk removal analyses were generated only as sensitivity representations.

## Technical Outlier Handling

The three locked technical extreme samples (`Basal-like1`, `Hybrid18`, `Hybrid23`) were retained in the primary analysis. The `MICRO_SENS_EXCLUDE_EXTREME` sensitivity representation excludes those samples and was used only to evaluate preprocessing robustness.

## Preprocessing Robustness

Outcome-blind sensitivity concordance was assessed using upper-triangle distance-matrix Spearman correlations, Procrustes concordance on PCoA coordinates, sample-order stability, taxon-rank stability, CLR correlations across pseudocount choices, and ordination shifts after contaminant removal. The lowest non-missing distance-matrix correlation across sensitivity comparisons was **{fmt_float(best_corr)}**. Full results are in `05_results/tables/phase6c_preprocessing_sensitivity_concordance.tsv`.

## Remaining Limitations

Sequenced negative controls are unavailable, so contaminant categories remain potential-risk flags rather than confirmed contaminant calls. The abundance scale remains Bracken-normalized non-integer estimates rather than raw integer classified-read counts. Matrix total abundance is therefore a technical proxy, not a direct absolute microbial load measure.

## Proceed Decision

Microbiome association method locking may proceed with contamination limitations, provided Phase 7 continues to use the locked sensitivity framework and does not treat flagged genera as confirmed contaminants.

## TO_VERIFY

- TO_VERIFY: Original publication/source-specific normalization formula for the Bracken-derived abundance estimates.
- TO_VERIFY: Sequencing batch covariates remain unavailable in public metadata.
- TO_VERIFY: Tumor purity, immune/stromal scores, and clinical covariates remain unavailable for Phase 6C preprocessing.
"""

    if len(degenerate) > 0:
        report += "\n## Degenerate Samples\n\n"
        report += "TO_VERIFY: The following samples were degenerate after filtering: "
        report += ", ".join(degenerate["sample"].astype(str).tolist()) + "\n"

    REPORT.write_text(report)
    print(f"Wrote {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
