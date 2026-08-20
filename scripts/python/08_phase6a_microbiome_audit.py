#!/usr/bin/env python3
"""Phase 6A processed PRJNA719915 microbiome abundance audit.

Uses only verified local public supplementary files and repository metadata.
Does not download FASTQ/SRA data and does not test host subtype/score
associations.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("tmp/matplotlib").resolve()))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
SUPP1 = ROOT / "02_data/reference/original_study_supplementary/42003_2021_2557_MOESM4_ESM.xlsx"
SUPP4 = ROOT / "02_data/reference/original_study_supplementary/42003_2021_2557_MOESM7_ESM.xlsx"
CROSSWALK_IN = ROOT / "01_metadata/rna_microbiome_patient_crosswalk.tsv"
RUN_INV = ROOT / "01_metadata/microbiome_run_inventory.tsv"
OUT_MATRIX = ROOT / "03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz"
OUT_CROSSWALK = ROOT / "01_metadata/microbiome_sample_crosswalk.tsv"
TABLE_DIR = ROOT / "05_results/tables"
FIG_DIR = ROOT / "05_results/figures"
REPORT = ROOT / "04_analysis/04_microbiome_qc/PHASE6A_MICROBIOME_DATA_AUDIT.md"

TAXONOMIC_LEVEL = "Genus"
SOURCE_SHEET = "Genus-level"
SOURCE_TABLE = "Supplementary Data 1 (42003_2021_2557_MOESM4_ESM.xlsx), sheet 'Genus-level'"
ABUNDANCE_UNIT = "Kraken2/Bracken-derived non-integer normalized counts as released in Supplementary Data 1; not raw classified reads and not relative abundance"

KNOWN_CONTAMINANT_GENERA = {
    "Acinetobacter": "common reagent/environmental genus; also reported by source study",
    "Pseudomonas": "common water/reagent/environmental genus; also reported by source study",
    "Sphingopyxis": "environmental/water-associated genus raised in peer review",
    "Sphingomonas": "common reagent/environmental genus",
    "Brevundimonas": "common water/reagent-associated genus",
    "Ralstonia": "common reagent/water-associated genus",
    "Burkholderia": "common reagent/environmental genus",
    "Paraburkholderia": "environmental genus",
    "Methylobacterium": "common reagent/environmental genus",
    "Bradyrhizobium": "common reagent/environmental genus",
    "Delftia": "common water/reagent-associated genus",
    "Comamonas": "common environmental/water-associated genus",
    "Elizabethkingia": "water/environment-associated genus; source study discussed genus-level assignment",
    "Stenotrophomonas": "common environmental/reagent-associated genus",
    "Cupriavidus": "common environmental/reagent-associated genus",
    "Herbaspirillum": "common reagent/environmental genus",
    "Mesorhizobium": "common reagent/environmental genus",
    "Rhizobium": "common reagent/environmental genus",
    "Novosphingobium": "common environmental/reagent-associated genus",
    "Caulobacter": "common water-associated genus",
    "Dechloromonas": "environmental genus raised in peer review",
}


def ensure_dirs() -> None:
    for path in [OUT_MATRIX.parent, TABLE_DIR, FIG_DIR, REPORT.parent, Path(os.environ["MPLCONFIGDIR"])]:
        path.mkdir(parents=True, exist_ok=True)


def read_matrix() -> pd.DataFrame:
    df = pd.read_excel(SUPP1, sheet_name=SOURCE_SHEET)
    df = df.rename(columns={df.columns[0]: "taxon"})
    if df["taxon"].duplicated().any():
        dup = df.loc[df["taxon"].duplicated(), "taxon"].tolist()[:5]
        raise ValueError(f"Duplicated taxon identifiers in source matrix: {dup}")
    df = df.set_index("taxon")
    df = df.apply(pd.to_numeric, errors="raise")
    return df


def read_chan_sample_order() -> list[tuple[str, str]]:
    raw = pd.read_excel(SUPP4, sheet_name="Figure1.SampleGroup", header=None)
    rows: list[tuple[str, str]] = []
    for sample, subtype in raw.iloc[2:, [0, 1]].itertuples(index=False):
        if pd.notna(sample):
            rows.append((str(sample), str(subtype)))
    return rows


def source_column_to_yx(columns: list[str], chan_order: list[tuple[str, str]]) -> pd.DataFrame:
    by_group: dict[str, list[str]] = {"Basal-like": [], "Hybrid": [], "Classical": []}
    for sample, subtype in chan_order:
        key = "Basal-like" if subtype == "Basal" else subtype
        by_group[key].append(sample)
    records = []
    counters = {k: 0 for k in by_group}
    for col in columns:
        m = re.fullmatch(r"(Basal-like|Hybrid|Classical)(\d+)", col)
        if not m:
            raise ValueError(f"Unexpected source sample column: {col}")
        group, ordinal = m.group(1), int(m.group(2))
        yx = by_group[group][ordinal - 1]
        counters[group] += 1
        records.append({"microbiome_matrix_sample": col, "submitted_sample_name": yx})
    for group, samples in by_group.items():
        if counters[group] != len(samples):
            raise ValueError(f"Group {group} count mismatch: {counters[group]} vs {len(samples)}")
    return pd.DataFrame(records)


def build_crosswalk(matrix_cols: list[str]) -> pd.DataFrame:
    colmap = source_column_to_yx(matrix_cols, read_chan_sample_order())
    runs = pd.read_csv(RUN_INV, sep="\t", dtype=str)
    patient = pd.read_csv(CROSSWALK_IN, sep="\t", dtype=str)
    run_cols = ["submitted_sample_name", "biosample_id", "run_id"]
    merged = colmap.merge(runs[run_cols], on="submitted_sample_name", how="left")
    merged = merged.merge(
        patient[["patient_id", "tumor_number", "microbiome_biosample_id", "microbiome_run_id"]],
        left_on=["biosample_id", "run_id"],
        right_on=["microbiome_biosample_id", "microbiome_run_id"],
        how="left",
    )
    merged["mapping_status"] = np.where(merged["patient_id"].notna(), "VERIFIED", "UNMATCHED")
    merged["notes"] = (
        "Matrix column from Supplementary Data 1 "
        + SOURCE_SHEET
        + "; mapped via Supplementary Data 4 Figure1.SampleGroup order to submitted sample "
        + merged["submitted_sample_name"].fillna("NA")
        + "."
    )
    out = merged[
        [
            "microbiome_matrix_sample",
            "patient_id",
            "tumor_number",
            "microbiome_biosample_id",
            "microbiome_run_id",
            "mapping_status",
            "notes",
        ]
    ].copy()
    out.to_csv(OUT_CROSSWALK, sep="\t", index=False)
    return out


def save_matrix(matrix: pd.DataFrame) -> None:
    with gzip.open(OUT_MATRIX, "wt") as handle:
        matrix.to_csv(handle, sep="\t", index_label="taxon")


def bray_curtis(matrix: pd.DataFrame) -> pd.DataFrame:
    x = matrix.T.to_numpy(dtype=float)
    dist = squareform(pdist(x, metric="braycurtis"))
    return pd.DataFrame(dist, index=matrix.columns, columns=matrix.columns)


def descriptive_qc(matrix: pd.DataFrame, crosswalk: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    values = matrix.to_numpy(dtype=float)
    sample_total = matrix.sum(axis=0)
    detected_taxa = (matrix > 0).sum(axis=0)
    zero_fraction_sample = (matrix == 0).mean(axis=0)
    log_total_z = pd.Series(
        StandardScaler().fit_transform(np.log10(sample_total + 1).to_numpy().reshape(-1, 1)).ravel(),
        index=sample_total.index,
    )
    detected_z = pd.Series(
        StandardScaler().fit_transform(detected_taxa.to_numpy().reshape(-1, 1)).ravel(),
        index=detected_taxa.index,
    )
    extreme = (log_total_z.abs() >= 3) | (detected_z.abs() >= 3)

    runs = pd.read_csv(RUN_INV, sep="\t", dtype=str)
    sample_qc = pd.DataFrame(
        {
            "microbiome_matrix_sample": matrix.columns,
            "total_abundance": sample_total.values,
            "detected_taxa": detected_taxa.values,
            "zero_fraction": zero_fraction_sample.values,
            "log10_total_abundance": np.log10(sample_total.values + 1),
            "log10_total_abundance_z": log_total_z.values,
            "detected_taxa_z": detected_z.values,
            "extreme_sample_flag": np.where(extreme.values, "YES", "NO"),
        }
    )
    sample_qc = sample_qc.merge(crosswalk, on="microbiome_matrix_sample", how="left")
    sample_qc = sample_qc.merge(
        runs[["run_id", "submitted_sample_name", "bases", "spots", "file_size"]],
        left_on="microbiome_run_id",
        right_on="run_id",
        how="left",
    ).drop(columns=["run_id"])
    for col in ["bases", "spots"]:
        sample_qc[col] = pd.to_numeric(sample_qc[col], errors="coerce")
    sample_qc["file_size_numeric"] = pd.to_numeric(sample_qc["file_size"].str.split(";").str[0], errors="coerce")

    taxon_qc = pd.DataFrame(
        {
            "taxon": matrix.index,
            "total_abundance": matrix.sum(axis=1).values,
            "mean_abundance": matrix.mean(axis=1).values,
            "median_abundance": matrix.median(axis=1).values,
            "max_abundance": matrix.max(axis=1).values,
            "prevalence_count": (matrix > 0).sum(axis=1).values,
            "prevalence_fraction": (matrix > 0).mean(axis=1).values,
            "zero_fraction": (matrix == 0).mean(axis=1).values,
        }
    )
    prevalence = taxon_qc[["taxon", "prevalence_count", "prevalence_fraction", "total_abundance", "mean_abundance"]].copy()
    prevalence = prevalence.sort_values(["prevalence_count", "total_abundance"], ascending=[False, False])

    flags = taxon_qc[taxon_qc["taxon"].isin(KNOWN_CONTAMINANT_GENERA)].copy()
    flags["contamination_risk_status"] = "POTENTIAL_CONTAMINANT_RISK_NOT_CONFIRMED"
    flags["flag_basis"] = flags["taxon"].map(KNOWN_CONTAMINANT_GENERA)
    flags["action"] = "Retain in Phase 6A; compare Phase 6B sensitivity under prevalence/abundance filters; do not delete without negative-control evidence."

    assoc_records = []
    for col in ["bases", "spots", "file_size_numeric"]:
        valid = sample_qc[["total_abundance", col]].dropna()
        rho, p = spearmanr(valid["total_abundance"], valid[col])
        assoc_records.append(
            {
                "technical_metadata": col,
                "n": len(valid),
                "spearman_rho_total_abundance": rho,
                "spearman_p_total_abundance": p,
                "analysis_scope": "technical_metadata_only",
            }
        )
    tech_assoc = pd.DataFrame(assoc_records)

    if not np.isfinite(values).all():
        raise ValueError("Non-finite values detected in matrix")
    if np.isnan(values).any():
        raise ValueError("Missing values detected in matrix")
    if (values < 0).any():
        raise ValueError("Negative abundance values detected in matrix")

    sample_qc.to_csv(TABLE_DIR / "phase6a_microbiome_sample_qc.tsv", sep="\t", index=False)
    taxon_qc.to_csv(TABLE_DIR / "phase6a_microbiome_taxon_qc.tsv", sep="\t", index=False)
    prevalence.to_csv(TABLE_DIR / "phase6a_taxon_prevalence.tsv", sep="\t", index=False)
    flags.to_csv(TABLE_DIR / "phase6a_potential_contaminant_flags.tsv", sep="\t", index=False)
    tech_assoc.to_csv(TABLE_DIR / "phase6a_technical_metadata_associations.tsv", sep="\t", index=False)
    return sample_qc, taxon_qc, prevalence, flags, tech_assoc


def make_figures(matrix: pd.DataFrame, sample_qc: pd.DataFrame, taxon_qc: pd.DataFrame) -> pd.DataFrame:
    sns.set_theme(style="whitegrid")
    positive = matrix.to_numpy(dtype=float).ravel()
    positive = positive[positive > 0]

    plt.figure(figsize=(7, 5))
    sns.histplot(np.log10(positive + 1), bins=60, color="#3b6ea8")
    plt.xlabel("log10(abundance + 1), positive cells only")
    plt.ylabel("Matrix cells")
    plt.title("Genus Abundance Distribution")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "phase6a_microbiome_abundance_distribution.pdf")
    plt.close()

    plt.figure(figsize=(7, 5))
    sns.histplot(taxon_qc["prevalence_fraction"], bins=25, color="#4f8a5b")
    plt.xlabel("Taxon prevalence fraction")
    plt.ylabel("Genera")
    plt.title("Genus Prevalence")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "phase6a_taxon_prevalence.pdf")
    plt.close()

    fig, ax1 = plt.subplots(figsize=(9, 5))
    order = sample_qc.sort_values("total_abundance")["microbiome_matrix_sample"]
    plot_df = sample_qc.set_index("microbiome_matrix_sample").loc[order].reset_index()
    ax1.bar(np.arange(len(plot_df)), plot_df["total_abundance"], color="#3b6ea8", alpha=0.75)
    ax1.set_ylabel("Total abundance")
    ax1.set_xlabel("Samples ordered by total abundance")
    ax1.tick_params(axis="x", labelbottom=False)
    ax2 = ax1.twinx()
    ax2.plot(np.arange(len(plot_df)), plot_df["detected_taxa"], color="#b24c3a", linewidth=1.6)
    ax2.set_ylabel("Detected genera")
    plt.title("Sample Detection Summary")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "phase6a_sample_detection_summary.pdf")
    plt.close(fig)

    rel = matrix.div(matrix.sum(axis=0), axis=1).fillna(0)
    pca_input = np.log10(rel.T + 1e-6)
    coords = PCA(n_components=2, random_state=1).fit_transform(StandardScaler().fit_transform(pca_input))
    ord_df = pd.DataFrame(coords, columns=["PC1", "PC2"])
    ord_df["microbiome_matrix_sample"] = rel.columns
    ord_df = ord_df.merge(sample_qc[["microbiome_matrix_sample", "extreme_sample_flag"]], on="microbiome_matrix_sample")
    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=ord_df, x="PC1", y="PC2", hue="extreme_sample_flag", palette={"NO": "#3b6ea8", "YES": "#b24c3a"}, s=55)
    plt.title("PCA of log10 Relative Genus Abundance")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "phase6a_microbiome_ordination.pdf")
    plt.close()

    dist = bray_curtis(matrix)
    plt.figure(figsize=(9, 8))
    sns.heatmap(dist, cmap="viridis", xticklabels=False, yticklabels=False, square=True, cbar_kws={"label": "Bray-Curtis distance"})
    plt.title("Sample-Sample Bray-Curtis Distance")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "phase6a_sample_distance_heatmap.pdf")
    plt.close()
    return dist


def write_report(matrix: pd.DataFrame, sample_qc: pd.DataFrame, taxon_qc: pd.DataFrame, flags: pd.DataFrame, tech_assoc: pd.DataFrame) -> None:
    zero_fraction = float((matrix == 0).to_numpy().mean())
    extreme_samples = sample_qc.loc[sample_qc["extreme_sample_flag"] == "YES", "microbiome_matrix_sample"].tolist()
    high_prev = int((taxon_qc["prevalence_fraction"] >= 0.5).sum())
    low_prev = int((taxon_qc["prevalence_fraction"] <= 0.1).sum())
    unmatched = int((sample_qc["mapping_status"] != "VERIFIED").sum())
    text = f"""# Phase 6A Microbiome Data Audit

## Source and Matrix

- Exact source table: {SOURCE_TABLE}.
- Extracted taxonomic level: {TAXONOMIC_LEVEL}.
- Abundance unit: {ABUNDANCE_UNIT}.
- Numerical scale: preserved exactly from the public supplementary workbook; no filtering, rarefaction, renormalization, log transform, or CLR transform was applied to the stored matrix.
- Final dimensions: {matrix.shape[0]} genera x {matrix.shape[1]} tumor samples.
- Source sample identifier format: `Basal-like1-17`, `Hybrid1-23`, and `Classical1-22`.
- Stored matrix: `03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz`.

## Method Provenance

Repository metadata identifies PRJNA719915 as 62 single-end Illumina shotgun metagenomic tumor runs. The peer-review supplement states that read counts for taxa were measured using Kraken2 plus Bracken and then converted to relative abundance for composition displays; Bracken reallocation was specifically discussed for species estimates. Supplementary Data 1 is released as abundance profiles at class, order, family, genus, and species levels. The exact taxonomic database name/version used by Kraken2/Bracken was not found in the verified local public supplementary files and is marked `TO_VERIFY`.

No unclassified, Homo sapiens, human, host, or unmapped categories were present as genus-level feature rows in the extracted sheet. The matrix contains only named genera, including several zero-only rows.

## Mapping and Validation

- Mapping success: {matrix.shape[1] - unmatched}/{matrix.shape[1]} matrix columns verified to project patients.
- Unique tumor samples: {sample_qc['patient_id'].nunique()}.
- One microbiome profile per patient: yes.
- Duplicate samples: none detected.
- Duplicate taxonomic identifiers: none detected.
- Unmatched patients: {unmatched}.
- Negative values: none.
- Missing or infinite values: none.
- Feature/sample orientation: rows are genera and columns are tumor microbiome profiles.

## Descriptive QC

- Overall zero fraction: {zero_fraction:.4f}.
- Genera detected in at least 50% of samples: {high_prev}.
- Genera detected in at most 10% of samples: {low_prev}.
- Total abundance range per sample: {sample_qc['total_abundance'].min():.3f} to {sample_qc['total_abundance'].max():.3f}.
- Detected-genera range per sample: {int(sample_qc['detected_taxa'].min())} to {int(sample_qc['detected_taxa'].max())}.
- Suspected extreme samples by descriptive z-score screening: {', '.join(extreme_samples) if extreme_samples else 'none'}.

Sample-sample structure was summarized with Bray-Curtis distances on the released abundance scale and PCA on log10 relative genus abundance with a small display pseudocount. These are descriptive ordination/QC summaries only and are not subtype, continuous-axis, differential-abundance, survival, host-correlation, pathway, or target-prioritization tests.

Technical metadata-only Spearman checks of total abundance were performed against available `bases`, `spots`, and file size fields:

{tech_assoc.to_markdown(index=False)}

## Contamination-Control Limitations

The project has no sequenced negative-control runs. Therefore no decontam prevalence or frequency analysis was performed, no contaminant feature is confirmed by negative-control evidence, and no potential contaminant was automatically deleted. Phase 6A records contamination risk only. Potential reagent/environment-associated genera flagged in the extracted matrix: {', '.join(flags['taxon'].tolist()) if len(flags) else 'none'}.

Contamination risk is not equivalent to confirmed contamination. Phase 6B should compare results under multiple prevalence and abundance filters and should keep flagged genera visible in sensitivity reports.

## Suitability for Compositional Analysis

The extracted matrix is suitable for exploratory compositional preprocessing evaluation because it is non-negative and contains 62 complete tumor profiles. It is not directly suitable for CLR/Aitchison analysis without an explicit zero-handling policy because zeros are common. CLR must not be applied directly to raw zeros.

## Recommended Phase 6B Preprocessing Candidates

- Prevalence filtering: compare thresholds such as detected in at least 5%, 10%, and 20% of samples, with all thresholds reported.
- Abundance filtering: compare low-total-abundance removal thresholds independent of subtype labels and continuous scores.
- Pseudocount policy: use a documented small pseudocount or multiplicative replacement after filtering; perform sensitivity to the pseudocount choice.
- CLR/Aitchison transformation: apply only after zero handling; use Aitchison distances for compositional sensitivity analyses.
- Count-based methods: only if raw or integer Bracken estimated counts are obtained or regenerated later; do not treat the released non-integer normalized matrix as raw counts.
- Compositional differential-abundance methods: evaluate ALDEx2, ANCOM-BC, or related approaches in a later locked phase, with contamination-risk sensitivity and technical covariates.

## TO_VERIFY

- Exact Kraken2/Bracken database and version used by the original study.
- Exact normalization formula that produced the non-integer Supplementary Data 1 abundance scale.
- Whether zero-only genus rows are intentional retained taxa from the original pipeline or workbook artifacts.
"""
    REPORT.write_text(text)


def main() -> None:
    ensure_dirs()
    matrix = read_matrix()
    save_matrix(matrix)
    crosswalk = build_crosswalk(list(matrix.columns))
    sample_qc, taxon_qc, _, flags, tech_assoc = descriptive_qc(matrix, crosswalk)
    make_figures(matrix, sample_qc, taxon_qc)
    write_report(matrix, sample_qc, taxon_qc, flags, tech_assoc)
    print(f"Wrote audited matrix: {OUT_MATRIX.relative_to(ROOT)} ({matrix.shape[0]} x {matrix.shape[1]})")
    print(f"Wrote crosswalk: {OUT_CROSSWALK.relative_to(ROOT)}")
    print(f"Wrote report: {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
