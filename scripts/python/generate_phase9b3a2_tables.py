#!/usr/bin/env python3
"""Generate Phase 9B3A.2 Spatial Validation Planning Tables.

This script creates/updates:
- 05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv
- 05_results/tables/phase9b3a2_spatial_model_hierarchy.tsv
- 01_metadata/phase9b3_spatial_dataset_inventory.tsv
- 01_metadata/phase9b3_spatial_parameter_inventory.tsv
- 05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv
- 05_results/tables/phase9b3a_spatial_dataset_qualification.tsv
- 05_results/tables/phase9b3a_spatial_resource_estimate.tsv
"""

import os
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")

def write_tsv(path: Path, header: list[str], rows: list[list[str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(val) for val in row) + "\n")
    print(f"Wrote {len(rows)} rows to {path}")

def generate_spatial_model_hierarchy():
    path = ROOT / "05_results/tables/phase9b3a2_spatial_model_hierarchy.tsv"
    header = [
        "canonical_dataset_id",
        "accession",
        "model_id",
        "model_class",
        "model_formula",
        "random_effect_structure",
        "replicate_unit",
        "inferential_unit",
        "purpose",
        "notes"
    ]
    rows = [
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "HWANG_NAIVE_MODEL_A",
            "Linear mixed-effects model",
            "feature_score ~ compartment + Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction",
            "(1 | patient_id) + (1 | patient_id:ROI_id)",
            "patient",
            "segment",
            "Compartment comparison",
            "Models paired tumor and stroma segments within each ROI. ROI random intercepts serve as paired blocks."
        ],
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "HWANG_NAIVE_MODEL_B",
            "Linear mixed-effects model",
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction",
            "(1 | patient_id)",
            "patient",
            "segment",
            "Tumor-segment-only axis model",
            "Restricted to tumor/malignant segments only to assess axis association."
        ],
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "HWANG_NAIVE_MODEL_C",
            "Linear mixed-effects model",
            "tumor_score_minus_stroma_score ~ Moffitt50_contrast",
            "(1 | patient_id)",
            "patient",
            "ROI",
            "Paired tumor-stroma contrast",
            "Models the paired difference within each ROI as a function of continuous Moffitt50 subtype contrast."
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "HWANG_TREATED_MODEL_A",
            "Linear mixed-effects model",
            "feature_score ~ compartment + Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction",
            "(1 | patient_id) + (1 | patient_id:ROI_id)",
            "patient",
            "segment",
            "Compartment comparison",
            "Models paired tumor and stroma segments in neoadjuvant-treated samples."
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "HWANG_TREATED_MODEL_B",
            "Linear mixed-effects model",
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction",
            "(1 | patient_id)",
            "patient",
            "segment",
            "Tumor-segment-only axis model",
            "Restricted to tumor/malignant segments in treated samples to evaluate treatment remodeling."
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "HWANG_TREATED_MODEL_C",
            "Linear mixed-effects model",
            "tumor_score_minus_stroma_score ~ Moffitt50_contrast",
            "(1 | patient_id)",
            "patient",
            "ROI",
            "Paired tumor-stroma contrast",
            "Models paired differences in neoadjuvant-treated samples."
        ],
        [
            "MONCADA_GSE111672",
            "GSE111672",
            "MONCADA_EXPLORATORY_MODEL",
            "Section-specific spatial analysis",
            "protein_secretion_score ~ Moffitt50_contrast + (1 | patient_id)",
            "None (within-section permutation / patient aggregation)",
            "patient",
            "spot",
            "Exploratory cross-platform spatial consistency",
            "Section-specific permutations, summaries within patients, and direction consistency verification. Not formal replication."
        ]
    ]
    write_tsv(path, header, rows)

def generate_spatial_analysis_unit_and_models():
    path = ROOT / "05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv"
    header = [
        "canonical_dataset_id",
        "accession",
        "platform",
        "patient_count",
        "section_count",
        "spatial_unit_type",
        "spatial_unit_count",
        "biological_replicate",
        "inferential_unit",
        "descriptive_unit",
        "primary_model",
        "reduced_model",
        "random_effect_structure",
        "composition_covariates",
        "treatment_handling",
        "minimum_complete_patients",
        "collinearity_rule",
        "model_failure_rule",
        "notes"
    ]
    rows = [
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "NanoString GeoMx DSP",
            "18",
            "18",
            "paired tumor-stroma segments within ROI",
            "256",
            "patient",
            "segment",
            "segment",
            "Model A: feature_score ~ compartment + Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction + (1 | patient_id) + (1 | patient_id:ROI_id); Model B: score ~ Moffitt50_contrast + CAF + myeloid + lymphoid + (1 | patient_id); Model C: tumor_minus_stroma ~ Moffitt50_contrast + (1 | patient_id)",
            "Omit lymphoid_fraction from Model A/B if unestimable or collinear.",
            "Model A: (1 | patient_id) + (1 | patient_id:ROI_id); Model B/C: (1 | patient_id)",
            "CAF_fraction, myeloid_fraction, lymphoid_fraction",
            "treatment-naïve samples only",
            "10",
            "variance_inflation_factor (VIF) > 10",
            "If model fails to converge, fit reduced model or fall back to patient-level regression and meta-analysis.",
            "ROIs are physically segmented into paired PanCK+ (epithelial) and PanCK- (stroma) areas. ROI random intercepts serve as paired blocks. Section random effects are not identifiable (1 section per patient)."
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "NanoString GeoMx DSP",
            "25",
            "25",
            "paired tumor-stroma segments within ROI",
            "352",
            "patient",
            "segment",
            "segment",
            "Model A: feature_score ~ compartment + Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction + (1 | patient_id) + (1 | patient_id:ROI_id); Model B: score ~ Moffitt50_contrast + CAF + myeloid + lymphoid + (1 | patient_id); Model C: tumor_minus_stroma ~ Moffitt50_contrast + (1 | patient_id)",
            "Omit lymphoid_fraction from Model A/B if unestimable or collinear.",
            "Model A: (1 | patient_id) + (1 | patient_id:ROI_id); Model B/C: (1 | patient_id)",
            "CAF_fraction, myeloid_fraction, lymphoid_fraction",
            "neoadjuvant-treated samples only; analyzed separately to assess treatment-sensitivity effects",
            "10",
            "variance_inflation_factor (VIF) > 10",
            "If model fails to converge, fit reduced model or fall back to patient-level regression and meta-analysis.",
            "Subjected to neoadjuvant therapy. Section random effects are not identifiable."
        ],
        [
            "MONCADA_GSE111672",
            "GSE111672",
            "Microarray Spatial Transcriptomics (ST)",
            "2",
            "6",
            "spot",
            "2248",
            "patient",
            "spot",
            "spot",
            "protein_secretion_score ~ Moffitt50_contrast + (1 | patient_id) (run within-section spatial analyses)",
            "None (within-section permutation / patient aggregation)",
            "None (section-specific and direction consistency analysis)",
            "None",
            "treatment-naïve samples only; classified as exploratory cross-platform spatial consistency",
            "2",
            "variance_inflation_factor (VIF) > 10",
            "Low patient count (n=2) precludes LMM. Perform section-specific spatial analysis and direction consistency check across sections and patients.",
            "Patient A has 4 sections; Patient B has 2 sections. Not treated as formal population-level replication."
        ]
    ]
    write_tsv(path, header, rows)

def generate_spatial_dataset_inventory():
    path = ROOT / "01_metadata/phase9b3_spatial_dataset_inventory.tsv"
    header = [
        "canonical_dataset_id",
        "accession",
        "secondary_accession",
        "repository",
        "publication",
        "PMID_or_DOI",
        "BioProject",
        "cohort_description",
        "patient_count",
        "section_count",
        "ROI_or_spot_count",
        "tumor_region_count",
        "control_or_adjacent_region_count",
        "treatment_group",
        "platform",
        "processed_expression_matrix_available",
        "spatial_coordinate_available",
        "tissue_image_available",
        "histology_annotations_available",
        "matched_single_cell_available",
        "malignant_region_annotations_available",
        "public_download_accessibility",
        "file_formats",
        "approximate_download_size",
        "suitability_status",
        "official_source"
    ]
    rows = [
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "GSE199102",
            "NCBI GEO",
            "Hwang et al. (2022)",
            "10.1038/s41588-023-01411-z",
            "PRJNA826084",
            "Treatment-naïve primary PDAC tumors spatial transcriptomics (18 patients, 18 sections).",
            "18",
            "18",
            "256",
            "256",
            "0",
            "treatment-naïve",
            "NanoString GeoMx DSP",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "RCC, PKC, CSV",
            "~1.2 GB",
            "PRIORITY_1",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199102"
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "GSE199102",
            "NCBI GEO",
            "Hwang et al. (2022)",
            "10.1038/s41588-023-01411-z",
            "PRJNA826084",
            "Neoadjuvant-treated primary PDAC tumors spatial transcriptomics (25 patients, 25 sections).",
            "25",
            "25",
            "352",
            "352",
            "0",
            "neoadjuvant-treated",
            "NanoString GeoMx DSP",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "RCC, PKC, CSV",
            "~1.5 GB",
            "PRIORITY_1",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199102"
        ],
        [
            "MONCADA_GSE111672",
            "GSE111672",
            "PRJNA437847",
            "NCBI GEO",
            "Moncada et al. (2020)",
            "10.1038/s41587-019-0392-8",
            "PRJNA437847",
            "Spatial transcriptomics of primary PDAC tumors from 2 patients (Patient A and Patient B, total 6 sections).",
            "2",
            "6",
            "2248",
            "2248",
            "0",
            "treatment-naïve",
            "Microarray Spatial Transcriptomics (ST)",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "yes",
            "TSV, CSV, JPG",
            "~80 MB",
            "EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672"
        ]
    ]
    write_tsv(path, header, rows)

def generate_spatial_parameter_inventory():
    path = ROOT / "01_metadata/phase9b3_spatial_parameter_inventory.tsv"
    header = [
        "analysis_id",
        "validation_layer",
        "dataset_id",
        "accession",
        "cohort_role",
        "primary_feature",
        "contextual_axis",
        "biological_comparator",
        "secondary_features",
        "spatial_unit",
        "scoring_method",
        "deconvolution_method",
        "statistical_model",
        "negative_control",
        "random_seed",
        "multiple_testing_family",
        "status",
        "notes"
    ]
    rows = [
        [
            "VAL_SPATIAL_CO_LOCAL_HWANG_NAIVE",
            "Layer 3",
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "primary_spatial",
            "HALLMARK_PROTEIN_SECRETION",
            "Moffitt50_contrast",
            "HALLMARK_SPERMATOGENESIS",
            "ELF1, MBD2, ZBTB7A, ZNF384, ZNF740",
            "segment",
            "decoupleR ssGSEA / VIPER (regulon-based)",
            "morphology_segmentation (PanCK+ tumor vs PanCK- stroma)",
            "Model A (Compartment comparison LMM with lymphoid_fraction), Model B (Tumor-only LMM with lymphoid_fraction), Model C (Paired contrast LMM)",
            "1000 label permutations, unrelated hallmark pathways (5), size-matched random gene sets (100)",
            "2026",
            "primary_hypotheses, secondary_localizations, spatial_autocorrelations",
            "ACTIVE_PRIMARY",
            "Primary spatial validation on treatment-naïve samples. Section random effects are not identifiable. Models paired segments."
        ],
        [
            "VAL_SPATIAL_CO_LOCAL_HWANG_TREATED",
            "Layer 3",
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "treatment_sensitivity",
            "HALLMARK_PROTEIN_SECRETION",
            "Moffitt50_contrast",
            "HALLMARK_SPERMATOGENESIS",
            "ELF1, MBD2, ZBTB7A, ZNF384, ZNF740",
            "segment",
            "decoupleR ssGSEA / VIPER (regulon-based)",
            "morphology_segmentation (PanCK+ tumor vs PanCK- stroma)",
            "Model A (Compartment comparison LMM with lymphoid_fraction), Model B (Tumor-only LMM with lymphoid_fraction), Model C (Paired contrast LMM)",
            "1000 label permutations, unrelated hallmark pathways (5), size-matched random gene sets (100)",
            "2026",
            "primary_hypotheses, secondary_localizations, spatial_autocorrelations",
            "ACTIVE_TREATMENT_SENSITIVITY",
            "Treatment sensitivity spatial validation. Section random effects are not identifiable. Models paired segments."
        ],
        [
            "VAL_SPATIAL_CO_LOCAL_MONCADA",
            "Layer 3",
            "MONCADA_GSE111672",
            "GSE111672",
            "exploratory_cross_platform_spatial_consistency",
            "HALLMARK_PROTEIN_SECRETION",
            "Moffitt50_contrast",
            "HALLMARK_SPERMATOGENESIS",
            "ELF1, MBD2, ZBTB7A, ZNF384, ZNF740",
            "spot",
            "decoupleR ssGSEA / VIPER (regulon-based)",
            "reference-based deconvolution (RCTD/MIA using matched scRNA-seq)",
            "Section-specific spatial analysis, within-section permutations, and direction consistency verification.",
            "1000 spot-coordinate permutations, unrelated hallmark pathways (5), size-matched random gene sets (100)",
            "2026",
            "primary_hypotheses, secondary_localizations, spatial_autocorrelations",
            "ACTIVE_EXPLORATORY",
            "Exploratory cross-platform spatial consistency. Low patient count (n=2) precludes formal replication."
        ]
    ]
    write_tsv(path, header, rows)

def generate_authoritative_cohort_set():
    path = ROOT / "05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv"
    header = [
        "canonical_dataset_id",
        "accession",
        "spatial_suitability_status",
        "spatial_analysis_role",
        "included_in_phase9b3_primary",
        "included_in_phase9b3_secondary",
        "treatment_sensitivity_only",
        "current_execution_authorized"
    ]
    rows = [
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "PRIORITY_1",
            "PRIMARY",
            "TRUE",
            "FALSE",
            "FALSE",
            "TRUE"
        ],
        [
            "MONCADA_GSE111672",
            "GSE111672",
            "EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY",
            "EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY",
            "FALSE",
            "FALSE",
            "FALSE",
            "TRUE"
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "PRIORITY_1",
            "TREATMENT_SENSITIVITY",
            "FALSE",
            "FALSE",
            "TRUE",
            "TRUE"
        ],
        [
            "GSE274103",
            "GSE274103",
            "EXPLORATORY_ONLY",
            "EXPLORATORY_ONLY",
            "FALSE",
            "FALSE",
            "FALSE",
            "FALSE"
        ],
        [
            "GSE272362",
            "GSE272362",
            "EXPLORATORY_ONLY",
            "EXPLORATORY_ONLY",
            "FALSE",
            "FALSE",
            "FALSE",
            "FALSE"
        ]
    ]
    write_tsv(path, header, rows)

def generate_dataset_qualification():
    path = ROOT / "05_results/tables/phase9b3a_spatial_dataset_qualification.tsv"
    header = [
        "canonical_dataset_id",
        "accession",
        "spatial_platform",
        "PDAC_patients",
        "tissue_sections",
        "spatial_suitability_status",
        "suitability_rationale",
        "exclusion_rules_applied"
    ]
    rows = [
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "NanoString GeoMx DSP",
            "18",
            "18",
            "PRIORITY_1",
            "Large cohort of 18 treatment-naïve patients with whole-transcriptome digital spatial profiling, matched snRNA-seq, and tumor/stroma morphology segments.",
            "none"
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "NanoString GeoMx DSP",
            "25",
            "25",
            "PRIORITY_1",
            "Large cohort of 25 neoadjuvant-treated patients; suitable as a treatment-sensitivity validation layer.",
            "none"
        ],
        [
            "MONCADA_GSE111672",
            "GSE111672",
            "Microarray Spatial Transcriptomics (ST)",
            "2",
            "6",
            "EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY",
            "Well-characterized matched scRNA-seq and spatial transcriptomics sections, but low sample size (n=2 patients) and lower platform resolution than Visium. Reclassified as exploratory.",
            "low_patient_count_precludes_formal_replication"
        ],
        [
            "GSE274103",
            "GSE274103",
            "10x Visium FFPE",
            "5",
            "5",
            "EXPLORATORY_ONLY",
            "Small cohort (n=5 patients) with 10x Visium profiling. Excluded from primary execution set to maintain focus on the two main cohorts.",
            "non_authorized_for_phase9b3_primary_execution"
        ],
        [
            "GSE272362",
            "GSE272362",
            "10x Visium",
            "10",
            "13",
            "EXPLORATORY_ONLY",
            "Secondary validation dataset, includes metastatic tissues. Excluded from primary execution set.",
            "non_authorized_for_phase9b3_primary_execution"
        ]
    ]
    write_tsv(path, header, rows)

def generate_resource_estimate():
    path = ROOT / "05_results/tables/phase9b3a_spatial_resource_estimate.tsv"
    header = [
        "canonical_dataset_id",
        "accession",
        "data_type",
        "patient_count",
        "section_count",
        "expected_spatial_units",
        "estimated_download_size",
        "disk_requirement",
        "peak_ram",
        "expected_runtime",
        "required_packages",
        "macbook_practical",
        "hpc_recommended",
        "notes"
    ]
    rows = [
        [
            "HWANG_GSE202051_NAIVE",
            "GSE202051",
            "GeoMx DSP WTA",
            "18",
            "18",
            "256 segments (paired tumor-stroma)",
            "~1.2 GB",
            "~5 GB",
            "16 GB",
            "~20 minutes",
            "R (decoupleR, GSVA, limma, StandR), Python (pandas, statsmodels)",
            "yes",
            "no",
            "ROIs are selected regions. Processing is relatively lightweight because data is not at single-cell density."
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "GeoMx DSP WTA",
            "25",
            "25",
            "352 segments (paired tumor-stroma)",
            "~1.5 GB",
            "~6 GB",
            "16 GB",
            "~25 minutes",
            "R (decoupleR, GSVA, limma, StandR), Python (pandas, statsmodels)",
            "yes",
            "no",
            "Treatment-sensitivity cohort."
        ],
        [
            "MONCADA_GSE111672",
            "GSE111672",
            "Microarray Spatial",
            "2",
            "6",
            "2248 spots",
            "~80 MB",
            "~500 MB",
            "8 GB",
            "~10 minutes",
            "R (decoupleR, GSVA, Seurat), Python (pandas, statsmodels)",
            "yes",
            "no",
            "First-generation ST data. Lightweight."
        ],
        [
            "GSE274103",
            "GSE274103",
            "10x Visium FFPE",
            "5",
            "5",
            "~20000 spots",
            "~150 MB",
            "~1 GB",
            "16 GB",
            "~15 minutes",
            "R (Seurat, decoupleR, GSVA)",
            "yes",
            "no",
            "Not authorized for current execution."
        ],
        [
            "GSE272362",
            "GSE272362",
            "10x Visium",
            "10",
            "13",
            "~40000 spots",
            "~500 MB",
            "~2 GB",
            "16 GB",
            "~20 minutes",
            "R (Seurat, decoupleR, GSVA)",
            "yes",
            "no",
            "Not authorized for current execution."
        ]
    ]
    write_tsv(path, header, rows)

def main():
    print("Generating Phase 9B3A.2 spatial validation planning tables...")
    generate_spatial_model_hierarchy()
    generate_spatial_analysis_unit_and_models()
    generate_spatial_dataset_inventory()
    generate_spatial_parameter_inventory()
    generate_authoritative_cohort_set()
    generate_dataset_qualification()
    generate_resource_estimate()
    print("All tables successfully generated/updated.")

if __name__ == "__main__":
    main()
