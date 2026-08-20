#!/usr/bin/env python3
"""Generate Phase 9B3A Spatial Transcriptomic Validation Planning Tables.

This script creates:
- 05_results/tables/phase9b3a_spatial_feature_hierarchy.tsv
- 05_results/tables/phase9b3a_spatial_dataset_qualification.tsv
- 05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv
- 05_results/tables/phase9b3a_spatial_resource_estimate.tsv
- 01_metadata/phase9b3_spatial_dataset_inventory.tsv
- 01_metadata/phase9b3_spatial_parameter_inventory.tsv
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

def generate_feature_hierarchy():
    path = ROOT / "05_results/tables/phase9b3a_spatial_feature_hierarchy.tsv"
    header = [
        "feature_layer",
        "feature_name",
        "bulk_evidence_category",
        "single_cell_evidence_category",
        "spatial_analysis_role",
        "primary_or_secondary",
        "spatial_coverage_required",
        "formal_inference_allowed",
        "exclusion_rule",
        "notes"
    ]
    rows = [
        # Primary Feature
        [
            "Hallmark",
            "HALLMARK_PROTEIN_SECRETION",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "MALIGNANT_CELL_INTRINSIC_SUPPORT",
            "PRIMARY_TARGET",
            "PRIMARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Primary spatial validation target. Malignant cell intrinsic support replicated in corrected Phase 9B2R."
        ],
        # Primary Contextual Axis
        [
            "Subtype_Axis",
            "Moffitt50_contrast",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "REPLICATED_MALIGNANT_AXIS",
            "PRIMARY_CONTEXTUAL_AXIS",
            "PRIMARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Primary contextual axis for continuous host transcriptomic state mapping."
        ],
        # Prespecified Biological Comparator
        [
            "Hallmark",
            "HALLMARK_SPERMATOGENESIS",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "PRESPECIFIED_BIOLOGICAL_COMPARATOR",
            "PRIMARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Prespecified biological comparator; negative control pathway."
        ],
        # Secondary Localization Features
        [
            "TF_regulon",
            "ELF1",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "STROMAL_OR_IMMUNE_SOURCE_SUPPORTED",
            "SECONDARY_LOCALIZATION_TARGET",
            "SECONDARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Secondary localization target. Localized to stromal/immune compartment."
        ],
        [
            "TF_regulon",
            "MBD2",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "STROMAL_OR_IMMUNE_SOURCE_SUPPORTED",
            "SECONDARY_LOCALIZATION_TARGET",
            "SECONDARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Secondary localization target. Localized to stromal/immune compartment."
        ],
        [
            "TF_regulon",
            "ZBTB7A",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "STROMAL_OR_IMMUNE_SOURCE_SUPPORTED",
            "SECONDARY_LOCALIZATION_TARGET",
            "SECONDARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Secondary localization target. Localized to stromal/immune compartment."
        ],
        [
            "TF_regulon",
            "ZNF384",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "STROMAL_OR_IMMUNE_SOURCE_SUPPORTED",
            "SECONDARY_LOCALIZATION_TARGET",
            "SECONDARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Secondary localization target. Localized to stromal/immune compartment."
        ],
        [
            "TF_regulon",
            "ZNF740",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "PARTIAL_CELLULAR_SUPPORT",
            "SECONDARY_LOCALIZATION_TARGET",
            "SECONDARY",
            "0.80",
            "TRUE",
            "coverage_below_80_percent",
            "Secondary localization target. Exhibits partial cellular compartment support."
        ],
        # WGCNA Modules (Sensitivity only, not authorized for formal inference)
        [
            "WGCNA_module",
            "MEblack",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "INSUFFICIENT_SINGLE_CELL_DATA",
            "SPATIAL_SENSITIVITY_LAYER",
            "SECONDARY",
            "0.80",
            "FALSE",
            "coverage_below_80_percent",
            "WGCNA module. Excluded in SC. Spatial inference blocked if spatial dataset-specific coverage < 80%."
        ],
        [
            "WGCNA_module",
            "MEblue",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "INSUFFICIENT_SINGLE_CELL_DATA",
            "SPATIAL_SENSITIVITY_LAYER",
            "SECONDARY",
            "0.80",
            "FALSE",
            "coverage_below_80_percent",
            "WGCNA module. Excluded in SC. Spatial inference blocked if spatial dataset-specific coverage < 80%."
        ],
        [
            "WGCNA_module",
            "MEgreen",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "INSUFFICIENT_SINGLE_CELL_DATA",
            "SPATIAL_SENSITIVITY_LAYER",
            "SECONDARY",
            "0.80",
            "FALSE",
            "coverage_below_80_percent",
            "WGCNA module. Excluded in SC. Spatial inference blocked if spatial dataset-specific coverage < 80%."
        ],
        [
            "WGCNA_module",
            "MEtan",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "INSUFFICIENT_SINGLE_CELL_DATA",
            "SPATIAL_SENSITIVITY_LAYER",
            "SECONDARY",
            "0.80",
            "FALSE",
            "coverage_below_80_percent",
            "WGCNA module. Excluded in SC. Spatial inference blocked if spatial dataset-specific coverage < 80%."
        ],
        [
            "WGCNA_module",
            "MEgreenyellow",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "INSUFFICIENT_SINGLE_CELL_DATA",
            "SPATIAL_SENSITIVITY_LAYER",
            "SECONDARY",
            "0.80",
            "FALSE",
            "coverage_below_80_percent",
            "WGCNA module. Excluded in SC. Spatial inference blocked if spatial dataset-specific coverage < 80%."
        ],
        # Composition-Context Features (CELL_COMPOSITION_EXPLAINED in single-cell)
        [
            "TF_regulon",
            "BHLHE40",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "CTCFL",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "E2F6",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "GRHL2",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "IRF3",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "JUNB",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "KLF1",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "KLF9",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "MBD1",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "MNT",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "MXI1",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "OTX2",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "SIX5",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "SNAI2",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "SNAPC4",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "TFAP4",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "TP63",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "ZBED1",
            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        [
            "TF_regulon",
            "ZNF24",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "CELL_COMPOSITION_EXPLAINED",
            "COMPOSITION_CONTEXT_EVALUATION",
            "SECONDARY",
            "0.80",
            "FALSE",
            "cell_composition_explained_not_mechanistic",
            "Mapped only for cell-type spatial distribution and consistency checks."
        ],
        # Un-supported TF
        [
            "TF_regulon",
            "KLF13",
            "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "NOT_SUPPORTED_AT_CELLULAR_LEVEL",
            "NOT_SUPPORTED_SPATIALLY",
            "SECONDARY",
            "0.80",
            "FALSE",
            "not_supported_at_cellular_level",
            "Not supported at cellular level in SC validation. Mapped only as negative biological baseline."
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
            "PRIORITY_2",
            "Well-characterized matched scRNA-seq and spatial transcriptomics sections, but low sample size (n=2 patients) and lower platform resolution than Visium.",
            "none"
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
            "PRIORITY_2",
            "SECONDARY",
            "FALSE",
            "TRUE",
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
            "256 segments",
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
            "352 segments",
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
            "PRIORITY_2",
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
            "LMM: Feature_Score ~ Moffitt50_contrast + CAF_Fraction + Myeloid_Fraction + Lymphoid_Fraction + (1 | patient_id)",
            "1000 label permutations, unrelated hallmark pathways (5), size-matched random gene sets (100)",
            "2026",
            "primary_hypotheses, secondary_localizations, spatial_autocorrelations",
            "ACTIVE_PRIMARY",
            "Primary spatial validation on treatment-naïve samples. Section random effects are not identifiable (1 section per patient)."
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
            "LMM: Feature_Score ~ Moffitt50_contrast + CAF_Fraction + Myeloid_Fraction + Lymphoid_Fraction + (1 | patient_id)",
            "1000 label permutations, unrelated hallmark pathways (5), size-matched random gene sets (100)",
            "2026",
            "primary_hypotheses, secondary_localizations, spatial_autocorrelations",
            "ACTIVE_TREATMENT_SENSITIVITY",
            "Treatment sensitivity spatial validation. Section random effects are not identifiable."
        ],
        [
            "VAL_SPATIAL_CO_LOCAL_MONCADA",
            "Layer 3",
            "MONCADA_GSE111672",
            "GSE111672",
            "secondary_spatial",
            "HALLMARK_PROTEIN_SECRETION",
            "Moffitt50_contrast",
            "HALLMARK_SPERMATOGENESIS",
            "ELF1, MBD2, ZBTB7A, ZNF384, ZNF740",
            "spot",
            "decoupleR ssGSEA / VIPER (regulon-based)",
            "reference-based deconvolution (RCTD/MIA using matched scRNA-seq)",
            "LMM: Feature_Score ~ Moffitt50_contrast + CAF_Fraction + Myeloid_Fraction + Lymphoid_Fraction + (1 | patient_id) + (1 | patient_id:section_id)",
            "1000 spot-coordinate permutations, unrelated hallmark pathways (5), size-matched random gene sets (100)",
            "2026",
            "primary_hypotheses, secondary_localizations, spatial_autocorrelations",
            "ACTIVE_SECONDARY",
            "Secondary spatial validation. Patient A has 4 sections; Patient B has 2 sections. Nested random effects are used."
        ]
    ]
    write_tsv(path, header, rows)

def main():
    print("Generating Phase 9B3A spatial validation planning tables...")
    generate_feature_hierarchy()
    generate_dataset_qualification()
    generate_authoritative_cohort_set()
    generate_resource_estimate()
    generate_spatial_dataset_inventory()
    generate_spatial_parameter_inventory()
    print("All tables successfully generated.")

if __name__ == "__main__":
    main()
