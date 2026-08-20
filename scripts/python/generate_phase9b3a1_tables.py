#!/usr/bin/env python3
"""Generate Phase 9B3A.1 Spatial Analysis Unit and Models TSV."""

import os
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
TABLE_PATH = ROOT / "05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv"

def main():
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
            "segment",
            "256",
            "patient",
            "segment",
            "segment",
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction + (1 | patient_id)",
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + (1 | patient_id)",
            "(1 | patient_id)",
            "CAF_fraction, myeloid_fraction, lymphoid_fraction",
            "treatment-naïve samples only",
            "10",
            "variance_inflation_factor (VIF) > 10",
            "If model fails to converge, fit reduced model. If still non-convergent, perform patient-level regression and aggregate coefficients via random-effects meta-analysis.",
            "ROIs are physically segmented into PanCK+ (epithelial) and PanCK- (stroma) areas. Section random effect is not identifiable because there is exactly 1 section per patient."
        ],
        [
            "HWANG_GSE202051_TREATED",
            "GSE202051",
            "NanoString GeoMx DSP",
            "25",
            "25",
            "segment",
            "352",
            "patient",
            "segment",
            "segment",
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction + (1 | patient_id)",
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + (1 | patient_id)",
            "(1 | patient_id)",
            "CAF_fraction, myeloid_fraction, lymphoid_fraction",
            "neoadjuvant-treated samples only; analyzed separately to assess treatment-sensitivity effects",
            "10",
            "variance_inflation_factor (VIF) > 10",
            "If model fails to converge, fit reduced model or fall back to patient-level regression and meta-analysis.",
            "Subjected to neoadjuvant therapy. Useful to assess treatment remodeling effects. Section random effect is not identifiable."
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
            "protein_secretion_score ~ Moffitt50_contrast + CAF_fraction + myeloid_fraction + lymphoid_fraction + (1 | patient_id) + (1 | patient_id:section_id)",
            "protein_secretion_score ~ Moffitt50_contrast + (1 | patient_id) + (1 | patient_id:section_id)",
            "(1 | patient_id) + (1 | patient_id:section_id)",
            "CAF_fraction, myeloid_fraction, lymphoid_fraction",
            "treatment-naïve samples only",
            "2",
            "variance_inflation_factor (VIF) > 10",
            "Due to low patient count (n=2), fit reduced model or fall back to patient-level regression followed by fixed-effects or random-effects meta-analysis of the coefficients.",
            "Patient A has 4 sections; Patient B has 2 sections. Nested random effects account for multiple sections per patient."
        ]
    ]
    
    TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with TABLE_PATH.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(val) for val in row) + "\n")
            
    print(f"Wrote {len(rows)} rows to {TABLE_PATH}")

if __name__ == "__main__":
    main()
