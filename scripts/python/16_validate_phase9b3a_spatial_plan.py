#!/usr/bin/env python3
"""Phase 9B3A.2 Spatial Validation Plan Validator.

This validator checks the consistency of the spatial validation planning tables,
ensuring all constraints on accessions, patient/section counts, features,
analysis units, model classes, covariates, ROI pairing, Moncada exploratory reclassification,
and execution roles are strictly enforced.
"""

import sys
import csv
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")

# Paths to files
FEATURE_HIERARCHY_PATH = ROOT / "05_results/tables/phase9b3a_spatial_feature_hierarchy.tsv"
DATASET_QUALIFICATION_PATH = ROOT / "05_results/tables/phase9b3a_spatial_dataset_qualification.tsv"
COHORT_SET_PATH = ROOT / "05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv"
RESOURCE_ESTIMATE_PATH = ROOT / "05_results/tables/phase9b3a_spatial_resource_estimate.tsv"
SPATIAL_INVENTORY_PATH = ROOT / "01_metadata/phase9b3_spatial_dataset_inventory.tsv"
SPATIAL_PARAMETER_PATH = ROOT / "01_metadata/phase9b3_spatial_parameter_inventory.tsv"
SPATIAL_MODELS_PATH = ROOT / "05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv"
MODEL_HIERARCHY_PATH = ROOT / "05_results/tables/phase9b3a2_spatial_model_hierarchy.tsv"

# Path to Markdown planning documents
METHOD_LOCK_PATH = ROOT / "04_analysis/09_external_validation/PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK.md"
PROTOCOL_PATH = ROOT / "09_docs/methods/PDAC_spatial_validation_protocol.md"

def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        print(f"Error: Required file not found: {path}", file=sys.stderr)
        sys.exit(1)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def check_term_in_file(path: Path, term: str) -> bool:
    if not path.exists():
        return False
    content = path.read_text(encoding="utf-8").lower()
    return term.lower() in content

def main() -> int:
    print("Executing Phase 9B3A.2 Spatial Validation Plan Validator...")
    errors: list[str] = []

    # 1. Load tables
    hierarchy = read_tsv(FEATURE_HIERARCHY_PATH)
    qualification = read_tsv(DATASET_QUALIFICATION_PATH)
    cohort_set = read_tsv(COHORT_SET_PATH)
    resource = read_tsv(RESOURCE_ESTIMATE_PATH)
    inventory = read_tsv(SPATIAL_INVENTORY_PATH)
    parameters = read_tsv(SPATIAL_PARAMETER_PATH)
    models_table = read_tsv(SPATIAL_MODELS_PATH)
    model_hierarchy = read_tsv(MODEL_HIERARCHY_PATH)

    # 2. Rule: Stated covariates and model formula disagree
    for row in models_table:
        covs = [c.strip().lower() for c in row.get("composition_covariates", "").split(",") if c.strip()]
        model = row.get("primary_model", "").lower()
        if "moncada" not in row.get("canonical_dataset_id", "").lower(): # Skip Moncada as it's exploratory section-specific
            for cov in covs:
                if cov not in model:
                    errors.append(f"Covariate and formula disagreement in models table for {row.get('canonical_dataset_id')}: covariate '{cov}' is not in model '{model}'")

    for row in parameters:
        model = row.get("statistical_model", "").lower()
        if "moncada" not in row.get("dataset_id", "").lower():
            if "lymphoid_fraction" not in model and "lymphoid" not in model:
                errors.append(f"Lymphoid covariate missing from parameter statistical model for {row.get('analysis_id')}: model='{row.get('statistical_model')}'")

    # 3. Rule: Lymphoid adjustment is required in prose but omitted from the full model
    for path in [METHOD_LOCK_PATH, PROTOCOL_PATH]:
        if not check_term_in_file(path, "lymphoid_fraction") and not check_term_in_file(path, "lymphoid"):
            errors.append(f"Lymphoid adjustment term missing from prose in {path.name}")
            
    # Check that model hierarchy Model A and B contain lymphoid
    for row in model_hierarchy:
        model_id = row.get("model_id", "")
        formula = row.get("model_formula", "")
        if ("MODEL_A" in model_id or "MODEL_B" in model_id) and "lymphoid" not in formula.lower():
            errors.append(f"Lymphoid adjustment missing from {model_id} formula: '{formula}'")

    # 4. Rule: Patient count and section count are stored in one ambiguous field
    for path in [SPATIAL_INVENTORY_PATH, DATASET_QUALIFICATION_PATH, RESOURCE_ESTIMATE_PATH, SPATIAL_MODELS_PATH, MODEL_HIERARCHY_PATH]:
        if path.exists():
            with path.open("r", encoding="utf-8") as f:
                header = f.readline().strip().split("\t")
                for col in header:
                    if "/" in col and ("patient" in col.lower() or "section" in col.lower()):
                        errors.append(f"Ambiguous patient/section column header '{col}' in {path.name}")

    # 5. Rule: ROI or spot is treated as the biological replicate / spots or sections are treated as patients
    for row in models_table:
        replicate = row.get("biological_replicate", "").lower()
        if "patient" not in replicate:
            errors.append(f"Biological replicate unit error for {row.get('canonical_dataset_id')}: '{replicate}', must be 'patient'")
        
        # Stated unit counts should not define spot/ROI as replicate
        unit = row.get("inferential_unit", "").lower()
        if unit == "patient":
            errors.append(f"Inferential unit is patient in models table, expected spot/segment nested within patient")
            
    for row in model_hierarchy:
        replicate = row.get("replicate_unit", "").lower()
        if "patient" not in replicate:
            errors.append(f"Biological replicate unit error in model hierarchy for {row.get('model_id')}: '{replicate}', must be 'patient'")

    # 6. Rule: Hwang compartment models omit ROI pairing
    # Model A (Compartment comparison LMM) must contain (1 | patient_id:ROI_id) or equivalent paired block
    hwang_model_a = [r for r in model_hierarchy if "HWANG" in r.get("canonical_dataset_id", "") and "MODEL_A" in r.get("model_id", "")]
    for row in hwang_model_a:
        re_str = row.get("random_effect_structure", "")
        if "roi_id" not in re_str.lower() or "patient_id" not in re_str.lower():
            errors.append(f"Hwang compartment comparison model '{row.get('model_id')}' omits ROI pairing or patient nesting: '{re_str}'")

    # 7. Rule: Tumor/stroma segments are treated as independent
    # In both Hwang naive and treated, ensure that segments are nested using (1 | patient_id)
    for row in model_hierarchy:
        if "HWANG" in row.get("canonical_dataset_id", ""):
            re_str = row.get("random_effect_structure", "")
            if "patient_id" not in re_str.lower():
                errors.append(f"Hwang model '{row.get('model_id')}' treats segments/ROIs as independent, missing patient random effect: '{re_str}'")

    # 8. Rule: Moncada is described as formal replication
    # Ensure Moncada's role is EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY (not SECONDARY or formal replication)
    moncada_cohorts = [r for r in cohort_set if "MONCADA" in r.get("canonical_dataset_id", "")]
    for row in moncada_cohorts:
        role = row.get("spatial_analysis_role", "")
        if "replication" in role.lower() or "secondary" in role.lower() or role == "SECONDARY":
            errors.append(f"Moncada is incorrectly labeled as a formal replication: '{role}'")
        if role != "EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY":
            errors.append(f"Moncada has incorrect exploratory role label: '{role}'")

    # 9. Rule: Hwang and Moncada matrices are pooled
    # Ensure deconvolution/pooling notes do not suggest direct merging
    for row in parameters:
        note = row.get("notes", "").lower()
        if "pool" in note and "matrix" in note:
            errors.append(f"Direct matrix pooling suggestion detected in parameter notes: '{row.get('notes')}'")
    for row in models_table:
        note = row.get("notes", "").lower()
        if "pool count matrices" in note or "direct merge" in note:
            errors.append(f"Direct matrix pooling suggestion detected in models table notes: '{row.get('notes')}'")

    # 10. Rule: Naïve and treated Hwang samples are pooled in the primary model
    hwang_cohorts = [r for r in cohort_set if "HWANG" in r.get("canonical_dataset_id", "")]
    if len(hwang_cohorts) < 2:
        errors.append("Hwang cohorts not split in authoritative cohort set")
    else:
        naive = [r for r in hwang_cohorts if "NAIVE" in r.get("canonical_dataset_id", "")]
        treated = [r for r in hwang_cohorts if "TREATED" in r.get("canonical_dataset_id", "")]
        if not naive or not treated:
            errors.append("Hwang cohort names do not split NAIVE and TREATED groups")
        
        # Verify execution roles are distinct
        naive_role = naive[0].get("spatial_analysis_role", "")
        treated_role = treated[0].get("spatial_analysis_role", "")
        if naive_role == treated_role:
            errors.append(f"Hwang Naive and Treated cohorts share the same role '{naive_role}' in cohort set")

    # 11. Rule: Stated model is incorrectly described as 'mixed-effects OLS'
    for path in [SPATIAL_PARAMETER_PATH, SPATIAL_MODELS_PATH, MODEL_HIERARCHY_PATH, METHOD_LOCK_PATH, PROTOCOL_PATH]:
        if check_term_in_file(path, "mixed-effects OLS") or check_term_in_file(path, "mixed effects OLS"):
            errors.append(f"Incorrect model terminology 'mixed-effects OLS' found in {path.name}")

    # 12. Rule: Active roles are distinct and include exploratory reclassification
    roles = set()
    for row in cohort_set:
        role = row.get("spatial_analysis_role", "")
        if row.get("current_execution_authorized", "").upper() == "TRUE":
            roles.add(role)
    expected_roles = {"PRIMARY", "EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY", "TREATMENT_SENSITIVITY"}
    for r in expected_roles:
        if r not in roles:
            errors.append(f"Active execution set is missing distinct role '{r}'")

    # 13. Reconcile Provenance Mappings (Moncada, Hwang)
    expected_publications = {
        "GSE111672": "Moncada et al. (2020)",
        "GSE202051": "Hwang et al. (2022)"
    }
    for row in inventory:
        acc = row.get("accession", "")
        pub = row.get("publication", "")
        if acc in expected_publications:
            expected = expected_publications[acc]
            if expected.lower() not in pub.lower():
                errors.append(f"Publication mismatch: accession {acc} mapped to '{pub}', expected '{expected}' in spatial dataset inventory")

    # 14. Feature Hierarchy Frozen Verification
    secretion_row = [r for r in hierarchy if r.get("feature_name") == "HALLMARK_PROTEIN_SECRETION"]
    if not secretion_row:
        errors.append("HALLMARK_PROTEIN_SECRETION is missing from feature hierarchy")
    else:
        row = secretion_row[0]
        if row.get("primary_or_secondary", "").upper() != "PRIMARY":
            errors.append("HALLMARK_PROTEIN_SECRETION must be classified as 'PRIMARY'")
        if row.get("spatial_analysis_role", "") != "PRIMARY_TARGET":
            errors.append("HALLMARK_PROTEIN_SECRETION spatial role must be 'PRIMARY_TARGET'")

    # 15. Summarize validation results
    if errors:
        print("Validation FAILED. Errors found:")
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1
    else:
        print("Validation PASSED. All Phase 9B3A.2 prospective rules are successfully locked.")
        return 0

if __name__ == "__main__":
    sys.exit(main())
