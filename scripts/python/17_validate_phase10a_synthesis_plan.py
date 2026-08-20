import os
import sys
import pandas as pd

def validate_phase10a():
    errors = []

    # 1. Check if files exist
    expected_files = [
        "04_analysis/10_target_prioritization/PHASE10A_CROSS_LAYER_SYNTHESIS_METHOD_LOCK.md",
        "09_docs/methods/PDAC_cross_layer_synthesis_protocol.md",
        "01_metadata/phase10a_evidence_synthesis_parameter_inventory.tsv",
        "05_results/tables/phase10a_cross_layer_evidence_inventory.tsv",
        "05_results/tables/phase10a_target_prioritization_framework.tsv"
    ]
    for f in expected_files:
        if not os.path.exists(f):
            errors.append(f"Missing file: {f}")

    # 2. Check TSV structures
    if os.path.exists("05_results/tables/phase10a_cross_layer_evidence_inventory.tsv"):
        df_inv = pd.read_csv("05_results/tables/phase10a_cross_layer_evidence_inventory.tsv", sep="\t")
        expected_cols_inv = ["feature_name", "discovery_evidence", "bulk_evidence", "sc_evidence", "spatial_evidence", "final_synthesis_category"]
        for col in expected_cols_inv:
            if col not in df_inv.columns:
                errors.append(f"Inventory missing column: {col}")
        
        if "HALLMARK_PROTEIN_SECRETION" in df_inv["feature_name"].values:
            secret_row = df_inv[df_inv["feature_name"] == "HALLMARK_PROTEIN_SECRETION"].iloc[0]
            if secret_row["final_synthesis_category"] != "MULTI_LAYER_SUPPORTED":
                errors.append("HALLMARK_PROTEIN_SECRETION must be MULTI_LAYER_SUPPORTED")
        else:
            errors.append("HALLMARK_PROTEIN_SECRETION missing from inventory")
            
        # check WGCNA modules
        wgcna_modules = ["MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"]
        for mod in wgcna_modules:
            if mod in df_inv["feature_name"].values:
                mod_row = df_inv[df_inv["feature_name"] == mod].iloc[0]
                if mod_row["final_synthesis_category"] != "INSUFFICIENT_DATA":
                    errors.append(f"{mod} must be INSUFFICIENT_DATA")

    if os.path.exists("05_results/tables/phase10a_target_prioritization_framework.tsv"):
        df_fw = pd.read_csv("05_results/tables/phase10a_target_prioritization_framework.tsv", sep="\t")
        expected_fw_cols = ["criteria", "description", "data_source", "weight", "threshold"]
        for col in expected_fw_cols:
            if col not in df_fw.columns:
                errors.append(f"Framework missing column: {col}")

    # 3. Check parameters
    if os.path.exists("01_metadata/phase10a_evidence_synthesis_parameter_inventory.tsv"):
        df_param = pd.read_csv("01_metadata/phase10a_evidence_synthesis_parameter_inventory.tsv", sep="\t")
        expected_param_cols = ["parameter_name", "parameter_value", "description", "source_phase"]
        for col in expected_param_cols:
            if col not in df_param.columns:
                errors.append(f"Parameters missing column: {col}")

    if errors:
        print("Validation Failed:")
        for e in errors:
            print(f" - {e}")
        sys.exit(1)
    else:
        print("Validation Passed: READY_FOR_PHASE10B_CROSS_LAYER_SYNTHESIS")
        sys.exit(0)

if __name__ == "__main__":
    validate_phase10a()
