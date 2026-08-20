import os
import pandas as pd
import sys

def run_validation():
    print("Validating Phase 10B cross-layer synthesis...")

    report_file = "04_analysis/10_target_prioritization/PHASE10B_TARGET_PRIORITIZATION_RESULTS.md"
    table_file = "05_results/tables/phase10b_candidate_target_scores.tsv"

    if not os.path.exists(report_file):
        print(f"FAIL: Report {report_file} not found.")
        sys.exit(1)

    if not os.path.exists(table_file):
        print(f"FAIL: Table {table_file} not found.")
        sys.exit(1)

    df = pd.read_csv(table_file, sep="\t")

    required_columns = ["feature_name", "final_synthesis_category", "druggability", "tumor_vs_normal", "priority_score"]
    for col in required_columns:
        if col not in df.columns:
            print(f"FAIL: Required column {col} missing from table.")
            sys.exit(1)

    if not df[df["feature_name"] == "CTCFL"].empty:
        ctcfl_row = df[df["feature_name"] == "CTCFL"].iloc[0]
        if "High" not in str(ctcfl_row["tumor_vs_normal"]):
            print("FAIL: CTCFL tumor_vs_normal evaluation incorrect.")
            sys.exit(1)

    with open(report_file, "r") as f:
        content = f.read()

    if "READY_FOR_MANUSCRIPT_DRAFTING" not in content:
        print("FAIL: Decision READY_FOR_MANUSCRIPT_DRAFTING not found in report.")
        sys.exit(1)

    print("Validation Passed: READY_FOR_MANUSCRIPT_DRAFTING")
    sys.exit(0)

if __name__ == "__main__":
    run_validation()
