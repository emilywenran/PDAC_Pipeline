import os
import pandas as pd
import sys

def main():
    audit_file = "05_results/tables/phase11f_claim_audit.tsv"
    scan_file = "05_results/tables/phase11f_prohibited_claim_scan.tsv"

    if not os.path.exists(audit_file):
        print(f"FAIL: {audit_file} not found.")
        sys.exit(1)

    if not os.path.exists(scan_file):
        print(f"FAIL: {scan_file} not found.")
        sys.exit(1)

    audit_df = pd.read_csv(audit_file, sep='\t')
    expected_cols = [
        "manuscript_section",
        "claim_text_or_location",
        "audit_item",
        "status",
        "required_action",
        "evidence_category_changed_yes_no",
        "target_ranking_changed_yes_no"
    ]
    for col in expected_cols:
        if col not in audit_df.columns:
            print(f"FAIL: {audit_file} missing column {col}")
            sys.exit(1)

    scan_df = pd.read_csv(scan_file, sep='\t')
    expected_scan_cols = [
        "term_or_pattern",
        "context",
        "allowed_or_problematic",
        "reason",
        "required_action"
    ]
    for col in expected_scan_cols:
        if col not in scan_df.columns:
            print(f"FAIL: {scan_file} missing column {col}")
            sys.exit(1)

    if not (audit_df["status"] == "PASS").all():
        print("FAIL: Not all items in audit_file have status PASS.")
        sys.exit(1)
        
    print("PASS: Phase 11F final claim audit validation successful.")

if __name__ == "__main__":
    main()
