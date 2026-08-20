import os
import pandas as pd
import sys

def main():
    citation_file = "05_results/tables/phase11g_citation_audit.tsv"
    callout_file = "05_results/tables/phase11g_figure_table_callout_audit.tsv"
    package_file = "05_results/tables/phase11g_submission_package_inventory.tsv"
    audit_md = "04_analysis/11_manuscript/PHASE11G_REFERENCE_SUBMISSION_AUDIT.md"

    files_to_check = [citation_file, callout_file, package_file, audit_md]
    for f in files_to_check:
        if not os.path.exists(f):
            print(f"FAIL: {f} not found.")
            sys.exit(1)

    # Validate Citation Audit
    cit_df = pd.read_csv(citation_file, sep='\t')
    cit_cols = [
        "manuscript_section",
        "citation_or_placeholder",
        "nearby_claim",
        "citation_status",
        "required_action",
        "source_file_if_available"
    ]
    for col in cit_cols:
        if col not in cit_df.columns:
            print(f"FAIL: {citation_file} missing column {col}")
            sys.exit(1)

    # Validate Callout Audit
    call_df = pd.read_csv(callout_file, sep='\t')
    call_cols = [
        "item_type",
        "item_id",
        "callout_present_yes_no",
        "legend_present_yes_no",
        "section_location",
        "consistency_status",
        "required_action"
    ]
    for col in call_cols:
        if col not in call_df.columns:
            print(f"FAIL: {callout_file} missing column {col}")
            sys.exit(1)
            
    if not (call_df["callout_present_yes_no"] == "No").all():
        print("FAIL: Expected all callouts to be 'No' based on the manuscript text.")
        sys.exit(1)

    # Validate Package Audit
    pkg_df = pd.read_csv(package_file, sep='\t')
    pkg_cols = [
        "package_item",
        "current_file_path",
        "status",
        "required_action",
        "notes"
    ]
    for col in pkg_cols:
        if col not in pkg_df.columns:
            print(f"FAIL: {package_file} missing column {col}")
            sys.exit(1)

    print("PASS: Phase 11G reference and submission audit validation successful.")

if __name__ == "__main__":
    main()
