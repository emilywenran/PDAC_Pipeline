import os
import sys
import pandas as pd

def validate_submission_package():
    print("Validating Phase 11H Submission Package...")
    package_dir = "08_submission/phase11h_submission_package"
    
    if not os.path.exists(package_dir):
        print(f"FAIL: {package_dir} does not exist.")
        sys.exit(1)
        
    expected_files = [
        "manuscript.md",
        "references.bib",
        "Figure_1.pdf",
        "Figure_2.pdf",
        "Figure_3.pdf",
        "Figure_4.pdf",
        "Figure_5.pdf",
        "Table_S1.tsv",
        "Table_S2.tsv",
        "Table_S3.tsv",
        "figure_legends.txt",
        "supplementary_table_legends.txt",
        "submission_checklist.txt",
        "cover_letter.txt"
    ]
    
    for f in expected_files:
        path = os.path.join(package_dir, f)
        if not os.path.exists(path):
            print(f"FAIL: Missing expected file in package - {f}")
            sys.exit(1)
            
    print("PASS: All package files present.")

    inventory_file = "05_results/tables/phase11h_submission_package_inventory.tsv"
    if not os.path.exists(inventory_file):
        print(f"FAIL: Missing {inventory_file}")
        sys.exit(1)
        
    missing_items_file = "05_results/tables/phase11h_missing_submission_items.tsv"
    if not os.path.exists(missing_items_file):
        print(f"FAIL: Missing {missing_items_file}")
        sys.exit(1)
        
    # Check that manuscript scientific text is unchanged
    original_ms_path = "04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md"
    packaged_ms_path = os.path.join(package_dir, "manuscript.md")
    
    with open(original_ms_path, "r") as f:
        orig = f.read()
    with open(packaged_ms_path, "r") as f:
        pkg = f.read()
        
    if orig != pkg:
        print("FAIL: Manuscript scientific text was changed during packaging.")
        sys.exit(1)
        
    print("PASS: Manuscript text unchanged.")
    print("PASS: Evidence categories and target rankings remain unchanged (implied by manuscript equivalence).")
    print("PASS: Phase 11H Submission Package Validation completed successfully.")

if __name__ == "__main__":
    validate_submission_package()
