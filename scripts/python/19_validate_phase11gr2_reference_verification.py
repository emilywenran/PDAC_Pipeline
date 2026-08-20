import os
import sys
import csv

def main():
    audit_file = "04_analysis/11_manuscript/PHASE11GR2_REFERENCE_VERIFICATION_AUDIT.md"
    tsv_file = "05_results/tables/phase11gr2_reference_verification.tsv"
    manuscript = "04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md"

    if not os.path.exists(audit_file):
        print(f"FAIL: {audit_file} not found.")
        sys.exit(1)
        
    if not os.path.exists(tsv_file):
        print(f"FAIL: {tsv_file} not found.")
        sys.exit(1)

    with open(manuscript, "r") as f:
        text = f.read()
        
    # Find all citation numbers in manuscript (basic extraction of [1] to [14])
    # The manuscript contains [1], [2,3], [4,5], [6–9], [10], [11], [12,13], [14]
    expected_citations = ["[1]", "[2]", "[3]", "[4]", "[5]", "[6]", "[7]", "[8]", "[9]", "[10]", "[11]", "[12]", "[13]"]
    
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        
    mapped_citations = [row["citation_number"] for row in rows]
    
    for c in expected_citations:
        if c not in mapped_citations:
            print(f"FAIL: Citation {c} used in manuscript is missing from mapping table.")
            sys.exit(1)
            
    has_fail = False
    for row in rows:
        if not row["verification_status"]:
            print(f"FAIL: Row {row['citation_number']} missing verification_status.")
            sys.exit(1)
        if row["verification_status"] == "FAIL":
            has_fail = True
            
    if not has_fail:
        # We expect FAILS because we found some. If the agent silently ignored them, this fails.
        print("FAIL: Known FAIL rows (like [3], [6], [7], [13]) were silently ignored.")
        sys.exit(1)

    constraints = [
        "CELL_COMPOSITION_EXPLAINED",
        "PARTIAL_SPATIAL_SUPPORT",
        "exploratory analyses in an independent spatial cohort (Moncada)",
        "positive concordance in only 1 of 6 sections"
    ]
    for c in constraints:
        if c not in text:
            print(f"FAIL: Constraint '{c}' not found in manuscript text. Text may have been inappropriately edited.")
            sys.exit(1)

    print("PASS: Phase 11G-R2 reference verification validation successful.")

if __name__ == "__main__":
    main()
