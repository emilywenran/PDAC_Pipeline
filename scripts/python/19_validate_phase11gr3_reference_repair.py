import os
import sys

def main():
    repair_md = "04_analysis/11_manuscript/PHASE11GR3_REFERENCE_REPAIR.md"
    repair_tsv = "05_results/tables/phase11gr3_reference_repair.tsv"
    map_file = "05_results/tables/phase11gr1_reference_mapping.tsv"
    manuscript = "04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md"
    bib_file = "09_docs/references/phase11g_references.bib"

    files_to_check = [repair_md, repair_tsv, map_file, manuscript, bib_file]
    for f in files_to_check:
        if not os.path.exists(f):
            print(f"FAIL: {f} not found.")
            sys.exit(1)

    with open(manuscript, "r") as f:
        text = f.read()

    with open(map_file, "r") as f:
        map_text = f.read()

    with open(bib_file, "r") as f:
        bib_text = f.read()

    # Verify PMIDs were updated
    if "32673551" in map_text:
        print("FAIL: PMID 32673551 for Bear2020 was not removed.")
        sys.exit(1)
    if "32471947" in map_text:
        print("FAIL: PMID 32471947 for Nejman2020 was not removed.")
        sys.exit(1)
    if "28916614" in map_text:
        print("FAIL: PMID 28916614 for Geller2017 was not removed.")
        sys.exit(1)

    if "32946773" not in map_text or "32467386" not in map_text or "28912244" not in map_text:
        print("FAIL: Updated PMIDs missing in map text.")
        sys.exit(1)

    # Verify Boehm2015 was removed
    if "Boehm2015" in bib_text:
        print("FAIL: Boehm2015 still present in .bib file.")
        sys.exit(1)
    if "Boehm2015" in map_text:
        print("FAIL: Boehm2015 still present in mapping table.")
        sys.exit(1)
    if "[12,13]" in text:
        print("FAIL: [12,13] still present in manuscript.")
        sys.exit(1)
    if "[14]" in text:
        print("FAIL: [14] still present in manuscript (should have been renumbered).")
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

    print("PASS: Phase 11G-R3 reference repair validation successful.")

if __name__ == "__main__":
    main()
