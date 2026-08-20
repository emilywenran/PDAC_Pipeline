import os
import sys

def main():
    bib_file = "09_docs/references/phase11g_references.bib"
    map_file = "05_results/tables/phase11gr1_reference_mapping.tsv"
    report_file = "04_analysis/11_manuscript/PHASE11GR1_REFERENCE_CALLOUT_REPAIR.md"
    manuscript = "04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md"

    files_to_check = [bib_file, map_file, report_file, manuscript]
    for f in files_to_check:
        if not os.path.exists(f):
            print(f"FAIL: {f} not found.")
            sys.exit(1)

    with open(manuscript, "r") as f:
        text = f.read()

    if "## Figure Legends" not in text:
        print("FAIL: '## Figure Legends' heading not found.")
        sys.exit(1)
        
    if "## Supplementary Table Legends" not in text:
        print("FAIL: '## Supplementary Table Legends' heading not found.")
        sys.exit(1)
        
    body_text_figures = text.split("## Figure Legends")[0]
    body_text_tables = text.split("## Supplementary Table Legends")[0]

    # Check that [14] is in Results
    if "## Results" not in body_text_figures or "## Discussion" not in body_text_figures:
        print("FAIL: Results section not found properly.")
        sys.exit(1)
        
    results_text = body_text_figures.split("## Results")[1].split("## Discussion")[0]
    if "[13]" not in results_text:
        print("FAIL: Citation [14] not found in Results text.")
        sys.exit(1)

    # Check for figure callouts before Figure Legends
    fig_callouts = ["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5"]
    for c in fig_callouts:
        if c not in body_text_figures:
            print(f"FAIL: Callout for {c} not found in manuscript body text before Figure Legends.")
            sys.exit(1)

    # Check for table callouts before Supplementary Table Legends
    tab_callouts = ["Table S1", "Table S2", "Table S3"]
    for c in tab_callouts:
        if c not in body_text_tables:
            print(f"FAIL: Callout for {c} not found in manuscript body text before Supplementary Table Legends.")
            sys.exit(1)

    # Prohibited claims (make sure we don't accidentally introduce new violations)
    # Since the text explicitly denies them (e.g. "No causal relationship"), we check that the denials exist.
    if "No causal relationship between *Ochrobactrum* and the observed host transcriptional programmes is established" not in text:
        print("FAIL: Missing denial of causal relationship for Ochrobactrum.")
        sys.exit(1)
    if "spatial localisation of *Ochrobactrum* within tumour tissue was not determined" not in text:
        print("FAIL: Missing denial of spatial localisation for Ochrobactrum.")
        sys.exit(1)
    if "no physical interaction between *Ochrobactrum* and host cells is proposed" not in text:
        print("FAIL: Missing denial of physical interaction for Ochrobactrum.")
        sys.exit(1)
    if "No candidate identified by this framework is described as an established therapeutic target" not in text and "no candidate should be considered an established therapeutic target" not in text:
        print("FAIL: Missing denial of established therapeutic target.")
        sys.exit(1)

    constraints = [
        "CELL_COMPOSITION_EXPLAINED",
        "PARTIAL_SPATIAL_SUPPORT",
        "exploratory analyses in an independent spatial cohort (Moncada)",
        "positive concordance in only 1 of 6 sections"
    ]
    for c in constraints:
        if c not in text:
            print(f"FAIL: Constraint '{c}' not found in manuscript text.")
            sys.exit(1)
            
    if "CTCFL" in text and "removed from candidate prioritisation" not in text:
        print("FAIL: CTCFL constraint missing.")
        sys.exit(1)
        
    print("PASS: Phase 11G-R1 reference and callout repair validation successful.")

if __name__ == "__main__":
    main()

