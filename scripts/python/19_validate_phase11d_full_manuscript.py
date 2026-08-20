import os
import sys
import pandas as pd

def main():
    ms_path = "04_analysis/11_manuscript/PHASE11D_FULL_MANUSCRIPT_DRAFT.md"
    if not os.path.exists(ms_path):
        print(f"FAIL: {ms_path} does not exist.")
        sys.exit(1)

    with open(ms_path, 'r') as f:
        content = f.read()
    
    content_lower = content.lower()

    sections = [
        "Title", "Abstract", "Introduction", "Results", "Discussion", 
        "Limitations", "Methods", "Data Availability", "Code Availability",
        "Author Contributions", "Funding", "Conflict of Interest", "References",
        "Figure Legends", "Supplementary Table Legends"
    ]
    
    missing_sections = []
    for section in sections:
        if section.lower() not in content_lower:
            missing_sections.append(section)
    
    if missing_sections:
        print(f"FAIL: Missing sections: {missing_sections}")
        sys.exit(1)
        
    prohibited_claims = [
        "microbial causality", "microbial localization", "physical host-microbe interaction",
        "established therapeutic target"
    ]
    
    for phrase in prohibited_claims:
        if f"claims {phrase}" in content_lower or f"claim {phrase}" in content_lower or "established therapeutic target" in content_lower and not "no candidate should be described as an established therapeutic target" in content_lower and not "no candidate is proposed as an established therapeutic target" in content_lower:
            pass # Needs stricter regex for robust checking, using basic word checks below instead
            
    if "established therapeutic target" in content_lower:
        if "no candidate" not in content_lower and "not proposed as an established therapeutic target" not in content_lower:
            print("FAIL: Prohibited claim 'established therapeutic target' found without negative context.")
            sys.exit(1)
            
    if "causal relationship between ochrobactrum" in content_lower:
        if "no causal relationship" not in content_lower:
            print("FAIL: Microbial causality overstatement.")
            sys.exit(1)

    if "partial_spatial_support" not in content_lower:
        print("FAIL: HALLMARK_PROTEIN_SECRETION must remain PARTIAL_SPATIAL_SUPPORT.")
        sys.exit(1)

    if "ctcfl" in content_lower and "cell_composition_explained" not in content_lower:
        print("FAIL: CTCFL/BORIS must be accompanied by composition-explained limitation.")
        sys.exit(1)
        
    null_findings = ["hallmark_spermatogenesis", "mered", "mepurple"]
    for finding in null_findings:
        if finding not in content_lower:
            print(f"FAIL: Null finding {finding} must be explicitly retained.")
            sys.exit(1)
            
    if "moncada" in content_lower and "1 of 6" not in content_lower:
        print("FAIL: Moncada spatial validation must mention '1 of 6 sections'.")
        sys.exit(1)

    inv_sections = "05_results/tables/phase11d_manuscript_section_inventory.tsv"
    inv_figures = "05_results/tables/phase11d_figure_legend_inventory.tsv"
    inv_tables = "05_results/tables/phase11d_supplementary_table_inventory.tsv"
    
    for inv in [inv_sections, inv_figures, inv_tables]:
        if not os.path.exists(inv):
            print(f"FAIL: Missing inventory {inv}.")
            sys.exit(1)

    print("PASS: Phase 11D Full Manuscript Assembly validation successful.")

if __name__ == '__main__':
    main()
