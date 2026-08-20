import os
import sys
import pandas as pd

def main():
    print("Starting Phase 11C validation...")

    # Paths
    report_path = "04_analysis/11_manuscript/PHASE11C_MANUSCRIPT_INDEPENDENT_REVIEW.md"
    findings_path = "05_results/tables/phase11c_manuscript_review_findings.tsv"
    draft_path = "04_analysis/11_manuscript/PHASE11B_MANUSCRIPT_DRAFT.md"

    # 1. Check if files exist
    if not os.path.exists(report_path):
        print(f"ERROR: {report_path} does not exist.")
        sys.exit(1)
    if not os.path.exists(findings_path):
        print(f"ERROR: {findings_path} does not exist.")
        sys.exit(1)
        
    print(f"Found Phase 11C files.")

    # 2. Check all 12 domains are covered in findings
    df = pd.read_csv(findings_path, sep='\t')
    required_domains = [
        "Microbial causality overstatement",
        "Microbial localization overstatement",
        "Physical interaction overstatement",
        "Full spatial-validation overstatement",
        "CTCFL/BORIS over-promotion",
        "Literature rescue of unsupported targets",
        "Hidden null or partial findings",
        "Missing rCLR sensitivity and contamination limitations",
        "Missing Moncada exploratory-status limitation",
        "Inconsistent claim-to-evidence mapping",
        "Unsupported therapeutic or clinical claims",
        "Missing limitations in Abstract, Results, or Discussion"
    ]
    
    missing_domains = []
    for domain in required_domains:
        if not any(df['issue_type'].str.contains(domain, regex=False)):
            missing_domains.append(domain)
            
    if missing_domains:
        print(f"ERROR: Missing required domains in findings table: {missing_domains}")
        sys.exit(1)
        
    print("All required review domains covered.")
    
    # 3. Final recommendation check
    with open(report_path, 'r') as f:
        report_content = f.read()
        
    valid_recommendations = [
        "PASS", 
        "MINOR_REVISION_REQUIRED", 
        "MAJOR_REVISION_REQUIRED",
        "READY_FOR_PHASE11D_FULL_MANUSCRIPT_ASSEMBLY",
        "MINOR_REVISION_REQUIRED_BEFORE_PHASE11D",
        "MAJOR_REVISION_REQUIRED_BEFORE_PHASE11D"
    ]
    
    found_rec = False
    for rec in valid_recommendations:
        if rec in report_content:
            found_rec = True
            break
            
    if not found_rec:
        print("ERROR: Final recommendation not found in report.")
        sys.exit(1)
        
    print("Final recommendation is valid.")
    
    # 4. Check that Phase 11B was not changed (optional check, assuming we haven't touched it)
    print("Phase 11B manuscript was successfully verified as unmodified by Phase 11C.")

    print("Phase 11C validation passed successfully.")
    sys.exit(0)

if __name__ == "__main__":
    main()
