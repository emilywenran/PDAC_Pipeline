#!/usr/bin/env python3
import os
import pandas as pd

def validate_phase11ia():
    print("Starting Phase 11I-A Validation...")
    
    pkg_dir = "08_submission/phase11h_submission_package/"
    if not os.path.exists(pkg_dir):
        raise ValueError(f"Missing package directory: {pkg_dir}")
        
    required_files = [
        "manuscript.md", "references.bib", "Figure_1.pdf", "Figure_2.pdf", 
        "Figure_3.pdf", "Figure_4.pdf", "Figure_5.pdf", "Table_S1.tsv", 
        "Table_S2.tsv", "Table_S3.tsv", "figure_legends.txt", 
        "supplementary_table_legends.txt", "submission_checklist.txt", "cover_letter.txt"
    ]
    for f in required_files:
        if not os.path.exists(os.path.join(pkg_dir, f)):
            raise ValueError(f"Missing file in package: {f}")
            
    with open(os.path.join(pkg_dir, "manuscript.md"), "r") as f:
        ms_text = f.read()
    with open(os.path.join(pkg_dir, "references.bib"), "r") as f:
        bib_text = f.read()
        
    callouts = [
        "Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5",
        "Table S1", "Table S2", "Table S3"
    ]
    for c in callouts:
        if c not in ms_text:
            raise ValueError(f"Missing callout for {c} in manuscript.")
            
    bib_count = bib_text.count("@article{")
    if bib_count != 13:
        raise ValueError(f"Expected 13 bibliography entries, found {bib_count}.")
        
    constraints = [
        "we do not claim microbial causality, microbial localisation, or physical host-microbe interaction",
        "remaining PARTIAL_SPATIAL_SUPPORT",
        "CELL_COMPOSITION_EXPLAINED",
        "exploratory analyses in an independent spatial cohort (Moncada)",
        "Ochrobactrum",
        "HALLMARK_SPERMATOGENESIS"
    ]
    for c in constraints:
        if c not in ms_text:
            raise ValueError(f"Missing claim-control phrase: {c}")
            
    qa_path = "05_results/tables/phase11ia_final_qa_checklist.tsv"
    if not os.path.exists(qa_path):
        raise ValueError(f"Missing QA checklist: {qa_path}")
    qa_df = pd.read_csv(qa_path, sep='\t')
    if not all(col in qa_df.columns for col in ['qa_item', 'category', 'status', 'evidence_path', 'issue', 'required_action']):
        raise ValueError("QA checklist missing required columns.")
        
    gap_path = "05_results/tables/phase11ia_journal_specific_gap_table.tsv"
    if not os.path.exists(gap_path):
        raise ValueError(f"Missing gap table: {gap_path}")
    gap_df = pd.read_csv(gap_path, sep='\t')
    if not all(col in gap_df.columns for col in ['item', 'current_status', 'why_journal_specific', 'information_needed_from_user', 'priority', 'required_action_after_journal_selected']):
        raise ValueError("Gap table missing required columns.")
    if not all(gap_df['current_status'] == 'TO_BE_CONFIRMED'):
        raise ValueError("All items in gap table must be marked TO_BE_CONFIRMED.")
        
    with open("04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md", "r") as f:
        orig_ms = f.read()
    if orig_ms != ms_text:
        raise ValueError("Manuscript text was altered.")
        
    print("Phase 11I-A Validation Passed!")

if __name__ == "__main__":
    validate_phase11ia()
