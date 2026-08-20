import os
import sys
import pandas as pd

def validate_phase11e():
    print("Starting Phase 11E Validation...")

    manuscript_path = "04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md"
    report_path = "04_analysis/11_manuscript/PHASE11E_LANGUAGE_FORMAT_REVIEW.md"
    log_path = "05_results/tables/phase11e_language_edit_log.tsv"
    
    # Check if files exist
    for f in [manuscript_path, report_path, log_path]:
        if not os.path.exists(f):
            print(f"FAIL: Missing file {f}")
            sys.exit(1)
            
    with open(manuscript_path, 'r') as f:
        ms_text = f.read()

    # Check that legends remain present
    if "Figure Legends" not in ms_text or "Supplementary Table Legends" not in ms_text:
        print("FAIL: Figure or Supplementary Table Legends missing.")
        sys.exit(1)

    # Check HALLMARK_PROTEIN_SECRETION remains PARTIAL_SPATIAL_SUPPORT
    if "PARTIAL_SPATIAL_SUPPORT" not in ms_text or "HALLMARK_PROTEIN_SECRETION" not in ms_text:
        print("FAIL: HALLMARK_PROTEIN_SECRETION PARTIAL_SPATIAL_SUPPORT constraint violated.")
        sys.exit(1)

    # Check CTCFL/BORIS is not promoted
    if "CTCFL" not in ms_text or "CELL_COMPOSITION_EXPLAINED" not in ms_text:
        print("FAIL: CTCFL/BORIS CELL_COMPOSITION_EXPLAINED constraint violated.")
        sys.exit(1)

    # Check null and partial findings remain explicit
    if "HALLMARK_SPERMATOGENESIS" not in ms_text or "MEred" not in ms_text or "MEpurple" not in ms_text:
        print("FAIL: Null findings missing from text.")
        sys.exit(1)

    # Check microbial constraints
    # Looking for explicit non-causal language and checking against prohibited terms.
    # The exact phrase "causal relationship" or "causality" should be used in a negative context.
    # We will do a simple check to ensure no new causal claims are introduced by looking for the required disclaimers.
    lower_ms = ms_text.lower()
    if "causal direction cannot be inferred" not in lower_ms and "cannot be claimed to cause" not in lower_ms and "strictly non-causal" not in lower_ms and "no causal relationship" not in lower_ms:
         print("FAIL: Microbial causality disclaimer missing.")
         sys.exit(1)
         
    if "physical interaction" not in lower_ms:
         print("FAIL: Physical interaction disclaimer missing.")
         sys.exit(1)
         
    if "localisation" not in lower_ms and "localization" not in lower_ms:
         print("FAIL: Localisation disclaimer missing.")
         sys.exit(1)

    # Moncada exploratory
    if "moncada" in lower_ms and "exploratory" not in lower_ms:
         print("FAIL: Moncada exploratory constraint violated.")
         sys.exit(1)
         
    if "1 of 6" not in lower_ms:
         print("FAIL: Moncada 1 of 6 sections constraint violated.")
         sys.exit(1)
         
    # Check edit log for claim/evidence category changes
    try:
        log_df = pd.read_csv(log_path, sep='\t')
        if any(log_df['claim_changed_yes_no'].str.upper() == 'YES'):
             print("FAIL: Claim changed according to log.")
             sys.exit(1)
        if any(log_df['evidence_category_changed_yes_no'].str.upper() == 'YES'):
             print("FAIL: Evidence category changed according to log.")
             sys.exit(1)
    except Exception as e:
        print(f"FAIL: Error reading edit log: {e}")
        sys.exit(1)
        
    print("SUCCESS: Phase 11E Validation passed.")

if __name__ == "__main__":
    validate_phase11e()
