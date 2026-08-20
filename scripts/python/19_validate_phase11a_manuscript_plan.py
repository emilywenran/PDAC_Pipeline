import sys
import re

def validate_plan():
    errors = []
    
    with open("04_analysis/11_manuscript/PHASE11A_MANUSCRIPT_CLAIM_MAP.md", "r") as f:
        claim_map = f.read().lower()
    
    with open("04_analysis/11_manuscript/PHASE11A_MANUSCRIPT_OUTLINE.md", "r") as f:
        outline = f.read().lower()
        
    with open("05_results/tables/phase11a_prohibited_claims.tsv", "r") as f:
        prohibited = f.read().lower()
        
    text = claim_map + "\n" + outline + "\n" + prohibited

    causal_words = [" causes ", " drives ", " induces ", " mediates ", " causality ", " causal "]
    for w in causal_words:
        if w in text:
            # We must be careful because the text might contain "Do not use causal language (e.g. "drives")"
            # Let's just check if it's used in a supportive context.
            pass

    # A better approach: check the exact constraints in the text.
    with open("04_analysis/11_manuscript/PHASE11A_MANUSCRIPT_CLAIM_MAP.md", "r") as f:
        claim_map_orig = f.read()
    with open("04_analysis/11_manuscript/PHASE11A_MANUSCRIPT_OUTLINE.md", "r") as f:
        outline_orig = f.read()

    combined_text = claim_map_orig + "\n" + outline_orig

    if not re.search(r'PARTIAL_SPATIAL_SUPPORT', combined_text):
        errors.append("Missing explicit classification of spatial validation as PARTIAL_SPATIAL_SUPPORT.")
        
    if "HALLMARK_SPERMATOGENESIS" not in combined_text or "MEred" not in combined_text or "MEpurple" not in combined_text:
        errors.append("Missing explicit reporting of null/unsupported evidence (e.g., HALLMARK_SPERMATOGENESIS, MEred, MEpurple).")

    if not re.search(r'Phase 10C2', combined_text, re.IGNORECASE):
        errors.append("Target prioritization does not explicitly reference Phase 10C2 ranking.")
        
    if re.search(r'(?i)promote\s+CTCFL', combined_text) or re.search(r'(?i)BORIS\s+is\s+prioritized', combined_text):
        errors.append("Unsupported promotion of CTCFL/BORIS detected.")

    if re.search(r'(?i)full spatial validation', combined_text):
        errors.append("Overstated spatial evidence detected ('full spatial validation').")

    if re.search(r'(?i)microbial localization', combined_text) and not re.search(r'(?i)(no|without|do not claim)\s+microbial localization', combined_text):
        errors.append("Microbial localization claimed incorrectly.")

    if re.search(r'(?i)causality', combined_text) and not re.search(r'(?i)(no|do not claim)\s+.*causality', combined_text):
        errors.append("Causal language used incorrectly.")

    if errors:
        for err in errors:
            print(f"VALIDATION FAILED: {err}")
        sys.exit(1)
    else:
        print("PHASE 11A VALIDATION PASSED")
        sys.exit(0)

if __name__ == "__main__":
    validate_plan()
