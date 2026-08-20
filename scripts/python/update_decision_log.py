#!/usr/bin/env python3
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
LOG_PATH = ROOT / "09_docs/planning/DECISION_LOG.md"

def main():
    content = LOG_PATH.read_text()
    
    # Let's locate the table row for D-27 and insert D-28, D-29, D-30, and D-31 after it.
    target_line = "| 2026-07-02 | **D-27** | Identify, evaluate, and lock the external-validation framework (Phase 9A) across four layers under a prospective method lock | `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_external_validation_protocol.md`, `01_metadata/external_validation_*`, `05_results/tables/phase9a_*` |"
    
    if target_line not in content:
        print("Error: Target D-27 table row not found")
        return
        
    table_additions = """| 2026-07-03 | **D-28** | Execute Locked Phase 9B1 Independent Bulk-Transcriptome Validation | `06_scripts/R/14_phase9b1_bulk_validation.R`, `04_analysis/09_external_validation/PHASE9B1_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1_*` |
| 2026-07-03 | **D-29** | Reject Phase 9B1 Bulk Validation Outputs due to Critical and Major Implementation Errors (Phase 9B1C) | `04_analysis/09_external_validation/PHASE9B1C_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c_*` |
| 2026-07-03 | **D-30** | Correct and Rerun Phase 9B1 Bulk Validation After Phase 9B1C Audit (Phase 9B1R) | `06_scripts/R/14_phase9b1r_corrected_bulk_validation.R`, `04_analysis/09_external_validation/PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1r_*` |
| 2026-07-03 | **D-31** | Complete Phase 9B1C2 independent review of the corrected Phase 9B1R bulk external validation results and report PASS_WITH_MINOR_CORRECTIONS decision | `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c2_*` |"""

    new_table_block = target_line + "\n" + table_additions
    content = content.replace(target_line, new_table_block)
    
    # Now let's append the detailed log for D-31 to the end of the file.
    detailed_log_addition = """
---

### D-31: Complete Phase 9B1C2 Independent Review of Corrected Phase 9B1R Bulk Validation
*   **Date:** 2026-07-03
*   **Decision:** Accept the corrected Phase 9B1R bulk-transcriptome external validation results under the final review decision `PASS_WITH_MINOR_CORRECTIONS` and approve proceeding to Phase 9B2 single-cell validation.
*   **Alternatives Considered:** Issue `FAIL_REQUIRES_REANALYSIS` if any critical errors remained; issue `PASS` without minor corrections and leave TFs hardcoded as TO_VERIFY.
*   **Scientific and Operational Justification:** Accepting the corrected results under `PASS_WITH_MINOR_CORRECTIONS` is scientifically justified because the calculations are now mathematically correct and all six findings from the previous audit have been successfully resolved. The minor correction involves the final evidence table's TF category reporting, which the reviewer independently resolved by reclassifying the 34 TFs using the successfully executed VIPER activity statistics. Proceeding to Phase 9B2 is approved since no critical or major implementation errors remain.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c2_*`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`.
"""
    content = content.rstrip() + "\n" + detailed_log_addition
    
    LOG_PATH.write_text(content)
    print("Successfully updated DECISION_LOG.md")

if __name__ == "__main__":
    main()
