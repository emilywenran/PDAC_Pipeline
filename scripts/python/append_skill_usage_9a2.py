#!/usr/bin/env python3
import os

file_path = "/Users/emily/thesis/PDAC/00_admin/SKILL_USAGE_LOG.tsv"

new_records = [
    "2026-07-03\tAntigravity\tPhase 9A.2\tdatabase-lookup\t/Users/emily/.agents/skills/database-lookup/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 01_metadata/external_validation_dataset_inventory.tsv, 05_results/tables/phase9a2_single_cell_dataset_provenance_audit.tsv, 05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv\tLocated and read SKILL.md. Queried GEO and CNCB-NGDC GSA to resolve accessions CRA001160 (Peng 2019) and GSE111672 (Moncada 2020), and re-audited GSE154778.\tUsed to correct dataset inventory metadata and resolve cohort mapping inconsistencies.",
    "2026-07-03\tAntigravity\tPhase 9A.2\tcitation-management\t/Users/emily/.agents/skills/citation-management/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 09_docs/references/phase9_external_validation_sources.bib, 09_docs/references/phase9_external_validation_source_audit.tsv\tLocated and read SKILL.md. Verified DOIs, PMIDs, and accessions for Peng 2019, Moncada 2020, Lin 2020, and Hwang 2022.\tUsed to update references and bibliography files for corrected single-cell cohorts.",
    "2026-07-03\tAntigravity\tPhase 9A.2\texperimental-design\t/Users/emily/.agents/skills/experimental-design/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 04_analysis/09_external_validation/PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION.md, 04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md, 09_docs/methods/PDAC_external_validation_protocol.md\tLocated and read SKILL.md. Evaluated experimental design, cell-type composition, treatment status, and compartment availability.\tUsed to design the expanded four-cohort single-cell validation set for Phase 9B2."
]

with open(file_path, "r") as f:
    content = f.read().splitlines()

# Append records
for r in new_records:
    if r not in content:
        content.append(r)

with open(file_path, "w") as f:
    f.write("\n".join(content) + "\n")

print("Appended Phase 9A.2 records to SKILL_USAGE_LOG.tsv")
