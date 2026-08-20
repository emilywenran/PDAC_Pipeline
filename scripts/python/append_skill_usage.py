import os

file_path = "/Users/emily/thesis/PDAC/00_admin/SKILL_USAGE_LOG.tsv"

new_records = [
    "2026-07-02\tAntigravity\tPhase 9A\tdatabase-lookup\t/Users/emily/.agents/skills/database-lookup/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 01_metadata/external_validation_dataset_inventory.tsv, 05_results/tables/phase9a_external_dataset_shortlist.tsv, 05_results/tables/phase9a_external_analysis_resource_estimate.tsv\tLocated and read SKILL.md using view_file. Audited public databases (GEO, SRA, ENA, EGA) to identify human PDAC expression and microbiome cohorts.\tUsed to find and qualify independent bulk (TCGA-PAAD, GSE71729, GSE62452), single-cell (GSE111672, GSE154778, GSE202051), spatial (GSE202051, GSE274103, GSM3405527), and microbiome (PRJNA542615, EGAS00001004572) datasets.",
    "2026-07-02\tAntigravity\tPhase 9A\tcitation-management\t/Users/emily/.agents/skills/citation-management/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 09_docs/references/phase9_external_validation_sources.bib, 09_docs/references/phase9_external_validation_source_audit.tsv\tLocated and read SKILL.md using view_file. Verified DOIs, PMIDs, and accessions for the independent validation studies.\tUsed to compile and format the bibliography and audit records for 12 public datasets.",
    "2026-07-02\tAntigravity\tPhase 9A\texperimental-design\t/Users/emily/.agents/skills/experimental-design/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md, 09_docs/methods/PDAC_external_validation_protocol.md\tLocated and read SKILL.md using view_file. Reviewed principles of local control, blocking, replication, and batch/covariate confounding for external cohorts.\tUsed to design the validation hierarchy across four layers, define inclusion criteria, and specify negative controls to prevent technical bias.",
    "2026-07-02\tAntigravity\tPhase 9A\tstatistical-analysis\t/Users/emily/thesis/PDAC/.agents/skills/statistical-analysis/SKILL.md\tTrue\tFalse\t00_admin/SKILL_USAGE_LOG.tsv, 01_metadata/external_validation_parameter_inventory.tsv, 05_results/tables/phase9a_signature_external_coverage_feasibility.tsv\tLocated and read SKILL.md using view_file. Checked guidelines for correlation, linear regression, multiple testing corrections, and evidence categories.\tUsed to prespecify bulk, single-cell (pseudobulk), spatial, and microbiome validation statistics and evidence-category thresholds."
]

with open(file_path, "r") as f:
    content = f.read().splitlines()

# Append records
for r in new_records:
    if r not in content:
        content.append(r)

with open(file_path, "w") as f:
    f.write("\n".join(content) + "\n")

print("Appended Phase 9A records to SKILL_USAGE_LOG.tsv")
