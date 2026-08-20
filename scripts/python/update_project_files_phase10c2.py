import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT = Path("/Users/emily/thesis/PDAC")

# Update SKILL_USAGE_LOG.tsv
skill_log = ROOT / "00_admin/SKILL_USAGE_LOG.tsv"
try:
    df = pd.read_csv(skill_log, sep="\t")
except FileNotFoundError:
    df = pd.DataFrame(columns=["timestamp", "agent_id", "skill_name", "task_description", "success_status", "missing_dependencies"])

new_skills = [
    {
        "timestamp": datetime.now().isoformat() + "Z",
        "agent_id": "antigravity",
        "skill_name": "database-lookup",
        "task_description": "Phase 10C2 independent review of external database logs",
        "success_status": "TRUE",
        "missing_dependencies": ""
    },
    {
        "timestamp": datetime.now().isoformat() + "Z",
        "agent_id": "antigravity",
        "skill_name": "experimental-design",
        "task_description": "Phase 10C2 independent review of target prioritization design",
        "success_status": "TRUE",
        "missing_dependencies": ""
    },
    {
        "timestamp": datetime.now().isoformat() + "Z",
        "agent_id": "antigravity",
        "skill_name": "statistical-analysis",
        "task_description": "Phase 10C2 independent review of composition penalty implementation",
        "success_status": "TRUE",
        "missing_dependencies": ""
    }
]

df = pd.concat([df, pd.DataFrame(new_skills)], ignore_index=True)
df.to_csv(skill_log, sep="\t", index=False)

# Update PROJECT_STATUS.md
status_file = ROOT / "00_admin/PROJECT_STATUS.md"
content = status_file.read_text()

# Add to the top
new_status = f"""# Project Status

Phase 10C2 independent review of corrected Phase 10B-R cross-layer synthesis completed on {datetime.now().strftime("%Y-%m-%d")}. The audit confirmed that Phase 10B-R successfully addressed all Phase 10C failures. Target scores were programmatically derived from the locked Phase 10A inventories, descriptive manual overrides were removed, composition-sensitive evidence (e.g. CTCFL/BORIS, HALLMARK_SPERMATOGENESIS) was appropriately penalized, and external database provenance was tracked properly. The final decision is **`PASS`**, authorizing the project to proceed to manuscript drafting.

"""

content = content.replace("# Project Status\n", new_status)

# Also update the Next Approved Task section
if "## Next Approved Task\n\n- Phase 10C2 Independent Review: independently review the corrected Phase 10B-R target prioritization reanalysis." in content:
    content = content.replace(
        "## Next Approved Task\n\n- Phase 10C2 Independent Review: independently review the corrected Phase 10B-R target prioritization reanalysis.",
        "## Next Approved Task\n\n- Phase 11: Manuscript Drafting."
    )

new_update = f"""## Latest Update: {datetime.now().strftime("%Y-%m-%d")} Phase 10C2 Independent Review

- Completed Phase 10C2 independent review of the corrected Phase 10B-R cross-layer evidence synthesis.
- Verified that all descriptive hardcoded overrides were removed and candidate scores derived programmatically from locked Phase 10A rules.
- Confirmed that every eligible candidate in the Phase 10A inventory was scored.
- Verified that CELL_COMPOSITION_EXPLAINED was appropriately penalized.
- Verified that CTCFL/BORIS was NOT promoted and was correctly penalized.
- Verified that unavailable external database queries correctly reduced confidence without unauthorized hardcoding.
- Verified that HALLMARK_PROTEIN_SECRETION and BHLHE40 were evaluated objectively.
- Generated Phase 10C2 review report `04_analysis/10_target_prioritization/PHASE10C2_CORRECTED_TARGET_PRIORITIZATION_INDEPENDENT_REVIEW.md` and audit tables.
- Final decision: `PASS`. Ready for manuscript drafting.

"""

content = content.replace("## Latest Update: 2026-07-04 Phase 10B-R Corrected Reanalysis", new_update + "## Latest Update: 2026-07-04 Phase 10B-R Corrected Reanalysis")

status_file.write_text(content)

# File manifest update (since update_file_manifest.py already ran and appended 9 files, we check if the new tables are there)
manifest = ROOT / "01_metadata/file_manifest.tsv"
try:
    man_df = pd.read_csv(manifest, sep="\t")
except FileNotFoundError:
    pass

import hashlib
def file_md5(path):
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

files_to_add = [
    "04_analysis/10_target_prioritization/PHASE10C2_CORRECTED_TARGET_PRIORITIZATION_INDEPENDENT_REVIEW.md",
    "05_results/tables/phase10c2_scoring_rule_audit.tsv",
    "05_results/tables/phase10c2_candidate_inclusion_audit.tsv",
    "05_results/tables/phase10c2_external_database_audit.tsv",
    "05_results/tables/phase10c2_penalty_application_audit.tsv",
    "05_results/tables/phase10c2_rank_reproduction_audit.tsv",
    "05_results/tables/phase10c2_review_findings.tsv"
]

manifest_rows = []
for f in files_to_add:
    p = ROOT / f
    if p.exists():
        if f not in man_df['file_path'].values:
            manifest_rows.append({
                "phase": "10C2",
                "category": "analysis_report" if "04_analysis" in f else "review_audit",
                "file_path": f,
                "file_format": "markdown" if f.endswith(".md") else "tsv",
                "description": "Phase 10C2 Independent Review Output",
                "md5_checksum": file_md5(p)
            })

if manifest_rows:
    man_df = pd.concat([man_df, pd.DataFrame(manifest_rows)], ignore_index=True)
    man_df.to_csv(manifest, sep="\t", index=False)

# Decision Log Update
decision_file = ROOT / "09_docs/planning/DECISION_LOG.md"
dec_content = decision_file.read_text()
if "D-51" not in dec_content:
    new_decision = """

### D-51: Phase 10C2 Independent Review Final Decision

*   **Date**: 2026-07-04
*   **Phase**: Phase 10C2
*   **Description**: Phase 10B-R results passed the independent review. All errors from Phase 10C were successfully corrected.
*   **Rationale**: Target prioritization was fully reproducible, objective, and compliant with Phase 10A method locks.
*   **Impact**: Project can proceed to Phase 11: Manuscript Drafting.
"""
    decision_file.write_text(dec_content + new_decision)

print("Update completed.")
