# GitHub Publication Package Plan — PDAC Microbiome–Transcriptome Project

**Purpose**: plan what to include when publishing this project on GitHub as an example of an LLM-assisted bioinformatics workflow — what to keep, what to leave out, the recommended repository layout, and what to disclose in the manuscript. This is a planning document only: it does not change any scientific result, the manuscript, or admin records, and it does not perform any git operations.

**Sources reviewed**: `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv` (116 rows), `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md` (decisions D-01 through D-56+), `04_analysis/`, `05_results/`, `06_scripts/`, `08_submission/`, `.gitignore`, the git log, and the pre-existing `09_docs/workflows/PDAC_PROJECT_WORKFLOW_AGENT_SKILL_SUMMARY.md` (Sections 1–3 below condense that document; see it for full detail).

---

## 1. Project Workflow Overview

The project is built around two primary data sources — GEO GSE172356 (host expression) and BioProject PRJNA719915 (tumor microbiome) — and proceeds through the following stages:

| Stage | Phase IDs | Content |
| :--- | :--- | :--- |
| 0. Data provenance audit | Phase 1A–1B | Accession audit; 1:1 crosswalk between RNA-seq and microbiome samples for 62 patients |
| 1. Primary bulk analysis (subtype + continuous axis) | Phase 2–5 | Expression QC/filtering, reproduction of the 94-gene CSY subtype labels, stability evaluation (INCONCLUSIVE), continuous basal–classical transcriptional axis scoring |
| 2. Primary microbiome analysis | Phase 6–8 | PRJNA719915 abundance audit, compositional (CLR) preprocessing, microbiome–host state associations (OLS/HC3, PERMANOVA), host-mechanism analysis (Hallmark/PROGENy/DoRothEA/WGCNA) |
| 3. External bulk validation | Phase 9A, 9B1→9B1C→9B1R→9B1C2 | Independent replication in TCGA-PAAD / GSE71729 / GSE62452; the first attempt was rejected for TF-proxy errors (9B1C), fixed in 9B1R, and passed independent review in 9B1C2 |
| 4. Single-cell validation | Phase 9A.1–9A.3, 9B2→9B2C→9B2R→9B2C2 | Pseudobulk validation on PENG_CRA001160; the first attempt was rejected for unexecuted negative controls (9B2C), fixed in 9B2R, and passed in 9B2C2 |
| 5. Spatial validation | Phase 9B3A–9B3A.2, 9B3B→9B3C→9B3R→9B3C2 | Spatial validation on Hwang GeoMx (naive/treated) + Moncada ST; the first attempt was rejected for hardcoded negative controls (9B3C), fixed in 9B3R, and passed in 9B3C2 |
| 6. Cross-layer evidence synthesis | Phase 10A→10B→10C→10B-R→10C2 | Multi-layer evidence scoring; the first attempt was rejected for undocumented external database queries and post hoc target selection (10C), fixed in 10B-R, and passed in 10C2 |
| 7. Target prioritisation | Alongside Phase 10 | OpenTargets / GTEx / ChEMBL queries; composite scoring on druggability, selectivity, and safety |
| 8. Manuscript drafting & claim-control audit | Phase 11A–11F | Claim map, full manuscript assembly, language editing, independent claim-control audit (15 constraints, 100% compliant) |
| 9. Reference / callout audit | Phase 11G(-R1/R2/R3) | Citation and figure/table callout repair and verification |
| 10. Submission package assembly | Phase 11H | Manuscript, figures, tables, cover letter, and checklist packaged together |
| 11. Final QA / journal gap audit | Phase 11I-A (followed by an 11J language-humanization pass) | Completeness and journal-gap audit; final status `READY_FOR_TARGET_JOURNAL_SELECTION` |

The governance pattern running through the whole pipeline is: **method lock (pre-registration) → execution → independent review (repeatedly producing `FAIL_REQUIRES_REANALYSIS` → fix → PASS cycles) → programmatic validator script → decision logged**, with 56+ formal decisions recorded in `DECISION_LOG.md`. This "pre-register, independently audit, and programmatically verify" loop is the project's most notable methodological feature and the strongest candidate for public showcasing.

---

## 2. Agent Usage Summary

Based on `SKILL_USAGE_LOG.tsv` (the `agent` field records only two role names, `Antigravity` and `Codex`) and the correspondence already documented in the pre-existing summary:

| Role name in the logs | Corresponding model (per existing project records) | Main responsibility | Typical phases |
| :--- | :--- | :--- | :--- |
| `Antigravity` | Gemini 3.1 Pro High | High-level decision-making and process control: writing method locks and analysis protocols, independent review (9B1C/9B2C/9B3C/10C, etc.), claim-control audits, citation/callout audits | Phase 1A, 3A, 4A, 5A, 7A/7C, 8A/8C, all of 9A, 9B1C/9B1C2, 9B2C/9B2C2, all of 9B3A/9B3C/9B3C2, 10A/10B/10C2, 11A/11C/11D/11F/11G, 11H, 11I-A |
| `Codex` | Gemini 3.5 Flash / Flash High | Concrete code implementation and computation: matrix cleaning, mixed-effects model fitting, association regressions, plotting | Phase 2A/2B, 4B, 5B, 6A/6C, 7A.5/7B, 8B, 9B1/9B1R, 9B2/9B2R, 9B3R, 10B-R, 11B |

Additionally, per the existing summary document (not present as structured `SKILL_USAGE_LOG.tsv` rows — only mentioned in free text inside individual review reports such as `PHASE9B3R0_DUAL_MODEL_REANALYSIS_AUDIT.md`, which is outside the scope of this review):

| Role | Use | Note |
| :--- | :--- | :--- |
| Claude | Statistical review of the Phase 9B3R0 spatial reanalysis | Not structurally logged in `SKILL_USAGE_LOG.tsv`; only referenced in the body text of individual audit reports |
| ChatGPT | Implementation-level code scan for the Phase 9B3R0 reanalysis; also reportedly used for interactive workflow/prompt design across the project | Same caveat as above, and its "interactive control / prompt design" role **is not documented anywhere verifiable in the project records** |

**Not documented in the project records** (should be disclosed honestly in the GitHub README rather than glossed over):
- The executing agent for Phase 1B and Phase 3B was not recorded.
- The specific scope of Claude / ChatGPT involvement, prompt content, and interaction records were never systematically archived — only scattered references in individual audit reports.
- The exact prompts/instructions used to hand off tasks between agents were never archived (only execution results and the `SKILL_USAGE_LOG.tsv` summaries exist).

---

## 3. K-Dense-AI Skills Usage Summary

Based on `SKILL_USAGE_LOG.tsv` (116 rows, mostly with `skill_loaded=True`):

| Skill | Phases used (summary) | Main purpose | Core? |
| :--- | :--- | :--- | :--- |
| `database-lookup` | 1A, 3A, 6A, 9A, 9A.1, 9A.2, 9B1, 9B2 (including restart/execution), 9B3A, 9B3B, 10C2 | Querying NCBI E-utilities / ENA / GEO / SRA metadata; verifying sample counts, platforms, and supplementary file structure | Yes |
| `exploratory-data-analysis` | 2A, 2B, 6A, 6B, 6C, 7A.5, 9B1, all of 9B2, 9B3B, 9B3C, 9B3R | EDA of expression/microbiome matrices: missingness, zero-inflation, outlier samples, distribution visualization | Yes |
| `citation-management` | 3A, 8A, 9A, 9A.2, 9B3A, 11G | BibTeX compilation, DOI/PMID verification | Yes |
| `experimental-design` | 3A, 4A, 6B, 7A, 7C, 8A, 8C, all of 9A, 9B1C, 9B2C/9B2R/9B2C2, all of 9B3A/9B3B/9B3C/9B3R, 10C2, 11F | Covariate design, nested random-effects structures, pseudoreplication avoidance, negative-control design | Yes |
| `statistical-analysis` | 4A, 4B, 5A, 5B, 6B, 6C, 7A/7A.5/7B/7C, 8A/8B/8C, 9A, all of 9B1, all of 9B2, all of 9B3A/9B3B/9B3C/9B3R/9B3C2, 10A, 10C2 | Fitting OLS(HC3), nested LMMs, Moran's I, BH-FDR correction, bootstrap/JT trend tests | Yes |
| `opentargets-database` | 10B | Tractability and disease-association queries for candidate targets | No (supporting) |
| `gtex-database` | 10B | Normal-tissue specificity queries for candidate targets | No (supporting) |
| `chembl-database` | 10B | Compound lookup for candidate targets (produced no substantive results in this project) | No (supporting) |
| `scientific-writing` | 11B | Manuscript writing; claim-trace and wording enforcement | Yes |

**Attempts recorded as unavailable/missing** (disclosed here for methodological transparency, not omitted):
- Phase 9B3C2: both `experimental-design` and `exploratory-data-analysis` are logged as `MISSING` (skill file absent from `.agents/skills/`).
- Phase 10A: `experimental-design`, `citation-management`, and `database-lookup` are logged as `UNAVAILABLE`.
- Phase 11A: a request for the overarching `K-Dense` controller skill is logged as `unavailable`.
- Phase 11C/11D/11E/11H/11I-A: the skill field in the log is `manuscript-review` / `manuscript-writing` / `manuscript-assembly` / `quality-assurance` with `exact_skill_path=N/A` — i.e., **no loadable K-Dense skill file is documented for these entries**, only the task description and output files.

---

## 4. Recommended for GitHub Publication

| Category | Path(s) | Rationale |
| :--- | :--- | :--- |
| README | New GitHub-facing `README.md` (distinct from the current lab-internal short version) | Project summary, workflow diagram, agent/skill usage summary, entry point for reproduction instructions |
| Workflow documentation | `09_docs/workflows/`, `09_docs/methods/*_protocol.md`, `04_analysis/*/PHASE*_METHOD_LOCK.md` and each phase's `*_RESULTS.md` / `*_INDEPENDENT_REVIEW.md` | Shows the full method-lock → execute → independent-review loop; this is the project's core methodological selling point |
| Decision log | `09_docs/planning/DECISION_LOG.md` | 56+ logged design decisions with rationale; demonstrates auditability |
| Analysis scripts | `06_scripts/python/`, `06_scripts/R/`, `06_scripts/shell/`, `06_scripts/workflow/` | All consistently named execution/summary scripts |
| Validators | `06_scripts/python/*_validate_*.py`, `06_scripts/python/test_*.py` | Programmatic validator for every phase — key evidence of reproducibility |
| K-Dense skill usage record | `00_admin/SKILL_USAGE_LOG.tsv` | Per-phase agent/skill invocation log, including failure/unavailable entries (should not be trimmed to look better) |
| Dataset/accession metadata | `01_metadata/*.tsv` (accession_inventory, sample_manifest, `*_crosswalk`, `*_parameter_inventory`, external_validation_dataset_inventory, etc.) | Public accessions and de-identified patient-ID crosswalks only — no raw sequence data |
| Selected result tables/figures | `05_results/tables/` (per-phase QC/association/audit/evidence-classification tables), `05_results/figures/*.pdf` | Recommend including in full (~49 MB total, manageable size) since these are exactly the "validation trail" evidence |
| Selected submission-package content | *(excluded from this release per current instructions — see the accompanying release folder's README)* | — |
| Reproducibility instructions | `renv.lock`, `07_envs/phase8_r_environment.yml`, `07_envs/phase8_R_sessionInfo.txt`, a Python dependency list (needs a new `requirements.txt`) | Environment lock files enabling others to rebuild the analysis environment |
| Project status / audit trail | `00_admin/PROJECT_STATUS.md` (after a privacy check), `01_metadata/file_manifest.tsv` (paths need to be made relative) | Demonstrates project governance and a complete output inventory |

---

## 5. Not Recommended for GitHub Publication

| Category | Path(s) | Rationale |
| :--- | :--- | :--- |
| Raw / large data | `02_data/raw/`, `02_data/external/` (4.0 GB; raw/processed files for GSE/TCGA/single-cell/spatial cohorts) | Too large, and mostly re-downloadable processed matrices from public databases; should be represented by accession numbers instead |
| R package caches / environment libraries | `renv/library/`, `renv/cache/`, `renv/staging/` (2.0 GB), `07_envs/R_bootstrap_lib/` (3.4 MB), `07_envs/R_user_cache/` (85 MB) | Large binary caches with no scientific content; already partially covered by `.gitignore` |
| Model / intermediate binary objects | `05_results/models/phase4b/`, `phase8b/`, `phase9b2/` | Already listed in `.gitignore`; large and not human-readable — should be recomputed on reproduction instead |
| Local caches / temporary files | `.cache/`, `tmp/`, `.matplotlib/`, `.pytest_cache/`, `__pycache__/`, `.DS_Store`, `08_logs/*.log` | Pure runtime artifacts with no publication value |
| Credentials / secrets | None found in this review (no `.env`, `*credential*`, `*.key`, etc.) | Still recommend a final scan with a tool such as `git secrets` / `trufflehog` before publishing |
| Private / local path references | Absolute-path columns in `01_metadata/file_manifest.tsv` (`/Users/emily/thesis/PDAC/...`); `file:///Users/emily/...` links in `09_docs/workflows/PDAC_PROJECT_WORKFLOW_AGENT_SKILL_SUMMARY.md` | Exposes the local username and directory layout; must be converted to relative paths before publishing |
| Redundant / private review files | `08_submission/phase11h_submission_package/manuscript_review.docx`, `manuscript_humanized_review*.docx`, `08_submission/phase11h_submission_package/archive/`, `comments_PDAC_AI_final.docx`, root-level `MS_PDAC.zip`, `PDAC_revised.zip` | Private reviewer annotations and duplicate packaging artifacts — working files, not publication material; the `.docx` comments may contain reviewer/collaborator private opinions |
| Unnecessary intermediate output | `10_manuscript/` (empty directory), `~$AC_Workflow_Summary.docx` (Word lock file) | No content, or a temporary lock file |
| Skill implementation code | `.agents/skills/` (already excluded by `.gitignore`) | K-Dense skill implementations are third-party/private assets — only the **usage log** (`SKILL_USAGE_LOG.tsv`) should be published, not the skill source itself |

---

## 6. Recommended GitHub Repository Layout

```
pdac-microbiome-transcriptome/
├── README.md                        # overview + workflow diagram + agent/skill summary + reproduction guide
├── LICENSE
├── CITATION.cff
├── requirements.txt                 # Python dependencies
├── renv.lock                        # R dependency lock
├── docs/
│   ├── workflow/
│   │   ├── PDAC_PROJECT_WORKFLOW_AGENT_SKILL_SUMMARY.md
│   │   └── GITHUB_PUBLICATION_PACKAGE_PLAN.md
│   ├── methods/                     # *_protocol.md (per-phase SOPs)
│   ├── decision_log/
│   │   └── DECISION_LOG.md
│   └── ai_usage/
│       └── SKILL_USAGE_LOG.tsv
├── metadata/                        # curated 01_metadata (paths made relative)
│   ├── accession_inventory.tsv
│   ├── sample_manifest.tsv
│   ├── *_crosswalk.tsv
│   └── *_parameter_inventory.tsv
├── analysis_reports/                # method-lock / result / independent-review reports from 04_analysis, by phase
│   ├── 01_data_audit/ ... 11_manuscript/
├── scripts/
│   ├── python/
│   ├── R/
│   ├── shell/
│   └── workflow/
├── validators/                      # optionally split out of scripts/ (*_validate_*.py / test_*.py), or keep in place with a README note
├── results/
│   ├── tables/                      # selected or full result tables
│   └── figures/                     # PDF figures
└── submission_package/
    ├── manuscript.md
    ├── references.bib
    ├── figure_legends.txt
    ├── supplementary_table_legends.txt
    ├── Figure_1-5.pdf
    └── Table_S1-S3.tsv
```

---

## 7. Recommended Additions to the Manuscript

1. **An `LLM-assisted bioinformatics workflow` Methods subsection**
   Briefly describe the "pre-registered method lock → agent execution → independent agent review → programmatic validator" governance cycle used throughout the analysis; describe the division of labor between agent roles (decision/review vs. code implementation) and the scope of K-Dense-AI skills used; explicitly note that all statistical parameters were locked before execution, and that independent review repeatedly triggered `FAIL_REQUIRES_REANALYSIS` and drove rework (e.g., Phase 9B1C, 9B2C, 9B3C, 10C) — cite this as evidence of methodological robustness.

2. **A `Use of AI-assisted technologies` disclosure statement**
   The current manuscript (`04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md`) does **not yet contain** this disclosure. Most target journals now require it. Recommend adding a statement naming the AI systems/agent roles involved (explicitly noting that some agent attributions are **not fully documented in the project records** — do not overstate completeness), describing which stages involved AI assistance (data processing, statistical modeling, writing/editing, independent review), and affirming that human authors take responsibility for the final scientific conclusions and manuscript content.

3. **A `Code Availability` / GitHub repository statement**
   The manuscript's `## Code Availability` section is currently a placeholder (`[Placeholder for Code Availability]`), as is `## Data Availability` (`[Placeholder for Data Availability]`). Once the target journal and GitHub repository are finalized, fill these in with the repository URL, license type, a version/DOI (e.g., via Zenodo archival), and the public accession numbers for the primary data and all Phase 9 external validation cohorts (GSE172356, PRJNA719915, and the external cohort accessions).

---

## Report Summary

- **Output file path**: `09_docs/workflows/GITHUB_PUBLICATION_PACKAGE_PLAN.md`
- **Main content recommended for publication**: workflow/method-lock/independent-review documentation (`04_analysis/`, `09_docs/`), all analysis and validator scripts (`06_scripts/`), `SKILL_USAGE_LOG.tsv`, `DECISION_LOG.md`, curated accession/metadata tables (`01_metadata/`), `05_results/tables` and `figures`, the finalized manuscript/figures/tables in `08_submission/phase11h_submission_package`, and environment lock files such as `renv.lock`.
- **Missing or undocumented agent/skill records**: the executing agent for Phase 1B and Phase 3B was never recorded; Claude and ChatGPT involvement appears only in the body text of individual audit reports, not in structured logs; the `experimental-design`/`exploratory-data-analysis` skills for Phase 9B3C2, the `experimental-design`/`citation-management`/`database-lookup` skills for Phase 10A, and the overarching `K-Dense` controller skill for Phase 11A are all logged as missing/unavailable; the skill file paths for Phase 11C/11D/11E/11H/11I-A are recorded as `N/A`.
- **Ready for packaging?**: **Largely yes, but three cleanup items should be completed first**: (1) remove or rewrite the local absolute paths / `file:///Users/emily/...` links in `01_metadata/file_manifest.tsv` and `09_docs/workflows/PDAC_PROJECT_WORKFLOW_AGENT_SKILL_SUMMARY.md`; (2) exclude `MS_PDAC.zip`, `PDAC_revised.zip`, `comments_PDAC_AI_final.docx`, `08_submission/phase11h_submission_package/archive/`, and the various `manuscript_*review*.docx` private review files from the publication scope; (3) add `LICENSE`, `CITATION.cff`, `requirements.txt`, and fill in the manuscript's Data/Code Availability and AI-assisted-technologies disclosure placeholders. Once these three items are complete, the content can be organized per Section 6 and published.
