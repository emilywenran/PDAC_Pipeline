# PDAC Microbiome–Transcriptome Plasticity — GitHub Release Package

This folder is a curated subset of the working project (`~/thesis/PDAC/`), selected for public release on GitHub. The original workspace was not modified; everything here is a copy. Selection criteria are documented in `docs/workflows/GITHUB_PUBLICATION_PACKAGE_PLAN.md`. **The submission package is intentionally excluded**, per project instructions.

## Project Summary

- Primary data sources: GEO `GSE172356` (host expression), NCBI BioProject `PRJNA719915` (intratumoral microbiome)
- Research goal: determine whether PDAC basal-like / classical / hybrid molecular subtypes are stable discrete states or positions on a continuous transcriptional gradient, and identify microbiome features robustly associated with host subtype state and candidate therapeutic vulnerabilities
- Methodological approach: an **LLM-assisted bioinformatics workflow** built on a "pre-registered method lock → agent execution → independent agent review → programmatic validator" cycle. Every design decision is logged in `docs/planning/DECISION_LOG.md`.

## Directory Structure

```
github_release/
├── README.md                    # this file
├── LICENSE                       # MIT
├── CITATION.cff                  # citation metadata
├── config.yaml                   # consolidated locked analysis parameters (primary runs)
├── environment.yml               # conda environment (Python + R)
├── requirements.txt              # Python dependencies only
├── run_pipeline.py                # top-level CLI: lists phases/scripts and dispatches to them
├── demo/
│   └── run_demo.py                # self-contained one-click demo (synthetic data, no downloads)
├── skills-lock.json               # K-Dense skill source lock (statistical-analysis)
├── docs/
│   ├── planning/                  # project charter, hypotheses, analysis plan, decision log, risk register
│   ├── methods/                   # per-phase method SOPs / protocols
│   ├── workflows/                 # workflow overview, agent/skill usage summary, GitHub publication plan, prompt-record audit
│   └── references/                # bibliography and source audit for external validation datasets
├── ai_usage/
│   ├── SKILL_USAGE_LOG.tsv        # per-phase K-Dense skill invocation log (including load failures / unavailable skills)
│   └── SKILLS_INVENTORY.tsv       # inventory of skills used in the project
├── project_status/
│   └── PROJECT_STATUS.md          # project status and completed-task log
├── metadata/                      # all accession / sample / parameter inventory tables from 01_metadata
├── reference_data/
│   └── PDAC_subtype_signatures/   # locked subtype signatures (Moffitt / PurIST / GSE172356 original)
├── analysis_reports/              # every method-lock / execution-result / independent-review report from 04_analysis, by phase
├── scripts/
│   ├── python/                    # analysis, summary, validator, and test scripts
│   ├── R/                         # R analysis and validation scripts
│   └── shell/                     # environment audit scripts
├── results/
│   ├── tables/                    # per-phase result and audit tables (299 files)
│   ├── figures/                   # result figures (PDF, 114 files)
│   └── reports/                   # initialization report
└── environment/                   # renv.lock, R environment info, software inventory (for R environment reproduction)
```

## Before Using This Package

1. **License and citation**: `LICENSE` (MIT) and `CITATION.cff` are included at the repository root. `CITATION.cff` will need its ORCID and publication-status fields updated once the manuscript is submitted.
2. **Hardcoded local paths inside scripts were left untouched**: roughly 40 scripts in `scripts/python/` and `scripts/R/` contain absolute paths from the original development machine (of the form `/Users/emily/thesis/PDAC/...`). To avoid silently changing script behavior, this release does **not** edit script contents — files were copied as-is. Anyone reproducing the analysis will need to adjust these paths manually (or run scripts from a checkout that mirrors the original project layout — see `run_pipeline.py` below). Documentation, log, and metadata files (non-code) already have local paths sanitized to `~/`.
3. **What's new in this release layer** (added on top of the existing project files, without modifying any existing script's internal logic):
   - `run_pipeline.py`: a lightweight CLI that lists all phases/scripts with descriptions and dispatches to the underlying phase script via subprocess. It does **not** rewrite any of the 92 existing analysis scripts into standalone CLI tools — those still expect to run from a full project checkout with `02_data/`, `03_processed/`, etc. present (which are intentionally excluded from this release; see below).
   - `config.yaml`: consolidates the primary/locked analysis parameters (random seeds, PERMANOVA permutations, CLR pseudocount, PurIST coefficients, FDR method, etc.) extracted directly from the `*_parameter_inventory.tsv` files in `metadata/`. Sensitivity-analysis parameter grids remain in the original TSVs; `config.yaml` is a human-readable summary of the primary configuration, not a replacement for them.
   - `environment.yml` / `requirements.txt`: consolidate the actual Python and R packages imported across `scripts/`, so the environment can be rebuilt with `conda env create -f environment.yml`.
   - `demo/run_demo.py`: a genuinely self-contained, one-command demo. It applies the real locked PurIST classifier coefficients (from `reference_data/PDAC_subtype_signatures/PurIST_signatures.tsv`) to a small **synthetic, seeded toy expression matrix** (explicitly labeled as synthetic — not real patient data, since the real GSE172356 matrix is not redistributed in this release). It exists to verify the environment is set up correctly and to illustrate the classifier logic end-to-end in a few seconds, not to reproduce the full validated pipeline.
4. **Not included**: raw/large sequencing data (`02_data/raw`, `02_data/external`), R package caches (`renv/library`, `renv/cache`), model binary objects (`05_results/models/`), local runtime caches (`.cache/`, `tmp/`, `__pycache__/`), `08_submission/` (submission package, excluded per current instructions), and private review artifacts (`.docx` annotations, archive folders, redundant zip files). See `docs/workflows/GITHUB_PUBLICATION_PACKAGE_PLAN.md`, Section 5, for the full rationale.
5. **Data/Code Availability and AI-usage disclosure**: the manuscript's Data Availability, Code Availability, and "Use of AI-Assisted Technologies" sections are filled in and reference this repository (`https://github.com/emilywenran/PDAC_Pipeline`). See `docs/workflows/GITHUB_PUBLICATION_PACKAGE_PLAN.md`, Section 7, for the original recommended wording this was based on.
6. **Prompt records**: no verbatim user prompts are stored anywhere in the project. `docs/workflows/PROMPT_AND_AGENT_EXECUTION_LOG.md` documents that finding and the closest available substitute (commit history + agent usage logs).

## Quick Start

```bash
# 1. Build the environment
conda env create -f environment.yml
conda activate pdac-github-release

# 2. Run the self-contained demo (no data download required, runs in seconds)
python demo/run_demo.py

# 3. List available analysis phases and their scripts
python run_pipeline.py --list

# 4. Show details for a specific phase (script path, language, description)
python run_pipeline.py --describe phase3b
```

## Suggested Next Steps

1. Review the contents of this folder to confirm nothing sensitive was inadvertently included.
2. Fill in the remaining `CITATION.cff` fields (ORCID, publication status) once available.
3. If you intend to re-run the full pipeline (not just the demo), restore the original project's data layout (`02_data/`, `03_processed/`) and update the hardcoded paths in `scripts/` as needed.
