#!/usr/bin/env python3
"""Top-level CLI for the PDAC Microbiome-Transcriptome Plasticity project.

This is a thin dispatcher, not a rewrite of the underlying analysis scripts.
It does two things:
  1. Lists every analysis phase in this repository with its scripts, report,
     and a short description (--list, --describe).
  2. Optionally invokes a phase's script via subprocess (--run), after checking
     whether the data layout the script expects is actually present.

Why a dispatcher instead of rewriting the 92 underlying scripts as standalone
CLI tools: those scripts were written and validated against the full project
data layout (02_data/, 03_processed/, etc.), which is intentionally excluded
from this GitHub release (see README.md and docs/workflows/
GITHUB_PUBLICATION_PACKAGE_PLAN.md, Section 5). Rewriting each script's
internals risked silently changing behavior that produced the published
results. This wrapper adds discoverability and a safe run path without
touching any existing script's logic.

For a fully self-contained, no-data-download demonstration, use demo/run_demo.py
instead of this CLI.

Usage:
    python run_pipeline.py --list
    python run_pipeline.py --describe phase9b1r
    python run_pipeline.py --run phase3b
    python run_pipeline.py --run phase3b --force   # run even if expected inputs are missing
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# Explicit phase -> script/report mapping, built from the repository's own
# script-naming conventions, the decision log, and the skill-usage log
# (see docs/workflows/PDAC_PROJECT_WORKFLOW_AGENT_SKILL_SUMMARY.md).
# `requires_full_data_layout` scripts expect 02_data/ and/or 03_processed/
# from the original project checkout, which are not part of this release.
PHASES: dict[str, dict] = {
    "phase1a": {
        "description": "Public data-source and accession provenance audit (GSE172356, PRJNA719915).",
        "python": [],
        "r": [],
        "report": "analysis_reports/01_data_audit/PHASE1A_SOURCE_AUDIT.md",
        "note": "Documentation/audit only in this repository — no dedicated numbered script.",
        "requires_full_data_layout": False,
    },
    "phase2a": {
        "description": "Audit the official GSE172356 processed expression matrix.",
        "python": ["scripts/python/03_phase2a_expression_audit.py", "scripts/python/03_validate_phase2a_expression.py"],
        "r": [],
        "report": "analysis_reports/03_expression_qc/PHASE2A_EXPRESSION_AUDIT.md",
        "requires_full_data_layout": True,
    },
    "phase2b": {
        "description": "Missingness audit, filtering, and analysis-ready expression matrix preparation.",
        "python": ["scripts/python/04_phase2b_prepare_expression.py", "scripts/python/04_validate_phase2b_expression.py"],
        "r": [],
        "report": "analysis_reports/03_expression_qc/PHASE2B_ANALYSIS_READY_EXPRESSION.md",
        "requires_full_data_layout": True,
    },
    "phase3b": {
        "description": "Reproduce the 94-gene CSY / Moffitt / PurIST subtype assignments.",
        "python": ["scripts/python/05_phase3b_reproduce_subtypes.py", "scripts/python/05_validate_phase3b_subtypes.py"],
        "r": [],
        "report": "analysis_reports/05_subtype_reproduction/PHASE3B_SUBTYPE_REPRODUCTION.md",
        "requires_full_data_layout": True,
    },
    "phase4b": {
        "description": "Molecular subtype clustering stability evaluation (K=2..6, 1000 resamples).",
        "python": ["scripts/python/06_summarize_phase4b_stability.py", "scripts/python/06_validate_phase4b_stability.py"],
        "r": ["scripts/R/06_phase4b_subtype_stability.R"],
        "report": "analysis_reports/06_subtype_stability/PHASE4B_SUBTYPE_STABILITY_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase5b": {
        "description": "Continuous basal-classical transcriptional axis scoring and trend tests.",
        "python": ["scripts/python/07_summarize_phase5b_axis.py", "scripts/python/07_validate_phase5b_axis.py"],
        "r": ["scripts/R/07_phase5b_continuous_axis.R"],
        "report": "analysis_reports/07_continuous_subtype_axis/PHASE5B_CONTINUOUS_AXIS_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase6a": {
        "description": "Audit the processed PRJNA719915 tumor microbiome abundance matrix.",
        "python": ["scripts/python/08_phase6a_microbiome_audit.py", "scripts/python/08_validate_phase6a_microbiome.py"],
        "r": [],
        "report": "analysis_reports/04_microbiome_qc/PHASE6A_MICROBIOME_DATA_AUDIT.md",
        "requires_full_data_layout": True,
    },
    "phase6c": {
        "description": "Microbiome preprocessing: prevalence filtering, CLR transform, Aitchison distance.",
        "python": [
            "scripts/python/09_phase6c_prepare_microbiome.py",
            "scripts/python/09_summarize_phase6c_microbiome.py",
            "scripts/python/09_validate_phase6c_microbiome.py",
        ],
        "r": ["scripts/R/09_phase6c_prepare_microbiome.R"],
        "report": "analysis_reports/04_microbiome_qc/PHASE6C_ANALYSIS_READY_MICROBIOME.md",
        "requires_full_data_layout": True,
    },
    "phase7a5": {
        "description": "ESTIMATE-derived tumor purity / immune / stromal host covariates.",
        "python": ["scripts/python/10_validate_phase7a5_host_covariates.py"],
        "r": ["scripts/R/10_phase7a5_host_covariates.R"],
        "report": "analysis_reports/08_host_microbiome_integration/PHASE7A5_HOST_COVARIATES.md",
        "requires_full_data_layout": True,
    },
    "phase7b": {
        "description": "Microbiome-host continuous-state association models (OLS/HC3, PERMANOVA).",
        "python": ["scripts/python/11_summarize_phase7b_associations.py", "scripts/python/11_validate_phase7b_associations.py"],
        "r": ["scripts/R/11_phase7b_microbiome_associations.R"],
        "report": "analysis_reports/08_host_microbiome_integration/PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase8b": {
        "description": "Host-mechanism analysis: Hallmark/PROGENy/DoRothEA/WGCNA for primary taxa.",
        "python": ["scripts/python/13_summarize_phase8b_mechanisms.py", "scripts/python/13_validate_phase8b_mechanisms.py"],
        "r": ["scripts/R/13_phase8b_host_mechanisms.R"],
        "report": "analysis_reports/08_host_microbiome_integration/PHASE8B_HOST_MECHANISM_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase9b1r": {
        "description": "Corrected independent bulk external validation (TCGA-PAAD, GSE71729, GSE62452).",
        "python": ["scripts/python/14_prepare_phase9b1r_bulk_data.py", "scripts/python/14_validate_phase9b1r_bulk_validation.py"],
        "r": ["scripts/R/14_phase9b1r_corrected_bulk_validation.R"],
        "report": "analysis_reports/09_external_validation/PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase9b2r": {
        "description": "Corrected single-cell cellular-source validation (PENG_CRA001160).",
        "python": ["scripts/python/15_prepare_phase9b2r_single_cell.py", "scripts/python/15_validate_phase9b2r_single_cell.py"],
        "r": ["scripts/R/15_phase9b2r_corrected_single_cell_validation.R"],
        "report": "analysis_reports/09_external_validation/PHASE9B2R_CORRECTED_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase9b3r": {
        "description": "Corrected spatial-transcriptomic validation (Hwang GeoMx + Moncada ST).",
        "python": ["scripts/python/16_phase9b3r_spatial_validation.py", "scripts/python/16_validate_phase9b3r_spatial.py"],
        "r": [],
        "report": "analysis_reports/09_external_validation/PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS.md",
        "requires_full_data_layout": True,
    },
    "phase10br": {
        "description": "Corrected cross-layer evidence synthesis and target-prioritisation scoring.",
        "python": ["scripts/python/18_phase10br_cross_layer_synthesis.py", "scripts/python/18_validate_phase10br_synthesis.py"],
        "r": [],
        "report": "analysis_reports/10_target_prioritization/PHASE10BR_CORRECTED_TARGET_PRIORITIZATION_RESULTS.md",
        "requires_full_data_layout": False,  # consumes results/tables/*, included in this release
    },
    "phase11d": {
        "description": "Full manuscript assembly from the locked claim map and Phase 11C review.",
        "python": ["scripts/python/19_validate_phase11d_full_manuscript.py"],
        "r": [],
        "report": "analysis_reports/11_manuscript/PHASE11D_FULL_MANUSCRIPT_DRAFT.md",
        "requires_full_data_layout": False,
    },
    "phase11h": {
        "description": "Submission package assembly (manuscript, figures, tables, cover letter).",
        "python": ["scripts/python/20_validate_phase11h_submission_package.py"],
        "r": [],
        "report": "analysis_reports/11_manuscript/PHASE11H_SUBMISSION_PACKAGE_ASSEMBLY.md",
        "requires_full_data_layout": False,
        "note": "The assembled submission package itself is excluded from this release.",
    },
}


def cmd_list() -> None:
    width = max(len(k) for k in PHASES)
    print(f"{'PHASE ID'.ljust(width)}  DESCRIPTION")
    print("-" * (width + 60))
    for phase_id, info in PHASES.items():
        print(f"{phase_id.ljust(width)}  {info['description']}")
    print(
        "\nRun 'python run_pipeline.py --describe <phase_id>' for scripts and report path.\n"
        "This list is a curated subset of the project's ~50 sub-phases, focused on the\n"
        "primary/corrected run for each major stage. See docs/planning/DECISION_LOG.md and\n"
        "ai_usage/SKILL_USAGE_LOG.tsv for the complete phase-by-phase history, including\n"
        "superseded/rejected attempts (e.g. phase9b1 before its 9B1R correction)."
    )


def cmd_describe(phase_id: str) -> None:
    info = PHASES.get(phase_id)
    if info is None:
        print(f"Unknown phase id: {phase_id!r}. Run --list to see available phases.", file=sys.stderr)
        sys.exit(1)
    print(f"Phase: {phase_id}")
    print(f"Description: {info['description']}")
    if info.get("note"):
        print(f"Note: {info['note']}")
    print(f"Report: {info['report']}")
    if info["python"]:
        print("Python scripts:")
        for p in info["python"]:
            print(f"  - {p}")
    if info["r"]:
        print("R scripts:")
        for r in info["r"]:
            print(f"  - {r}")
    if info["requires_full_data_layout"]:
        print(
            "\nThis phase's scripts read from 02_data/ and/or 03_processed/, which are NOT "
            "included in this release (see README.md). Restore the original project's data "
            "layout before running with --run, or pass --force to attempt it anyway."
        )


def cmd_run(phase_id: str, force: bool) -> None:
    info = PHASES.get(phase_id)
    if info is None:
        print(f"Unknown phase id: {phase_id!r}. Run --list to see available phases.", file=sys.stderr)
        sys.exit(1)

    if info["requires_full_data_layout"] and not force:
        missing = [f"{d}/" for d in ("02_data", "03_processed") if not (REPO_ROOT / d).exists()]
        if missing:
            print(
                f"Phase '{phase_id}' expects {', '.join(missing)} from the full project "
                "checkout, which this release does not include. Nothing was run.\n"
                "Restore those directories from the original project, or re-run with --force "
                "to attempt it anyway (it will likely fail with a file-not-found error)."
            )
            return

    scripts = [(REPO_ROOT / p, "python3") for p in info["python"]] + [
        (REPO_ROOT / r, "Rscript") for r in info["r"]
    ]
    if not scripts:
        print(f"Phase '{phase_id}' has no runnable script in this release; see its report instead.")
        return

    for script_path, interpreter in scripts:
        if not script_path.exists():
            print(f"Expected script not found: {script_path}", file=sys.stderr)
            continue
        print(f"--- Running: {interpreter} {script_path.relative_to(REPO_ROOT)} ---")
        result = subprocess.run([interpreter, str(script_path)], cwd=REPO_ROOT)
        if result.returncode != 0:
            print(f"Script exited with status {result.returncode}: {script_path}", file=sys.stderr)
            sys.exit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--list", action="store_true", help="List all available phases.")
    group.add_argument("--describe", metavar="PHASE_ID", help="Show scripts and report path for a phase.")
    group.add_argument("--run", metavar="PHASE_ID", help="Run a phase's script(s) via subprocess.")
    parser.add_argument("--force", action="store_true", help="With --run: attempt to run even if expected data inputs are missing.")
    args = parser.parse_args()

    if args.list:
        cmd_list()
    elif args.describe:
        cmd_describe(args.describe)
    elif args.run:
        cmd_run(args.run, args.force)


if __name__ == "__main__":
    main()
