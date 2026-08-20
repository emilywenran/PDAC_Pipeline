#!/usr/bin/env python3
"""Validate Phase 10B-R corrected synthesis outputs."""

from __future__ import annotations

import ast
import csv
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "04_analysis/10_target_prioritization/PHASE10BR_CORRECTED_TARGET_PRIORITIZATION_RESULTS.md"
EVIDENCE_IN = ROOT / "05_results/tables/phase10a_cross_layer_evidence_inventory.tsv"
SCRIPT = ROOT / "06_scripts/python/18_phase10br_cross_layer_synthesis.py"
CROSS_LAYER = ROOT / "05_results/tables/phase10br_cross_layer_evidence_scores.tsv"
TARGETS = ROOT / "05_results/tables/phase10br_candidate_target_scores.tsv"
EXTERNAL = ROOT / "05_results/tables/phase10br_external_database_query_audit.tsv"
PENALTIES = ROOT / "05_results/tables/phase10br_penalty_audit.tsv"
RANKS = ROOT / "05_results/tables/phase10br_rank_change_audit.tsv"

REQUIRED_OUTPUTS = [REPORT, CROSS_LAYER, TARGETS, EXTERNAL, PENALTIES, RANKS]
FORBIDDEN_REPORT_STRINGS = {
    "READY_FOR_MANUSCRIPT_DRAFTING",
    "most viable biological target",
    "Cancer-Testis Antigen profile",
    "exquisite tumor-versus-normal selectivity",
}
FORBIDDEN_SCRIPT_STRINGS = {
    "target_data =",
    "Cancer-Testis Antigen",
    "GTEx TPM >7",
    "High expression in normal tissues",
    "\"priority_score\"",
    "'priority_score'",
}


class ValidationError(RuntimeError):
    pass


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def fail(errors: list[str]) -> int:
    print("Phase 10B-R validation FAILED:", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    return 1


def validate_no_hardcoded_descriptive_overrides(errors: list[str]) -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    report = REPORT.read_text(encoding="utf-8")
    for text in FORBIDDEN_SCRIPT_STRINGS:
        if text in source:
            errors.append(f"Forbidden hardcoded descriptive override found in script: {text}")
    for text in FORBIDDEN_REPORT_STRINGS:
        if text in report:
            errors.append(f"Forbidden failed-Phase10B narrative/manuscript string found in report: {text}")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {"target_data", "priority_overrides"}:
                    errors.append(f"Forbidden override assignment found: {target.id}")


def validate_all_candidates_scored(errors: list[str]) -> None:
    inventory = read_tsv(EVIDENCE_IN)
    targets = read_tsv(TARGETS)
    missing = sorted(set(inventory["feature_name"]) - set(targets["feature_name"]))
    if missing:
        errors.append("Candidates skipped despite Phase 10A inventory membership: " + ", ".join(missing))
    if (targets["included_in_phase10br_scoring"] != "TRUE").any():
        errors.append("At least one candidate is not marked included_in_phase10br_scoring=TRUE")


def validate_external_database_reproducibility(errors: list[str]) -> None:
    external = read_tsv(EXTERNAL)
    required = {"feature_name", "database", "query_result_status", "query_result_path", "value_used_in_scoring"}
    if not required.issubset(external.columns):
        errors.append(f"External audit missing columns: {sorted(required - set(external.columns))}")
        return
    allowed_no_value = {
        "NOT_RUN_DATABASE_UNAVAILABLE",
        "NOT_APPLICABLE_NON_GENE_OR_PROCESS",
        "NOT_APPLICABLE_MICROBIAL_TAXON",
    }
    for _, row in external.iterrows():
        status = row["query_result_status"]
        value = row["value_used_in_scoring"]
        path = row["query_result_path"]
        if status == "LOCAL_QUERY_RESULT_AVAILABLE":
            if not path:
                errors.append(f"External database row lacks query result path: {row.to_dict()}")
            elif not (ROOT / path).exists():
                errors.append(f"External database query result path does not exist: {path}")
            if value == "NONE":
                errors.append(f"External database row has local result but no value marker: {row.to_dict()}")
        elif status in allowed_no_value:
            if value != "NONE":
                errors.append(f"External database field filled without reproducible evidence: {row.to_dict()}")
        else:
            errors.append(f"Unknown external database status: {status}")


def validate_cell_composition_penalty(errors: list[str]) -> None:
    targets = read_tsv(TARGETS)
    penalties = read_tsv(PENALTIES)
    comp_targets = targets[targets["sc_evidence"] == "CELL_COMPOSITION_EXPLAINED"]
    for feature in comp_targets["feature_name"]:
        row = penalties[
            (penalties["feature_name"] == feature)
            & (penalties["penalty_type"] == "cell_type_specificity")
            & (penalties["penalty_applied"] == "TRUE")
        ]
        if row.empty:
            errors.append(f"CELL_COMPOSITION_EXPLAINED was not penalized for {feature}")
    ctcfl = targets[targets["feature_name"] == "CTCFL"]
    if ctcfl.empty:
        errors.append("CTCFL missing from target scores")
    else:
        row = ctcfl.iloc[0]
        if row["sc_evidence"] != "CELL_COMPOSITION_EXPLAINED":
            errors.append("CTCFL single-cell evidence was altered from CELL_COMPOSITION_EXPLAINED")
        if "COMPOSITION_SENSITIVE" not in row["priority_decision"]:
            errors.append("CTCFL/BORIS was not blocked by composition-sensitive single-cell evidence")


def validate_literature_no_rescue(errors: list[str]) -> None:
    targets = read_tsv(TARGETS)
    penalties = read_tsv(PENALTIES)
    lit = penalties[penalties["penalty_type"] == "literature_rescue"]
    if lit.empty:
        errors.append("Literature rescue prohibition audit rows are missing")
    rescued = targets[
        (targets["sc_evidence"] == "CELL_COMPOSITION_EXPLAINED")
        & (~targets["priority_decision"].str.contains("COMPOSITION_SENSITIVE", regex=False))
    ]
    if not rescued.empty:
        errors.append("Composition-sensitive candidate appears rescued: " + ", ".join(rescued["feature_name"]))


def validate_preserved_evidence(errors: list[str]) -> None:
    cross = read_tsv(CROSS_LAYER)
    targets = read_tsv(TARGETS)
    hps = cross[cross["feature_name"] == "HALLMARK_PROTEIN_SECRETION"]
    if hps.empty or hps.iloc[0]["spatial_evidence"] != "PARTIAL_SPATIAL_SUPPORT":
        errors.append("HALLMARK_PROTEIN_SECRETION partial spatial support was not preserved")
    och = targets[targets["feature_name"] == "Ochrobactrum"]
    if och.empty or "DISCOVERY_ONLY" not in och.iloc[0]["derived_synthesis_category"]:
        errors.append("Ochrobactrum discovery-only status was not preserved")
    ineligible = targets[targets["feature_name"].isin(["MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"])]
    bad = ineligible[ineligible["priority_decision"] != "NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA"]
    if not bad.empty:
        errors.append("Ineligible WGCNA modules were not kept ineligible: " + ", ".join(bad["feature_name"]))


def validate_no_manuscript_drafting(errors: list[str]) -> None:
    report = REPORT.read_text(encoding="utf-8")
    if "READY_FOR_PHASE10C2_INDEPENDENT_REVIEW" not in report:
        errors.append("Final readiness decision missing")
    if "READY_FOR_MANUSCRIPT_DRAFTING" in report:
        errors.append("Phase 10B-R proceeds directly to manuscript drafting")


def main() -> int:
    errors: list[str] = []
    for path in REQUIRED_OUTPUTS + [SCRIPT, EVIDENCE_IN]:
        if not path.exists():
            errors.append(f"Required file missing: {path.relative_to(ROOT)}")
    if errors:
        return fail(errors)

    validate_no_hardcoded_descriptive_overrides(errors)
    validate_all_candidates_scored(errors)
    validate_external_database_reproducibility(errors)
    validate_cell_composition_penalty(errors)
    validate_literature_no_rescue(errors)
    validate_preserved_evidence(errors)
    validate_no_manuscript_drafting(errors)

    if errors:
        return fail(errors)
    print("Phase 10B-R validation passed: READY_FOR_PHASE10C2_INDEPENDENT_REVIEW")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
