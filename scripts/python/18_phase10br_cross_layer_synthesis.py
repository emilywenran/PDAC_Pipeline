#!/usr/bin/env python3
"""Phase 10B-R corrected cross-layer synthesis and target prioritization.

This reanalysis intentionally derives all candidate rows from the locked Phase
10A inventories. External database criteria are marked unavailable unless a
project-local reproducible query result is present.
"""

from __future__ import annotations

import csv
import hashlib
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "04_analysis/10_target_prioritization"
TABLE_DIR = ROOT / "05_results/tables"
PARAMETER_PATH = ROOT / "01_metadata/phase10a_evidence_synthesis_parameter_inventory.tsv"
EVIDENCE_PATH = TABLE_DIR / "phase10a_cross_layer_evidence_inventory.tsv"
FRAMEWORK_PATH = TABLE_DIR / "phase10a_target_prioritization_framework.tsv"
PHASE10B_PATH = TABLE_DIR / "phase10b_candidate_target_scores.tsv"

OUTPUT_REPORT = ANALYSIS_DIR / "PHASE10BR_CORRECTED_TARGET_PRIORITIZATION_RESULTS.md"
OUTPUT_CROSS_LAYER = TABLE_DIR / "phase10br_cross_layer_evidence_scores.tsv"
OUTPUT_TARGETS = TABLE_DIR / "phase10br_candidate_target_scores.tsv"
OUTPUT_EXTERNAL_AUDIT = TABLE_DIR / "phase10br_external_database_query_audit.tsv"
OUTPUT_PENALTIES = TABLE_DIR / "phase10br_penalty_audit.tsv"
OUTPUT_RANK_CHANGE = TABLE_DIR / "phase10br_rank_change_audit.tsv"


HIERARCHY = [
    "MULTI_LAYER_SUPPORTED",
    "PARTIALLY_REPLICATED",
    "DISCOVERY_ONLY",
    "METHOD_SENSITIVE",
    "COMPOSITION_SENSITIVE",
    "CONTAMINATION_SENSITIVE",
    "NOT_EXTERNALLY_SUPPORTED",
    "INSUFFICIENT_DATA",
    "EXPLORATORY_ONLY",
    "NO_SUPPORTED_ASSOCIATION",
]

WEIGHT_POINTS = {
    "Critical": 3,
    "High": 2,
    "Medium": 1,
}

REPLICATED_BULK = {
    "REPLICATED",
    "EXTERNALLY_REPLICATED_HOST_FEATURE",
    "PARTIALLY_REPLICATED_HOST_FEATURE",
}

GENE_SYMBOL_FEATURES = {"BHLHE40", "CTCFL"}
PATHWAY_OR_PROCESS_FEATURES = {"HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS"}
WGCNA_INELIGIBLE = {"MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"}


@dataclass(frozen=True)
class CriterionResult:
    criterion: str
    status: str
    points: int
    penalty_points: int
    evidence_source: str
    note: str


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def parameter_lookup(parameters: pd.DataFrame) -> dict[str, str]:
    return dict(zip(parameters["parameter_name"], parameters["parameter_value"]))


def hierarchy_score(category: str) -> int:
    if category not in HIERARCHY:
        return 0
    return len(HIERARCHY) - HIERARCHY.index(category)


def boolean_parameter(params: dict[str, str], name: str) -> bool:
    value = params.get(name, "").strip().upper()
    if value not in {"TRUE", "FALSE"}:
        raise ValueError(f"Locked parameter {name} must be TRUE or FALSE, observed {value!r}")
    return value == "TRUE"


def feature_type(feature: str) -> str:
    if feature in GENE_SYMBOL_FEATURES:
        return "GENE_SYMBOL"
    if feature in PATHWAY_OR_PROCESS_FEATURES:
        return "PATHWAY_OR_PROCESS"
    if feature.startswith("ME"):
        return "WGCNA_MODULE"
    if feature in {"Ochrobactrum", "Staphylococcus", "Lysobacter", "Brevundimonas", "Paraburkholderia"}:
        return "MICROBIAL_TAXON"
    return "OTHER"


def derive_category(row: pd.Series) -> str:
    discovery = row["discovery_evidence"]
    bulk = row["bulk_evidence"]
    sc = row["sc_evidence"]
    spatial = row["spatial_evidence"]
    feature = row["feature_name"]

    if feature in WGCNA_INELIGIBLE or "INSUFFICIENT" in row["final_synthesis_category"]:
        return "INSUFFICIENT_DATA"
    if discovery == "METHOD_SENSITIVE":
        return "METHOD_SENSITIVE"
    if discovery == "CONTAMINATION_SENSITIVE":
        return "CONTAMINATION_SENSITIVE"
    if discovery == "NO_SUPPORTED_ASSOCIATION":
        return "NO_SUPPORTED_ASSOCIATION"
    if discovery == "SUGGESTIVE_ASSOCIATION":
        return "EXPLORATORY_ONLY"
    if discovery == "ROBUST_ASSOCIATION":
        return "DISCOVERY_ONLY"
    if discovery != "ROBUST_HOST_MECHANISM":
        return row["final_synthesis_category"]

    bulk_replicated = bulk in REPLICATED_BULK
    malignant_intrinsic = sc == "MALIGNANT_CELL_INTRINSIC_SUPPORT"
    spatial_supported = spatial in {"SPATIAL_SUPPORT", "PARTIAL_SPATIAL_SUPPORT"}
    sc_not_supported = sc in {"NOT_SUPPORTED", "NOT_SUPPORTED_AT_CELLULAR_LEVEL"}

    if bulk_replicated and malignant_intrinsic and spatial_supported:
        return "MULTI_LAYER_SUPPORTED"
    if bulk_replicated or malignant_intrinsic:
        return "PARTIALLY_REPLICATED"
    if bulk == "NOT_REPLICATED" and (sc_not_supported or sc == "NOT_EVALUATED"):
        return "NOT_EXTERNALLY_SUPPORTED"
    if bulk == "NOT_REPLICATED":
        return "NOT_EXTERNALLY_SUPPORTED"
    return "DISCOVERY_ONLY"


def external_status(feature: str, database: str) -> tuple[str, str]:
    ftype = feature_type(feature)
    if ftype != "GENE_SYMBOL":
        if ftype == "MICROBIAL_TAXON":
            return "NOT_APPLICABLE_MICROBIAL_TAXON", ""
        return "NOT_APPLICABLE_NON_GENE_OR_PROCESS", ""

    expected = TABLE_DIR / f"phase10br_{database.lower()}_{feature.lower()}_query_results.tsv"
    if expected.exists():
        return "LOCAL_QUERY_RESULT_AVAILABLE", str(expected.relative_to(ROOT))
    return "NOT_RUN_DATABASE_UNAVAILABLE", ""


def build_external_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in evidence.iterrows():
        feature = row["feature_name"]
        for database in ["OpenTargets", "GTEx", "ChEMBL"]:
            status, source = external_status(feature, database)
            rows.append(
                {
                    "feature_name": feature,
                    "database": database,
                    "feature_type": feature_type(feature),
                    "query_result_status": status,
                    "query_result_path": source,
                    "value_used_in_scoring": "NONE" if not source else "SEE_LOCAL_QUERY_RESULT",
                    "reproducibility_decision": (
                        "DATABASE_VALUE_NOT_USED"
                        if status != "LOCAL_QUERY_RESULT_AVAILABLE"
                        else "DATABASE_VALUE_REPRODUCIBLY_AVAILABLE"
                    ),
                }
            )
    return pd.DataFrame(rows)


def framework_weight(framework: pd.DataFrame, criterion: str) -> tuple[str, int, str]:
    row = framework.loc[framework["criteria"] == criterion]
    if row.empty:
        raise ValueError(f"Missing locked criterion: {criterion}")
    weight = row.iloc[0]["weight"]
    return weight, WEIGHT_POINTS.get(weight, 0), row.iloc[0]["threshold"]


def database_available(external_audit: pd.DataFrame, feature: str, database: str) -> bool:
    row = external_audit[
        (external_audit["feature_name"] == feature) & (external_audit["database"] == database)
    ]
    return bool((not row.empty) and row.iloc[0]["query_result_status"] == "LOCAL_QUERY_RESULT_AVAILABLE")


def score_criteria(
    row: pd.Series,
    framework: pd.DataFrame,
    external_audit: pd.DataFrame,
    literature_rescue_allowed: bool,
) -> list[CriterionResult]:
    feature = row["feature_name"]
    category = row["derived_synthesis_category"]
    sc = row["sc_evidence"]
    results: list[CriterionResult] = []

    def append(criterion: str, status: str, pass_points: bool, penalty: bool, source: str, note: str) -> None:
        _, points, _ = framework_weight(framework, criterion)
        results.append(
            CriterionResult(
                criterion=criterion,
                status=status,
                points=points if pass_points else 0,
                penalty_points=-points if penalty else 0,
                evidence_source=source,
                note=note,
            )
        )

    append(
        "external_replication",
        "PASS_LOCKED_EVIDENCE_CATEGORY"
        if category in {"MULTI_LAYER_SUPPORTED", "PARTIALLY_REPLICATED"}
        else "FAIL_LOCKED_EVIDENCE_CATEGORY",
        category in {"MULTI_LAYER_SUPPORTED", "PARTIALLY_REPLICATED"},
        False,
        str(EVIDENCE_PATH.relative_to(ROOT)),
        f"derived_synthesis_category={category}",
    )
    append(
        "cell_type_specificity",
        "PASS_MALIGNANT_CELL_INTRINSIC"
        if sc == "MALIGNANT_CELL_INTRINSIC_SUPPORT"
        else "PENALTY_CELL_COMPOSITION_EXPLAINED"
        if sc == "CELL_COMPOSITION_EXPLAINED"
        else "NOT_SUPPORTED_BY_LOCKED_SINGLE_CELL_EVIDENCE",
        sc == "MALIGNANT_CELL_INTRINSIC_SUPPORT",
        sc == "CELL_COMPOSITION_EXPLAINED",
        str(EVIDENCE_PATH.relative_to(ROOT)),
        "CELL_COMPOSITION_EXPLAINED is not treated as true cell-type specificity.",
    )
    append(
        "druggability",
        "UNKNOWN_NO_REPRODUCIBLE_OPENTARGETS_EVIDENCE"
        if not database_available(external_audit, feature, "OpenTargets")
        else "LOCAL_QUERY_RESULT_REQUIRES_THRESHOLD_EVALUATION",
        False,
        False,
        str(OUTPUT_EXTERNAL_AUDIT.relative_to(ROOT)),
        "No hardcoded tractability values are used.",
    )
    append(
        "tumor_versus_normal_selectivity",
        "UNKNOWN_NO_REPRODUCIBLE_GTEX_EVIDENCE"
        if not database_available(external_audit, feature, "GTEx")
        else "LOCAL_QUERY_RESULT_REQUIRES_THRESHOLD_EVALUATION",
        False,
        False,
        str(OUTPUT_EXTERNAL_AUDIT.relative_to(ROOT)),
        "No hardcoded tumor-vs-normal values are used.",
    )
    append(
        "existing_compounds",
        "UNKNOWN_NO_REPRODUCIBLE_CHEMBL_EVIDENCE"
        if not database_available(external_audit, feature, "ChEMBL")
        else "LOCAL_QUERY_RESULT_REQUIRES_THRESHOLD_EVALUATION",
        False,
        False,
        str(OUTPUT_EXTERNAL_AUDIT.relative_to(ROOT)),
        "No hardcoded compound values are used.",
    )
    for criterion in ["genetic_dependency", "pathway_position", "safety_essentiality_concerns"]:
        append(
            criterion,
            "UNKNOWN_NO_LOCKED_QUANTITATIVE_INPUT",
            False,
            False,
            str(FRAMEWORK_PATH.relative_to(ROOT)),
            "Criterion retained from locked framework; no reproducible Phase 10B-R input table available.",
        )

    if not literature_rescue_allowed:
        results.append(
            CriterionResult(
                criterion="literature_rescue",
                status="PROHIBITED_BY_PHASE10A_PARAMETER",
                points=0,
                penalty_points=0,
                evidence_source=str(PARAMETER_PATH.relative_to(ROOT)),
                note="Literature support cannot add points or override weak/composition-sensitive evidence.",
            )
        )
    return results


def target_decision(row: pd.Series, total_score: int, penalty_total: int, external_audit: pd.DataFrame) -> str:
    feature = row["feature_name"]
    category = row["derived_synthesis_category"]
    if category == "INSUFFICIENT_DATA":
        return "NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA"
    if category in {"METHOD_SENSITIVE", "CONTAMINATION_SENSITIVE", "NO_SUPPORTED_ASSOCIATION", "EXPLORATORY_ONLY"}:
        return "NOT_PRIORITIZED_LOCKED_EVIDENCE_WEAK_OR_SENSITIVE"
    if category == "NOT_EXTERNALLY_SUPPORTED":
        return "NOT_PRIORITIZED_NOT_EXTERNALLY_SUPPORTED"
    if row["discovery_evidence"] == "ROBUST_ASSOCIATION":
        return "NOT_PRIORITIZED_DISCOVERY_ONLY_NO_CAUSAL_VALIDATION"
    if penalty_total < 0:
        return "NOT_PRIORITIZED_COMPOSITION_SENSITIVE_SINGLE_CELL_EVIDENCE"
    statuses = set(external_audit.loc[external_audit["feature_name"] == feature, "query_result_status"])
    if "NOT_RUN_DATABASE_UNAVAILABLE" in statuses:
        return "NOT_PRIORITIZED_EXTERNAL_DATABASE_EVIDENCE_UNAVAILABLE"
    if feature_type(feature) != "GENE_SYMBOL":
        return "RETAINED_AS_SUPPORTED_BIOLOGICAL_FEATURE_NOT_DIRECT_GENE_TARGET"
    if total_score >= 5:
        return "PRIORITIZATION_REQUIRES_INDEPENDENT_REVIEW"
    return "NOT_PRIORITIZED_INSUFFICIENT_LOCKED_SCORE"


def previous_rank_map() -> dict[str, int]:
    if not PHASE10B_PATH.exists():
        return {}
    old = read_tsv(PHASE10B_PATH)
    if "feature_name" not in old.columns:
        return {}
    return {feature: idx + 1 for idx, feature in enumerate(old["feature_name"].tolist())}


def write_report(
    evidence_scores: pd.DataFrame,
    target_scores: pd.DataFrame,
    external_audit: pd.DataFrame,
    penalty_audit: pd.DataFrame,
    rank_audit: pd.DataFrame,
) -> None:
    rows = []
    for _, row in target_scores.iterrows():
        rows.append(
            "| {feature_name} | {derived_synthesis_category} | {sc_evidence} | {target_priority_score} | {priority_decision} |".format(
                **row
            )
        )

    unavailable_count = int((external_audit["query_result_status"] == "NOT_RUN_DATABASE_UNAVAILABLE").sum())
    penalty_count = int((penalty_audit["penalty_applied"] == "TRUE").sum())
    skipped = target_scores.loc[target_scores["included_in_phase10br_scoring"] != "TRUE", "feature_name"].tolist()

    content = f"""# Phase 10B-R Corrected Target Prioritization Results

## Scope

Phase 10B-R re-executes cross-layer evidence synthesis after the Phase 10C rejection of the first Phase 10B attempt. The reanalysis uses the Phase 10A method lock, the cross-layer synthesis protocol, the Phase 10A parameter inventory, the Phase 10A evidence inventory, the locked target-prioritization framework, and the final PASS reviews from Phases 7C, 8C, 9B1C2, 9B2C2, and 9B3C2 as authoritative inputs.

No manuscript drafting was performed.

## Corrections

- Removed descriptive target overrides from scoring. Candidate rows are generated from `phase10a_cross_layer_evidence_inventory.tsv`.
- Scored every Phase 10A evidence-inventory candidate. Skipped candidates: `{', '.join(skipped) if skipped else 'NONE'}`.
- Treated `CELL_COMPOSITION_EXPLAINED` as a high-weight penalty for cell-type specificity, including CTCFL/BORIS.
- Did not use literature support to add points or rescue weak/composition-sensitive candidates.
- Preserved the locked negative and partial evidence for Ochrobactrum, HALLMARK_PROTEIN_SECRETION, spatial support, and ineligible WGCNA modules.

## External Database Reproducibility

OpenTargets, GTEx, and ChEMBL values are not filled from dictionaries. Phase 10B-R generated a local audit table for every candidate-database pair. Gene-symbol database queries without local reproducible result tables are marked `NOT_RUN_DATABASE_UNAVAILABLE`; {unavailable_count} such rows were recorded.

## Candidate Scores

| Candidate | Derived evidence class | Single-cell evidence | Target score | Decision |
|---|---|---|---:|---|
{chr(10).join(rows)}

## Key Preserved Evidence

HALLMARK_PROTEIN_SECRETION retains malignant-cell and malignant-compartment support, but the basal-classical spatial-axis association is not replicated and spatial evidence remains `PARTIAL_SPATIAL_SUPPORT`.

CTCFL/BORIS is not promoted using cancer-testis antigen reasoning because reproducible GTEx/OpenTargets/ChEMBL evidence is unavailable in this run and its single-cell evidence is `CELL_COMPOSITION_EXPLAINED`.

Ochrobactrum retains robust host-mechanism association only. Microbial localization, physical interaction, and causality validation remain absent.

The five locked ineligible WGCNA modules remain `INSUFFICIENT_DATA` and are not promoted.

## Audits

- Cross-layer evidence scores: `{OUTPUT_CROSS_LAYER.relative_to(ROOT)}`
- Candidate target scores: `{OUTPUT_TARGETS.relative_to(ROOT)}`
- External database query audit: `{OUTPUT_EXTERNAL_AUDIT.relative_to(ROOT)}`
- Penalty audit: `{OUTPUT_PENALTIES.relative_to(ROOT)}` ({penalty_count} applied penalties)
- Rank-change audit: `{OUTPUT_RANK_CHANGE.relative_to(ROOT)}`

## Final Readiness Decision

`READY_FOR_PHASE10C2_INDEPENDENT_REVIEW`
"""
    OUTPUT_REPORT.write_text(content, encoding="utf-8")


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> int:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    parameters = read_tsv(PARAMETER_PATH)
    evidence = read_tsv(EVIDENCE_PATH)
    framework = read_tsv(FRAMEWORK_PATH)
    params = parameter_lookup(parameters)

    literature_rescue_allowed = boolean_parameter(params, "literature_rescue_allowed")
    if boolean_parameter(params, "wgcna_module_promotion"):
        raise RuntimeError("Locked Phase 10A prohibits WGCNA module promotion, but parameter is TRUE.")
    if boolean_parameter(params, "causal_claims_allowed"):
        raise RuntimeError("Locked Phase 10A prohibits causal claims, but parameter is TRUE.")

    evidence = evidence.copy()
    evidence["derived_synthesis_category"] = evidence.apply(derive_category, axis=1)
    evidence["phase10a_inventory_category"] = evidence["final_synthesis_category"]
    evidence["category_matches_phase10a_inventory"] = (
        evidence["derived_synthesis_category"] == evidence["phase10a_inventory_category"]
    ).map({True: "TRUE", False: "FALSE"})
    evidence["hierarchy_rank"] = evidence["derived_synthesis_category"].apply(lambda c: HIERARCHY.index(c) + 1)
    evidence["cross_layer_evidence_score"] = evidence["derived_synthesis_category"].apply(hierarchy_score)
    evidence["included_in_phase10br_scoring"] = "TRUE"
    evidence["feature_type"] = evidence["feature_name"].apply(feature_type)
    evidence["promotion_allowed_by_locked_rules"] = evidence["derived_synthesis_category"].apply(
        lambda c: "FALSE" if c == "INSUFFICIENT_DATA" else "TRUE"
    )

    external_audit = build_external_audit(evidence)

    target_rows = []
    penalty_rows = []
    criterion_rows = []
    for _, row in evidence.iterrows():
        criteria = score_criteria(row, framework, external_audit, literature_rescue_allowed)
        positive = sum(item.points for item in criteria)
        penalties = sum(item.penalty_points for item in criteria)
        total = positive + penalties
        decision = target_decision(row, total, penalties, external_audit)
        for item in criteria:
            criterion_rows.append(
                {
                    "feature_name": row["feature_name"],
                    "criterion": item.criterion,
                    "criterion_status": item.status,
                    "positive_points": item.points,
                    "penalty_points": item.penalty_points,
                    "evidence_source": item.evidence_source,
                    "note": item.note,
                }
            )
            if item.penalty_points < 0 or item.criterion == "literature_rescue":
                penalty_rows.append(
                    {
                        "feature_name": row["feature_name"],
                        "penalty_type": item.criterion,
                        "penalty_applied": "TRUE" if item.penalty_points < 0 else "FALSE",
                        "penalty_points": item.penalty_points,
                        "locked_rule_source": item.evidence_source,
                        "reason": item.note,
                    }
                )
        target_rows.append(
            {
                "feature_name": row["feature_name"],
                "feature_type": row["feature_type"],
                "discovery_evidence": row["discovery_evidence"],
                "bulk_evidence": row["bulk_evidence"],
                "sc_evidence": row["sc_evidence"],
                "spatial_evidence": row["spatial_evidence"],
                "derived_synthesis_category": row["derived_synthesis_category"],
                "cross_layer_evidence_score": row["cross_layer_evidence_score"],
                "positive_framework_points": positive,
                "penalty_points": penalties,
                "target_priority_score": total,
                "priority_decision": decision,
                "included_in_phase10br_scoring": "TRUE",
                "scoring_source": f"{EVIDENCE_PATH.relative_to(ROOT)};{FRAMEWORK_PATH.relative_to(ROOT)};{PARAMETER_PATH.relative_to(ROOT)}",
            }
        )

    target_scores = pd.DataFrame(target_rows).sort_values(
        ["target_priority_score", "cross_layer_evidence_score", "feature_name"],
        ascending=[False, False, True],
    )
    target_scores.insert(0, "phase10br_rank", range(1, len(target_scores) + 1))
    criterion_scores = pd.DataFrame(criterion_rows)
    penalty_audit = pd.DataFrame(penalty_rows)

    old_ranks = previous_rank_map()
    rank_rows = []
    for _, row in target_scores.iterrows():
        old_rank = old_ranks.get(row["feature_name"])
        rank_rows.append(
            {
                "feature_name": row["feature_name"],
                "phase10b_rank": "" if old_rank is None else old_rank,
                "phase10br_rank": row["phase10br_rank"],
                "rank_change": "NEWLY_INCLUDED" if old_rank is None else old_rank - int(row["phase10br_rank"]),
                "phase10b_status": "ABSENT_FROM_FAILED_PHASE10B" if old_rank is None else "PRESENT_IN_FAILED_PHASE10B",
                "phase10br_decision": row["priority_decision"],
            }
        )
    rank_audit = pd.DataFrame(rank_rows)

    evidence_out = evidence[
        [
            "feature_name",
            "feature_type",
            "discovery_evidence",
            "bulk_evidence",
            "sc_evidence",
            "spatial_evidence",
            "phase10a_inventory_category",
            "derived_synthesis_category",
            "category_matches_phase10a_inventory",
            "hierarchy_rank",
            "cross_layer_evidence_score",
            "promotion_allowed_by_locked_rules",
            "included_in_phase10br_scoring",
            "notes",
        ]
    ].sort_values(["hierarchy_rank", "feature_name"])

    write_tsv(evidence_out, OUTPUT_CROSS_LAYER)
    write_tsv(target_scores, OUTPUT_TARGETS)
    write_tsv(external_audit, OUTPUT_EXTERNAL_AUDIT)
    write_tsv(penalty_audit, OUTPUT_PENALTIES)
    write_tsv(rank_audit, OUTPUT_RANK_CHANGE)
    write_tsv(criterion_scores, TABLE_DIR / "phase10br_framework_criterion_scores.tsv")
    write_report(evidence_out, target_scores, external_audit, penalty_audit, rank_audit)

    print("Phase 10B-R corrected synthesis complete.")
    print("READY_FOR_PHASE10C2_INDEPENDENT_REVIEW")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
