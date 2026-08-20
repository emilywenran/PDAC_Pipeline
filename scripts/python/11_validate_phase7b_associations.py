#!/usr/bin/env python3
"""Validate Phase 7B locked microbiome association outputs."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "05_results" / "tables"
REPORT = ROOT / "04_analysis" / "08_host_microbiome_integration" / "PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md"


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def main() -> int:
    failures: list[str] = []
    runtime = pd.read_csv(TABLE / "phase7b_runtime_validation.tsv", sep="\t")
    primary = pd.read_csv(TABLE / "phase7b_primary_genus_associations.tsv", sep="\t")
    cov = pd.read_csv(TABLE / "phase7b_covariate_model_sensitivity.tsv", sep="\t")
    secondary = pd.read_csv(TABLE / "phase7b_secondary_outcome_associations.tsv", sep="\t")
    preproc = pd.read_csv(TABLE / "phase7b_preprocessing_sensitivity.tsv", sep="\t")
    maaslin = pd.read_csv(TABLE / "phase7b_primary_maaslin2.tsv", sep="\t")
    presence = pd.read_csv(TABLE / "phase7b_presence_absence_associations.tsv", sep="\t")
    report = REPORT.read_text()

    require(runtime["passed"].all(), "runtime validation contains failed checks", failures)
    require((runtime.loc[runtime["validation_check"] == "exactly_62_unique_patients", "passed"] == True).all(), "primary patient count is not validated as 62", failures)
    require(len(primary) == 122, f"primary genus rows != 122 ({len(primary)})", failures)
    require(primary["genus"].nunique() == 122, "primary genus names are not unique", failures)
    require(primary["model"].eq("Model_0").all(), "primary table contains non-Model 0 rows", failures)
    require(primary["bh_q_value"].notna().all(), "primary BH q values missing", failures)
    require(cov["model"].isin(["Model_0", "Model_1", "Model_3P", "Model_3I", "Model_3S"]).all(), "unexpected covariate model present", failures)
    require(not cov["model"].eq("Model_2").any(), "clinical Model 2 results were generated", failures)
    require(not cov["formula"].str.contains("estimate_score", na=False).any(), "ESTIMATE score was included as a covariate model", failures)
    require(cov.groupby("model")["host_score_q_value_within_model"].apply(lambda s: s.notna().all()).all(), "covariate q values are not complete within model families", failures)
    require(secondary.groupby("host_outcome")["genus"].nunique().eq(122).all(), "secondary outcome families do not each contain 122 genera", failures)
    require(secondary.groupby("host_outcome")["bh_q_value"].apply(lambda s: s.notna().all()).all(), "secondary q values missing within outcome families", failures)
    require(maaslin["normalization"].eq("NONE").all() and maaslin["transform"].eq("NONE").all(), "MaAsLin2 locked NONE/NONE settings not preserved", failures)
    require(not maaslin["status"].str.contains("NORMALIZATION|TRANSFORMED", case=False, na=False).any(), "MaAsLin2 table suggests second normalization/transform", failures)
    require(preproc.loc[preproc["analysis_id"].isin(["MICRO_SENS_NO_HIGH_RISK", "MICRO_SENS_NO_CONTAMINANTS"]), "matrix_path"].str.contains("sensitivity/").all(), "contaminant sensitivity did not use locked sensitivity matrices", failures)
    require(presence.loc[~presence["eligible_locked_prevalence"], "p_value"].isna().all(), "presence/absence model was applied to ineligible genera", failures)

    causal_patterns = [
        r"\bcauses?\b",
        r"\bcaused\b",
        r"\bdrives?\b",
        r"\bmediates?\b",
        r"\binduces?\b",
    ]
    for pat in causal_patterns:
        require(re.search(pat, report, flags=re.IGNORECASE) is None, f"causal language found in report: {pat}", failures)

    expected_figs = [
        "phase7b_global_community_effects.pdf",
        "phase7b_primary_genus_effects.pdf",
        "phase7b_primary_pvalue_qvalue_plot.pdf",
        "phase7b_primary_genus_scatterplots.pdf",
        "phase7b_covariate_sensitivity_forest.pdf",
        "phase7b_preprocessing_sensitivity_heatmap.pdf",
        "phase7b_contamination_sensitivity.pdf",
        "phase7b_influence_summary.pdf",
        "phase7b_secondary_outcome_concordance.pdf",
        "phase7b_public_subtype_descriptive.pdf",
    ]
    for fig in expected_figs:
        path = ROOT / "05_results" / "figures" / fig
        require(path.exists() and path.stat().st_size > 1000, f"missing or empty figure: {fig}", failures)

    if failures:
        print("VALIDATION_FAIL")
        for f in failures:
            print(f"- {f}")
        return 1
    print("VALIDATION_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
