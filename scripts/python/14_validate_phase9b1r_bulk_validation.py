#!/usr/bin/env python3
from pathlib import Path
import sys
import pandas as pd

ROOT = Path("/Users/emily/thesis/PDAC")
TABLE = ROOT / "05_results/tables"
RSCRIPT = ROOT / "06_scripts/R/14_phase9b1r_corrected_bulk_validation.R"

REQUIRED = [
    "phase9b1r_purist_runtime_validation.tsv",
    "phase9b1r_module_transfer_coverage.tsv",
    "phase9b1r_hallmark_runtime_validation.tsv",
    "phase9b1r_hallmark_scores.tsv.gz",
    "phase9b1r_tf_runtime_validation.tsv",
    "phase9b1r_tf_activity_scores.tsv.gz",
    "phase9b1r_cohort_replication_results.tsv",
    "phase9b1r_module_replication_results.tsv",
    "phase9b1r_negative_control_results.tsv",
    "phase9b1r_cross_cohort_synthesis.tsv",
    "phase9b1r_host_feature_replication_evidence.tsv",
]

EXPECTED_TF_COUNTS = {
    "EXTERNALLY_REPLICATED_HOST_FEATURE": 12,
    "PARTIALLY_REPLICATED_HOST_FEATURE": 13,
    "NOT_REPLICATED": 9,
    "TO_VERIFY": 0,
}


def fail(msg):
    raise SystemExit(f"Phase 9B1R validation failed: {msg}")


def classify_tf_feature(group: pd.DataFrame) -> str:
    eligible = len(group)
    supported = int((group["replication_status"] == "SUPPORTED").sum())
    not_supported = int((group["replication_status"] == "NOT_SUPPORTED").sum())
    if eligible == 0:
        return "INSUFFICIENT_EXTERNAL_DATA"
    if supported >= 2:
        return "EXTERNALLY_REPLICATED_HOST_FEATURE"
    if supported == 1:
        return "PARTIALLY_REPLICATED_HOST_FEATURE"
    if supported == 0 and not_supported == 0:
        return "PARTIALLY_REPLICATED_HOST_FEATURE"
    return "NOT_REPLICATED"


def main() -> int:
    for name in REQUIRED:
        path = TABLE / name
        if not path.exists() or path.stat().st_size == 0:
            fail(f"missing required output {name}")

    text = RSCRIPT.read_text()
    forbidden = ["TF_PROXY_", "single-gene expression", "15-gene proxy", "hallmark_proxy"]
    for token in forbidden:
        if token in text:
            fail(f"forbidden proxy implementation token remains in R script: {token}")

    pur = pd.read_csv(TABLE / "phase9b1r_purist_runtime_validation.tsv", sep="\t")
    if not pur["intercept_included"].all() or not (pur["intercept_value"].round(3) == -6.815).all():
        fail("PurIST intercept is absent or incorrect")
    if not ((pur["probability_min"] >= 0) & (pur["probability_max"] <= 1)).all():
        fail("PurIST probabilities outside [0, 1]")
    if (pur["probability_sd"] <= 0).any():
        fail("PurIST probabilities have no sample variation")

    cov = pd.read_csv(TABLE / "phase9b1r_module_transfer_coverage.tsv", sep="\t")
    allowed = {"ELIGIBLE", "INELIGIBLE_LOW_COVERAGE", "TO_VERIFY"}
    if not set(cov["eligibility_status"]).issubset(allowed):
        fail("module coverage eligibility has nonlocked values")
    if (cov.loc[cov["coverage_fraction"] < 0.80, "eligibility_status"] == "ELIGIBLE").any():
        fail("module coverage threshold >=80% is not enforced")

    mod = pd.read_csv(TABLE / "phase9b1r_module_replication_results.tsv", sep="\t")
    bad = mod[(mod["gene_coverage"] < 0.80) & (mod["eligible_for_validation"].astype(str).str.upper() == "TRUE")]
    if len(bad):
        fail("ineligible low-coverage modules entered formal replication")

    neg = pd.read_csv(TABLE / "phase9b1r_negative_control_results.tsv", sep="\t")
    required_controls = {
        "size-matched randomized module gene sets",
        "expression-matched randomized module gene sets",
        "gene-label permutation",
        "unrelated Hallmark pathway",
        "patient-label permutation",
    }
    if not required_controls.issubset(set(neg["control_type"])):
        fail("not all required negative controls were executed")
    if neg["status"].eq("EXECUTED").sum() == 0:
        fail("negative controls are not marked executed")

    repl = pd.read_csv(TABLE / "phase9b1r_cohort_replication_results.tsv", sep="\t")
    for col in ["cohort", "feature_layer", "feature_name", "gene_coverage", "eligible_for_validation",
                "coefficient", "ci_low", "ci_high", "p_value", "q_value", "discovery_direction",
                "external_direction", "replication_status", "exclusion_reason", "notes"]:
        if col not in repl.columns:
            fail(f"cohort replication table missing {col}")
    if repl.groupby(["cohort", "feature_layer"])["q_value"].transform("size").isna().any():
        fail("FDR families could not be evaluated")

    evidence = pd.read_csv(TABLE / "phase9b1r_host_feature_replication_evidence.tsv", sep="\t")
    allowed_evidence = {
        "EXTERNALLY_REPLICATED_HOST_FEATURE",
        "PARTIALLY_REPLICATED_HOST_FEATURE",
        "NOT_REPLICATED",
        "INSUFFICIENT_EXTERNAL_DATA",
        "TO_VERIFY",
    }
    if not set(evidence["evidence_category"]).issubset(allowed_evidence):
        fail("evidence categories are outside Phase 9A rules")
    tf = evidence[evidence["feature_layer"].eq("tf_activity")]
    if len(tf) == 0:
        fail("TF evidence table is missing")
    if tf["evidence_category"].eq("TO_VERIFY").all():
        fail("all eligible TFs were assigned TO_VERIFY without statistical evaluation")
    tf_runtime = pd.read_csv(TABLE / "phase9b1r_tf_runtime_validation.tsv", sep="\t")
    if not tf_runtime["activity_calculation_status"].astype(str).str.contains("EXECUTED").any():
        fail("VIPER activity was not executed for any TF cohort")
    eligible_tf_runtime = tf_runtime[tf_runtime["eligibility"].eq("ELIGIBLE")]
    if len(eligible_tf_runtime) and eligible_tf_runtime["activity_calculation_status"].astype(str).ne("EXECUTED").any():
        fail("TF-symbol expression or other proxy logic appears to have replaced VIPER activity for at least one adequately covered TF")
    repl_tf = repl[repl["feature_layer"].eq("tf_activity")].copy()
    derived_rows = []
    for feature_name, group in repl_tf.groupby("feature_name", sort=True):
        derived_rows.append({
            "feature_name": feature_name,
            "derived_category": classify_tf_feature(group),
            "eligible_cohorts": int(group["eligible_for_validation"].sum()),
            "supported_cohorts": int((group["replication_status"] == "SUPPORTED").sum()),
        })
    derived = pd.DataFrame(derived_rows)
    observed = tf[["feature_name", "evidence_category"]].drop_duplicates()
    merged = derived.merge(observed, on="feature_name", how="left", validate="one_to_one")
    if merged["evidence_category"].isna().any():
        missing = merged.loc[merged["evidence_category"].isna(), "feature_name"].tolist()
        fail(f"TF evidence table is missing categories for: {', '.join(missing)}")
    if (merged["derived_category"] != merged["evidence_category"]).any():
        bad = merged.loc[merged["derived_category"] != merged["evidence_category"], ["feature_name", "derived_category", "evidence_category"]]
        fail("executor-derived TF category counts disagree with the locked evidence rules: " + "; ".join(
            f"{r.feature_name}:{r.derived_category}!={r.evidence_category}" for r in bad.itertuples(index=False)
        ))
    counts = merged["derived_category"].value_counts().to_dict()
    for cat, expected in EXPECTED_TF_COUNTS.items():
        observed_count = int(counts.get(cat, 0))
        if observed_count != expected:
            fail(f"locked TF count mismatch for {cat}: expected {expected}, observed {observed_count}")
    if tf["evidence_category"].eq("TO_VERIFY").any():
        fail("successfully calculated and adequately covered TFs were left TO_VERIFY without a documented reason")

    print("Phase 9B1R validator PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
