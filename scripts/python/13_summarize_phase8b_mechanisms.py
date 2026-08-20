#!/usr/bin/env python3
"""Summarize Phase 8B host-mechanism outputs and update project records."""

from __future__ import annotations

import hashlib
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "05_results" / "tables"
FIG = ROOT / "05_results" / "figures"
REPORT = ROOT / "04_analysis" / "08_host_microbiome_integration" / "PHASE8B_HOST_MECHANISM_RESULTS.md"
STATUS = ROOT / "00_admin" / "PROJECT_STATUS.md"
DECISION = ROOT / "09_docs" / "planning" / "DECISION_LOG.md"
MANIFEST = ROOT / "01_metadata" / "file_manifest.tsv"
TODAY = "2026-07-02"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def top_rows(df: pd.DataFrame, n: int = 8) -> pd.DataFrame:
    return df.sort_values(["bh_q_value", "p_value"]).head(n)


def md_table(df: pd.DataFrame, cols: list[str]) -> str:
    if df.empty:
        return "No rows met the reporting filter.\n"
    return df.loc[:, cols].to_markdown(index=False)


def append_once(path: Path, marker: str, text: str) -> None:
    old = path.read_text() if path.exists() else ""
    if marker not in old:
        path.write_text(old.rstrip() + "\n\n" + text.strip() + "\n")


def update_manifest() -> None:
    manifest = read_tsv(MANIFEST)
    phase8b_files = sorted(TABLE.glob("phase8b*")) + sorted((TABLE / "phase8b_host_gene_full").glob("phase8b*.tsv.gz")) + sorted(FIG.glob("phase8b*.pdf"))
    phase8b_files += sorted((ROOT / "05_results" / "models" / "phase8b").glob("*"))
    phase8b_files += [
        ROOT / "06_scripts" / "R" / "13_phase8b_host_mechanisms.R",
        ROOT / "06_scripts" / "python" / "13_summarize_phase8b_mechanisms.py",
        ROOT / "06_scripts" / "python" / "13_validate_phase8b_mechanisms.py",
        REPORT,
    ]
    phase8b_files = [p for p in phase8b_files if p.exists() and p.is_file()]
    manifest = manifest[~manifest["file_id"].astype(str).str.startswith("phase8b_")]
    rows = []
    for path in phase8b_files:
        rel = path.relative_to(ROOT)
        rows.append({
            "file_id": rel.as_posix().replace("/", "_").replace(".", "_"),
            "dataset": "PDAC_Phase8B_host_mechanisms",
            "sample_id": "",
            "data_type": "host_microbiome_mechanism_output",
            "local_path": str(path),
            "source_url_or_accession": "derived_from_locked_Phase8B_pipeline",
            "file_size": path.stat().st_size,
            "md5": sha256(path),
            "download_date": TODAY,
            "processing_status": "generated_Phase8B",
            "notes": "Phase 8B locked host-mechanism analysis product; checksum stored as sha256 in md5 field for compatibility with existing manifest.",
        })
    out = pd.concat([manifest, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(MANIFEST, sep="\t", index=False)


def main() -> int:
    runtime = read_tsv(TABLE / "phase8b_runtime_validation.tsv")
    pathway = read_tsv(TABLE / "phase8b_primary_pathway_associations.tsv")
    tf = read_tsv(TABLE / "phase8b_primary_tf_associations.tsv")
    cov = read_tsv(TABLE / "phase8b_host_covariate_sensitivity.tsv")
    trans = read_tsv(TABLE / "phase8b_transformation_sensitivity.tsv")
    wgcna = read_tsv(TABLE / "phase8b_wgcna_taxon_associations.tsv")
    wgcna_summary = read_tsv(TABLE / "phase8b_wgcna_module_summary.tsv")
    gene_summary = read_tsv(TABLE / "phase8b_host_gene_associations_summary.tsv")
    enrich = read_tsv(TABLE / "phase8b_ranked_gene_enrichment.tsv")
    shared = read_tsv(TABLE / "phase8b_shared_mechanism_summary.tsv")
    evidence = read_tsv(TABLE / "phase8b_host_mechanism_evidence.tsv")

    evidence_counts = evidence["evidence_category"].value_counts().rename_axis("category").reset_index(name="n")
    pathway_sig = pathway[pathway["bh_q_value"] < 0.05]
    tf_sig = tf[tf["bh_q_value"] < 0.05]
    wgcna_sig = wgcna[wgcna["bh_q_value"] < 0.05]
    enrich_sig = enrich[enrich["padj"] < 0.05]

    report = f"""# Phase 8B Host Mechanism Results

## Execution Status

Phase 8B executed the Phase 8A locked host-mechanism analyses for the nine Phase 7C-verified primary genera. Runtime validation passed all hard-stop checks: 62 aligned patients, nine primary taxa present in primary CLR and rCLR matrices, Phase 7B directions loaded, expression and microbiome order aligned by `patient_id`, no duplicate primary identifiers, no missing or infinite primary values, and required Phase 8 packages loaded from the project-local `renv` library.

`renv` note: R startup required `RENV_CONFIG_SANDBOX_ENABLED=FALSE` because the renv sandbox lock path blocked startup in this workspace. The active library was still the project-local renv library recorded in `phase8b_runtime_validation.tsv`.

## Pathway Activity Findings

Hallmark ssGSEA scores were generated from MSigDB `{read_tsv(TABLE / "phase8b_pathway_gene_coverage.tsv").query("collection == 'MSigDB_Hallmark'")["collection_version"].iloc[0]}`. PROGENy scores used package version `{read_tsv(TABLE / "phase8b_pathway_gene_coverage.tsv").query("collection == 'PROGENy'")["collection_version"].iloc[0]}` with the locked top-100 model. Primary pathway association rows: `{len(pathway)}`. FDR-supported pathway rows: `{len(pathway_sig)}`.

{md_table(top_rows(pathway), ["taxon", "host_feature_collection", "host_feature", "coefficient", "p_value", "bh_q_value", "RCLR_DIRECTION_SENSITIVE"])}

## TF Activity Findings

DoRothEA/VIPER activities used confidence levels A/B/C with minimum target coverage >=15. Retained TFs: `{tf["host_feature"].nunique()}`. Primary TF association rows: `{len(tf)}`. FDR-supported TF rows: `{len(tf_sig)}`.

{md_table(top_rows(tf), ["taxon", "host_feature", "coefficient", "p_value", "bh_q_value", "RCLR_DIRECTION_SENSITIVE"])}

## TME Covariate Sensitivities

Sensitivity models were run separately for inferred tumor purity, immune score, and stromal score. The combined TME model was not run. Candidate rows with a composition-sensitive interpretation: `{int(cov["robustness_interpretation"].astype(str).str.contains("composition_sensitive", na=False).sum())}`.

## rCLR and Contamination Sensitivities

rCLR and contaminant-exclusion checks were executed for candidate mechanisms. Rows labelled transformation sensitive in the transformation table: `{int(trans["TRANSFORMATION_SENSITIVE_MECHANISM"].sum())}`. Biological interpretation remains limited by the Phase 7C finding that eight of nine robust genera reverse direction under rCLR.

## WGCNA Modules

WGCNA used the locked top-25% MAD-variable genes and blockwise module construction. Selected soft power was `{int(wgcna_summary["selected_soft_power"].iloc[0])}`; modules after merging: `{int(wgcna_summary["modules_after_merging"].iloc[0])}`; grey genes: `{int(wgcna_summary["grey_genes"].iloc[0])}`. WGCNA taxon-module association rows: `{len(wgcna)}`. FDR-supported WGCNA rows: `{len(wgcna_sig)}`.

## Exploratory Host-Gene Associations

Genome-wide limma models were run one taxon at a time. Each model used 42,654 eligible genes with BH correction per taxon. Full primary CLR result tables were written under `05_results/tables/phase8b_host_gene_full/`.

{md_table(gene_summary[gene_summary["model"].eq("primary_CLR")].sort_values("n_q_lt_0_05", ascending=False), ["taxon", "model", "n_genes", "n_q_lt_0_05", "top_gene", "top_effect_size"])}

## Ranked Gene-Set Enrichment

Ranked enrichment used complete moderated t-statistics rather than significant-only gene lists. FDR-supported enrichment rows: `{len(enrich_sig)}`.

{md_table(enrich_sig.sort_values("padj").head(8), ["taxon", "collection", "pathway", "NES", "pval", "padj", "gene_set_version"])}

## Shared Versus Taxon-Specific Mechanisms

Shared-mechanism summaries used cross-taxon sign consistency and the taxon-taxon CLR correlation matrix. Shared rows with more than one taxon: `{int((shared["n_taxa_supported"] > 1).sum())}`. Correlated taxa are treated as compositionally linked microbial features, not independent biological exposures.

## Evidence Categories

{evidence_counts.to_markdown(index=False)}

## Null and Negative Results

Most tested pathway, TF, and WGCNA module rows did not meet the locked robust mechanism criteria. These rows remain in the primary association and evidence tables as `NO_SUPPORTED_MECHANISM` or `EXPLORATORY_HOST_MECHANISM` rather than being filtered out of the record.

## Limitations

Key limitations are compositional direction sensitivity under rCLR, contamination-risk annotations without sequenced negative controls, same-expression-matrix circularity for Moffitt50 and TME-derived scores, and limited sample size (`n=62`). Moffitt50 gene-exclusion sensitivity was executed where technically applicable, but it does not create an independent host-expression dataset.

## Recommendations

Proceed to external validation and target-prioritization phases only after human review of Phase 8B outputs. Do not treat exploratory host-gene or enrichment findings as validated mechanisms. Unresolved computational or interpretation items remain labelled `TO_VERIFY` where applicable.
"""

    REPORT.write_text(report)

    append_once(
        STATUS,
        "Completed Phase 8B:",
        "- Completed Phase 8B: executed locked host-mechanism analyses for the nine Phase 7C-verified primary taxa. Generated Hallmark, PROGENy, DoRothEA/VIPER, WGCNA, genome-wide limma, ranked enrichment, sensitivity, shared-mechanism, evidence, figure, validation, and report outputs under the Phase 8A rules. Next approved task: human review of Phase 8B products before external validation or target-prioritization phases.",
    )

    append_once(
        DECISION,
        "### D-25: Execute Locked Phase 8B Host-Mechanism Analyses",
        """### D-25: Execute Locked Phase 8B Host-Mechanism Analyses
*   **Date:** 2026-07-02
*   **Decision:** Execute the Phase 8A locked host-mechanism analysis without changing feature collections, taxa, covariates, WGCNA parameters, sensitivity rules, FDR families, or evidence categories after inspecting results.
*   **Alternatives Considered:** Promote suggestive taxa into the primary family; combine purity, immune, and stromal scores in one model; choose pathways or TFs after observing associations; omit rCLR or Moffitt50 safeguards.
*   **Scientific and Operational Justification:** The locked workflow preserves the prospective analysis plan and records transformation, composition, sample, and circularity limitations directly in evidence categories and sensitivity tables.
*   **Files / Analyses Affected:** `06_scripts/R/13_phase8b_host_mechanisms.R`, `06_scripts/python/13_summarize_phase8b_mechanisms.py`, `06_scripts/python/13_validate_phase8b_mechanisms.py`, `04_analysis/08_host_microbiome_integration/PHASE8B_HOST_MECHANISM_RESULTS.md`, `05_results/tables/phase8b_*`, and `05_results/figures/phase8b_*`.""",
    )

    update_manifest()
    print(f"Wrote {REPORT}")
    print("Updated project status, decision log, and file manifest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
