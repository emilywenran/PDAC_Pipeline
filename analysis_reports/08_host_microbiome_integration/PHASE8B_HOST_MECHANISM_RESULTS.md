# Phase 8B Host Mechanism Results

## Execution Status

Phase 8B executed the Phase 8A locked host-mechanism analyses for the nine Phase 7C-verified primary genera. Runtime validation passed all hard-stop checks: 62 aligned patients, nine primary taxa present in primary CLR and rCLR matrices, Phase 7B directions loaded, expression and microbiome order aligned by `patient_id`, no duplicate primary identifiers, no missing or infinite primary values, and required Phase 8 packages loaded from the project-local `renv` library.

`renv` note: R startup required `RENV_CONFIG_SANDBOX_ENABLED=FALSE` because the renv sandbox lock path blocked startup in this workspace. The active library was still the project-local renv library recorded in `phase8b_runtime_validation.tsv`.

## Pathway Activity Findings

Hallmark ssGSEA scores were generated from MSigDB `2026.1.Hs`. PROGENy scores used package version `1.32.0` with the locked top-100 model. Primary pathway association rows: `576`. FDR-supported pathway rows: `50`.

| taxon         | host_feature_collection   | host_feature                 |   coefficient |     p_value |   bh_q_value | RCLR_DIRECTION_SENSITIVE   |
|:--------------|:--------------------------|:-----------------------------|--------------:|------------:|-------------:|:---------------------------|
| Ochrobactrum  | MSigDB_Hallmark           | HALLMARK_SPERMATOGENESIS     |    -0.0126216 | 1.09672e-06 |  5.48358e-05 | False                      |
| Azoarcus      | MSigDB_Hallmark           | HALLMARK_SPERMATOGENESIS     |     0.0115646 | 2.70314e-05 |  0.00135157  | True                       |
| Cutibacterium | MSigDB_Hallmark           | HALLMARK_SPERMATOGENESIS     |     0.0108036 | 2.77437e-05 |  0.00138718  | True                       |
| Azoarcus      | MSigDB_Hallmark           | HALLMARK_PANCREAS_BETA_CELLS |     0.0738242 | 0.000484525 |  0.00801675  | True                       |
| Azoarcus      | MSigDB_Hallmark           | HALLMARK_ALLOGRAFT_REJECTION |     0.0325606 | 0.00051183  |  0.00801675  | True                       |
| Azoarcus      | MSigDB_Hallmark           | HALLMARK_IL2_STAT5_SIGNALING |     0.0155441 | 0.00064134  |  0.00801675  | True                       |
| Burkholderia  | PROGENy                   | p53                          |    -0.351127  | 0.000639049 |  0.00894669  | True                       |
| Candida       | PROGENy                   | PI3K                         |    -0.381627  | 0.000709545 |  0.00993363  | True                       |

## TF Activity Findings

DoRothEA/VIPER activities used confidence levels A/B/C with minimum target coverage >=15. Retained TFs: `220`. Primary TF association rows: `1980`. FDR-supported TF rows: `247`.

| taxon        | host_feature   |   coefficient |     p_value |   bh_q_value | RCLR_DIRECTION_SENSITIVE   |
|:-------------|:---------------|--------------:|------------:|-------------:|:---------------------------|
| Ochrobactrum | MBD1           |     -0.308662 | 6.9952e-09  |  1.53894e-06 | False                      |
| Ochrobactrum | SNAPC4         |      0.21492  | 3.46101e-07 |  3.44879e-05 | False                      |
| Ochrobactrum | IRF3           |      0.306746 | 4.8966e-07  |  3.44879e-05 | False                      |
| Ochrobactrum | MXI1           |     -0.225462 | 6.27052e-07 |  3.44879e-05 | False                      |
| Ochrobactrum | MNT            |     -0.362559 | 1.56295e-06 |  6.87699e-05 | False                      |
| Ochrobactrum | SNAI2          |     -0.451449 | 2.97174e-06 |  0.000108964 | False                      |
| Candida      | ZNF24          |      0.190717 | 6.90143e-07 |  0.000151831 | True                       |
| Ochrobactrum | TEAD4          |      0.330905 | 1.27447e-05 |  0.000359699 | False                      |

## TME Covariate Sensitivities

Sensitivity models were run separately for inferred tumor purity, immune score, and stromal score. The combined TME model was not run. Candidate rows with a composition-sensitive interpretation: `1281`.

## rCLR and Contamination Sensitivities

rCLR and contaminant-exclusion checks were executed for candidate mechanisms. Rows labelled transformation sensitive in the transformation table: `649`. Biological interpretation remains limited by the Phase 7C finding that eight of nine robust genera reverse direction under rCLR.

## WGCNA Modules

WGCNA used the locked top-25% MAD-variable genes and blockwise module construction. Selected soft power was `5`; modules after merging: `16`; grey genes: `2886`. WGCNA taxon-module association rows: `144`. FDR-supported WGCNA rows: `39`.

## Exploratory Host-Gene Associations

Genome-wide limma models were run one taxon at a time. Each model used 42,654 eligible genes with BH correction per taxon. Full primary CLR result tables were written under `05_results/tables/phase8b_host_gene_full/`.

| taxon            | model       |   n_genes |   n_q_lt_0_05 | top_gene   |   top_effect_size |
|:-----------------|:------------|----------:|--------------:|:-----------|------------------:|
| Ochrobactrum     | primary_CLR |     42654 |          4845 | RF01956    |          0.574003 |
| Cutibacterium    | primary_CLR |     42654 |          3682 | RF01956    |         -0.533976 |
| Azoarcus         | primary_CLR |     42654 |          3373 | RNA5SP429  |         -0.682229 |
| Burkholderia     | primary_CLR |     42654 |          2978 | HR         |         -0.924587 |
| Candida          | primary_CLR |     42654 |          2631 | CHCHD2P11  |         -0.872914 |
| Rhizobium        | primary_CLR |     42654 |           305 | UPK3B      |         -0.983965 |
| Ensifer          | primary_CLR |     42654 |            81 | LINC00396  |          0.698749 |
| Chryseobacterium | primary_CLR |     42654 |             8 | SPINK6     |          1.25836  |
| Herbaspirillum   | primary_CLR |     42654 |             0 | AC010265.1 |          0.564922 |

## Ranked Gene-Set Enrichment

Ranked enrichment used complete moderated t-statistics rather than significant-only gene lists. FDR-supported enrichment rows: `3376`.

| taxon        | collection   | pathway                          |      NES |        pval |        padj | gene_set_version   |
|:-------------|:-------------|:---------------------------------|---------:|------------:|------------:|:-------------------|
| Burkholderia | Hallmark     | HALLMARK_MYC_TARGETS_V2          | -2.13696 | 0.000188608 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_DNA_REPAIR              | -2.2125  | 0.000189573 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_E2F_TARGETS             | -3.02727 | 0.000188324 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_G2M_CHECKPOINT          | -2.96306 | 0.000188964 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_GLYCOLYSIS              | -2.2897  | 0.000188324 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_IL2_STAT5_SIGNALING     |  1.72525 | 0.000213129 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_IL6_JAK_STAT3_SIGNALING |  2.25527 | 0.000213995 | 0.000534988 | 2026.1.Hs          |
| Burkholderia | Hallmark     | HALLMARK_INFLAMMATORY_RESPONSE   |  2.37286 | 0.000212269 | 0.000534988 | 2026.1.Hs          |

## Shared Versus Taxon-Specific Mechanisms

Shared-mechanism summaries used cross-taxon sign consistency and the taxon-taxon CLR correlation matrix. Shared rows with more than one taxon: `176`. Correlated taxa are treated as compositionally linked microbial features, not independent biological exposures.

## Evidence Categories

| category                           |    n |
|:-----------------------------------|-----:|
| NO_SUPPORTED_MECHANISM             | 1894 |
| EXPLORATORY_HOST_MECHANISM         |  470 |
| TRANSFORMATION_SENSITIVE_MECHANISM |  273 |
| ROBUST_HOST_MECHANISM              |   43 |
| COMPOSITION_SENSITIVE_MECHANISM    |   20 |

## Null and Negative Results

Most tested pathway, TF, and WGCNA module rows did not meet the locked robust mechanism criteria. These rows remain in the primary association and evidence tables as `NO_SUPPORTED_MECHANISM` or `EXPLORATORY_HOST_MECHANISM` rather than being filtered out of the record.

## Limitations

Key limitations are compositional direction sensitivity under rCLR, contamination-risk annotations without sequenced negative controls, same-expression-matrix circularity for Moffitt50 and TME-derived scores, and limited sample size (`n=62`). Moffitt50 gene-exclusion sensitivity was executed where technically applicable, but it does not create an independent host-expression dataset.

## Recommendations

Proceed to external validation and target-prioritization phases only after human review of Phase 8B outputs. Do not treat exploratory host-gene or enrichment findings as validated mechanisms. Unresolved computational or interpretation items remain labelled `TO_VERIFY` where applicable.
