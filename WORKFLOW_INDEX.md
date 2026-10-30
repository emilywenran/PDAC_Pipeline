# AI-assisted PDAC Research — Workflow Index

> Scope: index of recovered human-to-agent instructions from ChatGPT history, **not** a verified record of successful execution. Entries are listed in the order in which they appear in the recovered archive.

## Provenance and interpretation

Publication copy updated on 2026-10-30: the source filename has been aligned with the supplied `PHASE_PROMPTS.txt`. Line ranges and links use GitHub line numbering based on LF characters; the supplied index counted 72 additional Unicode line separators. Prompt bytes, archive order and recorded evidence statuses are unchanged.

- Source: [`PHASE_PROMPTS.txt`](PHASE_PROMPTS.txt), recovered retrospectively from ChatGPT conversations and converted from an RTF-formatted file to UTF-8 plain text. The conversion retains text, but visual formatting may differ.
- Repository paths and output files mentioned inside prompts describe **requested** actions. Their creation, completion, validation, timing, and final acceptance must be verified independently using scripts, review reports, Git history and output manifests.
- `Antigravity` / `Codex` labels are preserved from the archive; the use of an independent AI reviewer does **not** by itself establish independent human review.
- Source line ranges below refer to the UTF-8 [`PHASE_PROMPTS.txt`](PHASE_PROMPTS.txt) file as archived; they permit direct lookup but are not commit identifiers.
- SHA-256 of the converted UTF-8 source file: `555d6f417989233f3b7064612d4c3efaddfb486630af3b9522b59e4428cf9dd8`.

## Recovered prompt sequence

| ID | Environment | Task (summary) | Instruction category | Source lines | Evidence status |
|---|---|---|---|---:|---|
| 1A | Antigravity | 1A: authoritative public-data source audit and accession mapping. | Data / provenance audit | [4–391](PHASE_PROMPTS.txt#L4-L391) | Prompt recovered; execution not independently verified |
| 1B | Codex | 1B: finalize the patient-level sample, clinical, and file manifests. | Research task / setup | [392–460](PHASE_PROMPTS.txt#L392-L460) | Prompt recovered; execution not independently verified |
| 2A | Codex | 2A: acquire and audit the processed GSE172356 host gene-expression matrix. | Data / provenance audit | [461–572](PHASE_PROMPTS.txt#L461-L572) | Prompt recovered; execution not independently verified |
| 2B | Codex | 2B: resolve missing expression values, define expression filtering and transformation, and generate an analysis-ready matrix for PDAC sub… | Research task / setup | [573–714](PHASE_PROMPTS.txt#L573-L714) | Prompt recovered; execution not independently verified |
| 3A | Antigravity | 3A: identify, verify, and lock the PDAC subtype reproduction methods before applying them to GSE172356. | Method design / lock | [715–852](PHASE_PROMPTS.txt#L715-L852) | Prompt recovered; execution not independently verified |
| 3B | Codex | 3B: reproduce PDAC molecular subtype assignments using the methods locked in Phase 3A, and quantify agreement with the verified public su… | Method design / lock | [853–1051](PHASE_PROMPTS.txt#L853-L1051) | Prompt recovered; execution not independently verified |
| 4A | Antigravity | 4A: prespecify and lock the statistical framework for evaluating whether basal-like, classical, and hybrid PDAC labels represent stable d… | Method design / lock | [1052–1322](PHASE_PROMPTS.txt#L1052-L1322) | Prompt recovered; execution not independently verified |
| 4B | Codex | 4B: execute the locked PDAC subtype stability analyses. | Analysis execution | [1323–1677](PHASE_PROMPTS.txt#L1323-L1677) | Prompt recovered; execution not independently verified |
| 5A | Antigravity | 5A: prespecify and lock the continuous basal–classical transcriptional-axis analysis. | Method design / lock | [1678–1812](PHASE_PROMPTS.txt#L1678-L1812) | Prompt recovered; execution not independently verified |
| 5B | Codex | 5B: execute the locked continuous basal–classical transcriptional-axis analyses. ## Mandatory K-Dense skill Use the installed K-Dense `st… | Analysis execution | [1813–2178](PHASE_PROMPTS.txt#L1813-L2178) | Prompt recovered; execution not independently verified |
| 6A | Codex | 6A: identify, extract, and audit the processed tumor microbiome abundance data associated with PRJNA719915. | Data / provenance audit | [2179–2312](PHASE_PROMPTS.txt#L2179-L2312) | Prompt recovered; execution not independently verified |
| 6B | Antigravity | 6B: prespecify and lock the tumor microbiome preprocessing, compositional transformation, and contamination-sensitivity framework. | Method design / lock | [2313–2546](PHASE_PROMPTS.txt#L2313-L2546) | Prompt recovered; execution not independently verified |
| 6C | Codex | 6C: execute the locked tumor microbiome preprocessing and generate analysis-ready compositional matrices. ## Mandatory K-Dense skills Use: | Analysis execution | [2547–2731](PHASE_PROMPTS.txt#L2547-L2731) | Prompt recovered; execution not independently verified |
| 7A | Antigravity | 7A: prespecify and lock the statistical framework for associations between the tumor microbiome and PDAC transcriptional states. | Method design / lock | [2732–3003](PHASE_PROMPTS.txt#L2732-L3003) | Prompt recovered; execution not independently verified |
| 7A.5 | Codex | 7A.5: calculate and validate host-derived tumor purity, immune, and stromal covariates before executing the locked microbiome-association… | Research task / setup | [3004–3267](PHASE_PROMPTS.txt#L3004-L3267) | Prompt recovered; execution not independently verified |
| 7B | Codex | 7B: execute the locked associations between the tumor microbiome and PDAC transcriptional states. | Analysis execution | [3268–3832](PHASE_PROMPTS.txt#L3268-L3832) | Prompt recovered; execution not independently verified |
| 7C | Antigravity | 7C: independent statistical and implementation review of the completed Phase 7B microbiome–host transcriptional-state association analysi… | AI-assisted independent review | [3833–4155](PHASE_PROMPTS.txt#L3833-L4155) | Prompt recovered; execution not independently verified |
| 8A | Antigravity | 8A: prespecify and lock the host-mechanism analysis framework for the verified tumor-microbiome associations. | Method design / lock | [4156–4578](PHASE_PROMPTS.txt#L4156-L4578) | Prompt recovered; execution not independently verified |
| 8A.5 | Codex | 8A.5: prepare and validate the R environment required for the locked Phase 8B host-mechanism analysis. | Analysis execution | [4579–4701](PHASE_PROMPTS.txt#L4579-L4701) | Prompt recovered; execution not independently verified |
| 8B | Codex | 8B: execute the locked host-mechanism analyses for the verified tumor-microbiome associations. ## Mandatory K-Dense skill Use the install… | Analysis execution | [4702–5133](PHASE_PROMPTS.txt#L4702-L5133) | Prompt recovered; execution not independently verified |
| 8C | Antigravity | 8C: independent statistical, implementation, and evidence review of the completed Phase 8B host–microbiome mechanism analysis. This is a … | AI-assisted independent review | [5134–5384](PHASE_PROMPTS.txt#L5134-L5384) | Prompt recovered; execution not independently verified |
| 9A | Antigravity | 9A: identify, evaluate, and lock the external-validation framework for the PDAC transcriptional-state and host–microbiome findings. This … | Method design / lock | [5385–5977](PHASE_PROMPTS.txt#L5385-L5977) | Prompt recovered; execution not independently verified |
| 9B1 | Codex | 9B1: execute the locked independent bulk-transcriptome validation analyses. | Analysis execution | [5978–6167](PHASE_PROMPTS.txt#L5978-L6167) | Prompt recovered; execution not independently verified |
| 9B1C | Antigravity | 9B1C: independent statistical and implementation review of the completed bulk-transcriptome external validation. This is a review phase. | AI-assisted independent review | [6168–6392](PHASE_PROMPTS.txt#L6168-L6392) | Prompt recovered; execution not independently verified |
| 9B1R | Codex | 9B1R: correct and rerun the independent bulk-transcriptome external validation following the Phase 9B1C audit. This is a targeted reanaly… | Correction / reanalysis | [6393–6703](PHASE_PROMPTS.txt#L6393-L6703) | Prompt recovered; execution not independently verified |
| 9B1C2 | Antigravity | 9B1C2: independent review of the corrected Phase 9B1R bulk-transcriptome external-validation analysis. This is a verification phase. | AI-assisted independent review | [6704–7001](PHASE_PROMPTS.txt#L6704-L7001) | Prompt recovered; execution not independently verified |
| 9B2 | Codex | 9B2: execute patient-aware single-cell cellular-source evaluation of the partially replicated and discovery-supported PDAC host programs.… | Analysis execution | [7002–7474](PHASE_PROMPTS.txt#L7002-L7474) | Prompt recovered; execution not independently verified |
| 9B1R-minor | Codex | the Phase 9B1R minor correction required by the Phase 9B1C2 review. | Correction / reanalysis | [7475–7571](PHASE_PROMPTS.txt#L7475-L7571) | Prompt recovered; execution not independently verified |
| 9B1C2-closure | Antigravity | a limited closure review for the remaining Phase 9B1C2 minor correction concerning FIND_05 / FIND_07 TF evidence classification. | AI-assisted independent review | [7572–7618](PHASE_PROMPTS.txt#L7572-L7618) | Prompt recovered; execution not independently verified |
| 9A.1 | Antigravity | 9A.1: reconcile the locked PRIORITY_1 single-cell cohort set before Phase 9B2 execution. This is a metadata and planning correction only. | Correction / reanalysis | [7619–7780](PHASE_PROMPTS.txt#L7619-L7780) | Prompt recovered; execution not independently verified |
| 9A.2 | Antigravity | 9A.2: correct the single-cell dataset identities, accessions, provenance, and Phase 9B2 cohort definitions. This is a metadata and planni… | Correction / reanalysis | [7781–8672](PHASE_PROMPTS.txt#L7781-L8672) | Prompt recovered; execution not independently verified |
| 9A.3 | Antigravity | a limited Phase 9A.3 execution-scope correction for Phase 9B2. The Phase 9B2 restart stopped correctly before data acquisition because th… | Correction / reanalysis | [8673–8848](PHASE_PROMPTS.txt#L8673-L8848) | Prompt recovered; execution not independently verified |
| 9B2C | Antigravity | 9B2C: independent statistical, implementation, annotation, and evidence review of the completed Phase 9B2 primary single-cell cellular-so… | AI-assisted independent review | [8849–9223](PHASE_PROMPTS.txt#L8849-L9223) | Prompt recovered; execution not independently verified |
| 9B2R | Codex | 9B2R: correct and rerun the primary single-cell cellular-source analysis following the Phase 9B2C independent review. This is a consolida… | Correction / reanalysis | [9224–9621](PHASE_PROMPTS.txt#L9224-L9621) | Prompt recovered; execution not independently verified |
| 9B2C2 | Antigravity | 9B2C2: a complete independent statistical, implementation, annotation, biological-replicate, negative-control, and evidence review of the… | AI-assisted independent review | [9622–10113](PHASE_PROMPTS.txt#L9622-L10113) | Prompt recovered; execution not independently verified |
| 9B3A | Antigravity | 9B3A: spatial-transcriptomic validation planning, dataset qualification, feature hierarchy locking, and prospective statistical-design lo… | Method design / lock | [10114–10604](PHASE_PROMPTS.txt#L10114-L10604) | Prompt recovered; execution not independently verified |
| 9B3A.1 | Antigravity | 9B3A.1: limited spatial-design consistency correction before Phase 9B3 execution. This is a planning, method-lock, metadata, and validato… | Correction / reanalysis | [10605–10829](PHASE_PROMPTS.txt#L10605-L10829) | Prompt recovered; execution not independently verified |
| 9B3A.2 | Antigravity | 9B3A.2: final ROI-pairing and low-patient-count spatial-design correction before Phase 9B3 execution. This is a narrow planning, model-lo… | Correction / reanalysis | [10830–11092](PHASE_PROMPTS.txt#L10830-L11092) | Prompt recovered; execution not independently verified |
| 9B3B | Codex | 9B3B: execute the prospectively locked spatial-transcriptomic validation. | Analysis execution | [11093–11314](PHASE_PROMPTS.txt#L11093-L11314) | Prompt recovered; execution not independently verified |
| 9B3C | Antigravity | 9B3C: complete independent review of the executed Phase 9B3B spatial validation. | AI-assisted independent review | [11315–11561](PHASE_PROMPTS.txt#L11315-L11561) | Prompt recovered; execution not independently verified |
| 9B3R0 | Antigravity | 9B3R0: dual-model pre-reanalysis audit of the failed Phase 9B3B spatial validation. | Correction / reanalysis | [11562–11653](PHASE_PROMPTS.txt#L11562-L11653) | Prompt recovered; execution not independently verified |
| 9B3R | Codex | 9B3R: corrected reanalysis of the failed Phase 9B3B spatial validation. | Correction / reanalysis | [11654–11800](PHASE_PROMPTS.txt#L11654-L11800) | Prompt recovered; execution not independently verified |
| 9B3C2 | Antigravity | 9B3C2: complete independent review of the corrected Phase 9B3R spatial validation. | AI-assisted independent review | [11801–11928](PHASE_PROMPTS.txt#L11801-L11928) | Prompt recovered; execution not independently verified |
| 10A | Antigravity | 10A: prospective cross-layer evidence-synthesis and target- prioritization planning. This is a planning and method-lock phase only. | Method design / lock | [11929–12018](PHASE_PROMPTS.txt#L11929-L12018) | Prompt recovered; execution not independently verified |
| 10B-R | Codex | 10B-R: corrected cross-layer evidence synthesis and target prioritization reanalysis after Phase 10C failed the first Phase 10B attempt. | Correction / reanalysis | [12019–12113](PHASE_PROMPTS.txt#L12019-L12113) | Prompt recovered; execution not independently verified |
| 10C2 | Antigravity | 10C2: complete independent review of the corrected Phase 10B-R cross-layer evidence synthesis and target-prioritization reanalysis. This … | AI-assisted independent review | [12114–12231](PHASE_PROMPTS.txt#L12114-L12231) | Prompt recovered; execution not independently verified |
| 11A | Antigravity | 11A: manuscript structure, claim-control, and figure/table planning. This is a manuscript-planning phase only. | Method design / lock | [12232–12299](PHASE_PROMPTS.txt#L12232-L12299) | Prompt recovered; execution not independently verified |

## Major study stages

| Stages | Role in the study |
|---|---|
| 1A–2B | Public dataset accession checks, patient-level mapping, RNA-expression QC, missingness handling, and creation of an analysis-ready expression matrix |
| 3A–5B | Subtype method verification/reproduction, cluster stability, and continuous basal–classical axis construction |
| 6A–7C | Tumour microbiome provenance and preprocessing, host–microbiome association testing and audit |
| 8A–8C | Host pathway, regulatory and mechanistic programme analyses and audit |
| 9A–9B3C2 | External bulk, single-cell and spatial study planning, execution, multi-stage review, correction and reanalysis |
| 10A–10C2 | Cross-layer evidence synthesis and exploratory target prioritisation, with documented correction/review prompts |
| 11A | Manuscript planning and claim-control instructions (the archive extends beyond Phase 10) |

## Important correction loops and gaps

- **Bulk external validation:** `9B1 → 9B1C → 9B1R → 9B1C2`, followed by recovered minor correction and closure instructions. Verify the resulting PASS/FAIL outcome from saved review reports rather than treating the requested final decision as a result.
- **Single-cell provenance/execution:** An initial `9B2` execution instruction appears before subsequent `9A.1`, `9A.2` and `9A.3` reconciliation/correction prompts; this reflects recovered archive order, not the ideal phase order. Confirm stopped/restarted runs in the repository.
- **Single-cell review:** `9B2C → 9B2R → 9B2C2`. Check acceptance and which revised analyses entered the manuscript.
- **Spatial:** `9B3A → 9B3A.1 → 9B3A.2 → 9B3B → 9B3C → 9B3R0 → 9B3R → 9B3C2`. Preserve original and corrected artifacts, including adverse review findings.
- **Evidence synthesis:** `10A → 10B-R → 10C2` is recovered, but this archive does not include a clearly identifiable verbatim initial `10B` instruction or initial `10C` instruction. Their existence is mentioned indirectly; do not reconstruct their verbatim text.
- Separate research instructions may remain missing, including later manuscript drafting and the exploratory external RNA-seq-derived microbial-profile calculation. Cross-check against the final manuscript and repository.

## Suggested next audit (minimal manual effort)

Add references, where available, to the execution report and Git commit for the critical correction loops rather than annotating every recovered prompt. At minimum verify the source-to-result chain for 7C, 8C, 9B1C2, 9B2C2, 9B3C2 and 10C2. Record `NOT_RECOVERED` where original prompts cannot be found.

## Public disclosure notes

- Remove or deliberately retain local home-directory paths according to your privacy policy; if you redact them, mark the redactions explicitly and keep an untouched internal original.
- Review for secrets and confidential information manually. A preliminary scan found local filesystem paths, but no obvious conventional API-key patterns.
- If prompts are edited for publication, describe them as *redacted recovered prompts*, not verbatim prompts. Do not silently tidy wording or reorder history.
