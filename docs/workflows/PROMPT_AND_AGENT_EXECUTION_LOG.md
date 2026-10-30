# Prompt and Executing-Agent Record

## 📌 Current status on 2026-10-30

Recovered phase prompts are now available in [PHASE_PROMPTS.txt](../../PHASE_PROMPTS.txt), with navigation and source-line links in [WORKFLOW_INDEX.md](../../WORKFLOW_INDEX.md). The archive contains 47 indexed instructions from Phase 1A through Phase 11A, including correction and closure instructions. This supersedes the initial packaging audit's finding that no verbatim prompt files were available in the workspace at that time.

The archive was supplied retrospectively from ChatGPT conversations and converted from RTF into UTF-8 plain text, as described in the supplied index. The public prompt file is byte-identical to that supplied text. Its SHA-256 is `555d6f417989233f3b7064612d4c3efaddfb486630af3b9522b59e4428cf9dd8`.

The archive preserves requested instructions, including expected outputs and requested review decisions. Those instructions are not evidence that the requested work completed or passed review. The source-line ranges are navigation aids, not execution timestamps. Local development paths are retained, and no missing prompt text has been reconstructed. Exact model versions and human review must be supported separately by project records.

## 📚 How the records relate

| Record | Interpretation |
| --- | --- |
| Recovered prompts and workflow index | Original instructions as supplied after retrospective recovery; incomplete coverage is documented in the index |
| Skill-use log | Agent-authored metadata about recorded skill invocation, including unavailable skills |
| Scripts, outputs and review reports | Saved implementation and execution/review evidence; verify each claim against the relevant artifact |
| Original workspace commit table below | Historical analysis-workspace chronology assembled during the initial packaging audit; not a transcript or the commit history of this curated GitHub repository |

The following table is retained as the earlier execution-metadata summary. Its commit hashes refer to the originating analysis workspace, whose Git history is distinct from the curated public-release repository. A recovered role label does not independently confirm the executing model or final outcome.

## Phase-Level Execution Record (not original prompts — agent-authored summaries only)

| # | Date | Commit | Phase | Executing agent (per SKILL_USAGE_LOG.tsv) | Recorded task description (commit message, supplemented by SKILL_USAGE_LOG notes) |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | 2026-06-30 | `1f69ecb` | Phase 0 | Not documented in project records | Initialize PDAC subtype stability project |
| 2 | 2026-06-30 | `1cb3d24` | Phase 0 (planning) | Not documented in project records | Add scientific planning and analysis guardrails (added PROJECT_CHARTER / ANALYSIS_PLAN / HYPOTHESES / RISK_REGISTER, etc.) |
| 3 | 2026-06-30 | `24ee042` | Phase 0 (correction, R-LOG-01) | Antigravity (per the author field in DECISION_LOG R-LOG-01) | Correct shotgun metagenomics analysis plan (revised from a 16S design to a shotgun metagenomics workflow) |
| 4 | 2026-06-30 | `8b31252` | Phase 1A | Antigravity | Finalized Phase 1A (accession audit, `database-lookup` skill) |
| 5 | 2026-06-30 | `0e5cdf6` | Phase 1B | Not documented in project records | Finalize PDAC patient-level sample mapping |
| 6 | 2026-06-30 | `f2cdd6a` | Phase 2A | Codex | Audit GSE172356 processed expression matrix (`exploratory-data-analysis` skill) |
| 7 | 2026-06-30 | `58e355a` | Phase 2B | Codex | Prepare analysis-ready PDAC expression matrix |
| 8 | 2026-07-01 | `33a71cf` | Phase 3A | Antigravity | Lock PDAC subtype reproduction methods (`database-lookup`/`citation-management`/`experimental-design`) |
| 9 | 2026-07-01 | `13f1bef` | Phase 3B | Not documented in project records | Reproduce PDAC molecular subtype assignments |
| 10 | 2026-07-01 | `7d832e8` | Phase 4A | Antigravity | Lock PDAC subtype stability analysis plan |
| 11 | 2026-07-01 | `f32596e` | Phase 4B | Codex | Evaluate PDAC subtype clustering stability |
| 12 | 2026-07-01 | `17dbf9d` | Phase 5A | Antigravity | Lock PDAC continuous subtype axis analysis plan |
| 13 | 2026-07-01 | `8cccd29` | Phase 5A correction (D-14) | Antigravity | Reconcile Moffitt continuous axis gene sets (LEMD1 gene-set correction) |
| 14 | 2026-07-01 | `b51aee0` | Phase 5B | Codex | Evaluate PDAC continuous transcriptional axis |
| 15 | 2026-07-01 | `29bfd27` | Phase 6A | Codex | Audit processed PDAC tumor microbiome data |
| 16 | 2026-07-01 | `003fa58` | Phase 6B | Antigravity | Lock PDAC microbiome preprocessing framework |
| 17 | 2026-07-01 | `97f2e78` | Phase 6B correction | Antigravity | Amend PDAC microbiome preprocessing framework |
| 18 | 2026-07-01 | `cc3d735` | Phase 6C | Codex | Prepare analysis-ready PDAC microbiome matrices |
| 19 | 2026-07-01 | `57d548b` | Phase 7A | Antigravity | Lock PDAC microbiome host-state association models |
| 20 | 2026-07-01 | `4745442` | Phase 7A.5 | Codex | Add host tumor microenvironment covariates (ESTIMATE covariates) |
| 21 | 2026-07-02 | `7f9d1de` | Phase 7B | Codex | Analyze microbiome associations with PDAC transcriptional state |
| 22 | 2026-07-02 | `26e4d05` | Phase 7C | Antigravity | Review PDAC microbiome association results (independent review) |
| 23 | 2026-07-02 | `b23dc27` | Phase 8A | Antigravity | Lock PDAC host microbiome mechanism analysis |
| 24 | 2026-07-02 | `e4346d1` | Phase 8A.5 | Not documented in project records | Prepare Phase 8 host mechanism R environment |
| 25 | 2026-07-02 | `b7f6b11` | Phase 8A.5 (continued) | Not documented in project records | Prepare Phase 8 host mechanism R environment |
| 26 | 2026-07-02 | `63bb872` | Phase 8B | Codex | Checkpoint Phase 8B host microbiome mechanism analysis |
| 27 | 2026-07-02 | `1d2828c` | Phase 8C | Antigravity | Review PDAC host microbiome mechanism results (independent review) |
| 28 | 2026-07-03 | `0b4c5e0` | Phase 9A | Antigravity | Lock PDAC external validation framework |
| 29 | 2026-07-03 | `ee3db99` | Phase 9B1 → 9B1R | Codex (9B1 execution) / Antigravity (9B1C review, FAIL) / Codex (9B1R fix) | Correct PDAC bulk external validation analysis (the original Phase 9B1 execution and the 9B1C rejection do not appear to have their own separate commits, and may not have been committed individually) |
| 30 | 2026-07-03 | `72fa7cc` | Phase 9B1C2 | Antigravity | Review corrected PDAC bulk external validation |
| 31 | 2026-07-03 | `069c32d` | Phase 9B1C2 (minor correction) | Antigravity | Correct TF evidence classification in bulk validation |
| 32 | 2026-07-03 | `3c8bd92` | Phase 9B1C2 (closure) | Antigravity | Close corrected bulk validation review |
| 33 | 2026-07-03 | `ddb7b9f` | Phase 9B1R (output update) | Codex | Update corrected bulk validation outputs |
| 34 | 2026-07-03 | `753448a` | Phase 9 (tooling) | Not documented in project records | Add Phase 9 workflow utility scripts |
| 35 | 2026-07-03 | `c0b87cc` | Engineering maintenance | Not documented in project records | Ignore external processed data and local caches (`.gitignore` update) |
| 36 | 2026-07-03 | `221cbbd` | Phase 9A.1/9A.2 | Antigravity | Correct Phase 9 single-cell dataset provenance |
| 37 | 2026-07-03 | `359a0b2` | Phase 9A.2 | Antigravity | Add Phase 9A.2 provenance correction utilities |
| 38 | 2026-07-03 | `7b053d3` | Phase 9A.3 | Antigravity | Define Phase 9B2 primary execution scope |
| 39 | 2026-07-03 | `c8f8f27` | Phase 9B2R/9B2C2 | Codex (execution) / Antigravity (review) | Finalize corrected single-cell cellular-source validation |
| 40 | 2026-07-03 | `fd3bf4a` | Phase 9B2/9B2C (archived) | Codex (execution) / Antigravity (review, FAIL) | Archive initial Phase 9B2 analysis and independent review |
| 41 | 2026-07-03 | `2fbfbb6` | Phase 9B3A/9B3A.1/9B3A.2 | Antigravity | Lock PDAC spatial validation design and cohort hierarchy |
| 42 | 2026-07-04 | `ff919fd` | Phase 9B3R + Phase 10 (combined commit) | Codex (9B3R execution) / Antigravity (9B3C2 review, Phase 10A/10B/10C2) | Finalize corrected PDAC spatial validation and target prioritization |
| 43 | 2026-07-04 | `49c473f` | Phase 11A | Antigravity | Lock PDAC manuscript claims and figure plan |
| 44 | 2026-07-04 | `5d90d49` | Phase 11B | Codex | Finalize Phase 11B manuscript draft validation |
| 45 | 2026-07-04 | `a9a1ef6` | Phase 11C | Antigravity | Complete Phase 11C manuscript independent review |
| 46 | 2026-07-04 | `81322e9` | Phase 11D | Antigravity | Assemble Phase 11D full manuscript draft |
| 47 | 2026-07-04 | `4cf5a2e` | Phase 11E | Antigravity | Complete Phase 11E language and format review |
| 48 | 2026-07-04 | `188e4b8` | Phase 11F | Antigravity | Complete Phase 11F final claim audit |
| 49 | 2026-07-04 | `1ff8de2` | Phase 11G(-R1/R2/R3) | Antigravity | Complete Phase 11G reference and callout repair |
| 50 | 2026-07-04 | `fc2c992` | Post-Phase 11G | Antigravity | Finalize post-Phase 11G workspace cleanup |
| 51 | 2026-07-04 | `b16788d` | Phase 11H | Antigravity | Assemble Phase 11H submission package |
| 52 | 2026-07-04 | `106ac5a` | Phase 11I-A | Antigravity | Complete Phase 11I-A final QA and journal gap audit |

**Note**: `Antigravity` / `Codex` are the role names recorded in `SKILL_USAGE_LOG.tsv`. Per the existing document `09_docs/workflows/PDAC_PROJECT_WORKFLOW_AGENT_SKILL_SUMMARY.md`, they are said to correspond to Gemini 3.1 Pro High and Gemini 3.5 Flash/Flash High respectively — but that mapping itself **is not independently verifiable from any file in this repository**; it is only recorded in the existing summary document, so cite it with appropriate caution.

---

## 📝 Remaining record gaps

The supplied index identifies missing or incompletely recovered instructions, including initial Phase 10B/10C and later manuscript or supplementary-analysis work. Do not fill those gaps with reconstructed text described as verbatim. For future work, save the instruction when issued and link it to the resulting script, output and review record. Retrospective recovery of the present archive does not establish that every earlier interaction was recorded contemporaneously.
