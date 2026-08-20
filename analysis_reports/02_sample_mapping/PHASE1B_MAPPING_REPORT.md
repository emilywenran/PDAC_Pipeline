# Phase 1B Mapping Report

## Final Mapping Counts

- Unique patients: 62
- RNA samples: 62 GEO samples / 62 RNA BioSamples / 62 RNA experiments / 62 RNA runs
- Microbiome samples: 62 microbiome BioSamples / 62 microbiome experiments / 62 microbiome runs
- Mapping basis: Phase 1A verified match between GEO tumor number and microbiome `source_material_id` / `potential_patient_id`.
- Mapping status: all 62 patient rows are `VERIFIED`.

## Uniqueness Checks

| Identifier class | Rows | Unique IDs | Duplicates |
| :--- | ---: | ---: | :--- |
| GEO sample | 62 | 62 | None |
| RNA BioSample | 62 | 62 | None |
| RNA experiment | 62 | 62 | None |
| RNA run | 62 | 62 | None |
| Microbiome BioSample | 62 | 62 | None |
| Microbiome experiment | 62 | 62 | None |
| Microbiome run | 62 | 62 | None |

## Subtype Counts

| Subtype | Patients |
| :--- | ---: |
| Basal | 17 |
| Classical | 22 |
| Hybrid | 23 |

Subtype source: Guo et al. 2021 Supplementary Data 4, worksheet `Figure1.SampleGroup`, Chan-Seng-Yue subtype labels. Phase 1A found these labels consistent with microbiome BioSample subtype descriptions under the tumor-number/source-material mapping; no subtype conflicts remain unresolved.

## Clinical Missingness

| Variable | Available | Missing | Notes |
| :--- | ---: | ---: | :--- |
| overall_survival_time | 53 | 9 | Days from Supplementary Data 4 `Figure3. Survival`. |
| overall_survival_event | 53 | 9 | Published status code retained exactly as supplied (`1`/`2`); no event recoding performed in Phase 1B. |
| age | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |
| sex | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |
| stage | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |
| grade | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |
| treatment | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |
| disease_free_survival_time | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |
| disease_free_survival_event | 0 | 62 | Not present in verified Phase 1A public/supplementary metadata. |

Patients missing overall survival variables: PDAC_001, PDAC_007, PDAC_012, PDAC_017, PDAC_032, PDAC_033, PDAC_042, PDAC_043, PDAC_044.

Published survival status counts among available rows: {'2': 37, '1': 16}.

## Unresolved Conflicts

None. No duplicated GEO, BioSample, experiment, or run IDs were detected in the patient crosswalk. Supplementary subtype labels are complete for all 62 mapped microbiome aliases.

## Eligibility

- Subtype analyses: all 62 patients are eligible.
- Survival analysis: 53 patients have both overall survival time and published survival status available. The status code is intentionally retained as published and must be interpreted before survival modeling.
- Phase 2 expression-data processing: may proceed for all 62 RNA samples. Phase 2 should not perform subtype reclassification or survival modeling until the downstream phase explicitly authorizes those tasks and confirms survival event coding.

## Outputs Created

- `01_metadata/sample_manifest.tsv`
- `01_metadata/clinical_metadata.tsv`
- `01_metadata/file_manifest.tsv`
- `01_metadata/rna_microbiome_patient_crosswalk.tsv`
