# Phase 11G-R2: Reference Verification Audit

## 1. Overview
This audit verifies the internal consistency and mapped claims of the references cited in the Phase 11G-R1b manuscript (`04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md`) against `09_docs/references/phase11g_references.bib` and `05_results/tables/phase11gr1_reference_mapping.tsv`. 

## 2. Methodology
- Extracted references 1 through 14 from the mapping table and the `.bib` file.
- Fetched metadata from PubMed/CrossRef using the supplied DOIs to verify PMIDs, titles, and author lists.
- Checked whether each reference supports the mapped claims.
- Assessed internal consistency of PMIDs versus DOIs for all references.

## 3. Findings

### Total References Checked: 14

- **PASS**: 10 references
- **FAIL**: 4 references
- **UNVERIFIED**: 0 references

### Problematic References

1. **[3] Bear et al. (2020)**
   - **Issue**: PMID `32673551` in the mapping file does not match the DOI `10.1016/j.ccell.2020.08.004`. The provided PMID actually points to an unrelated study ("Telomere attrition with age in a wild amphibian population"). The true PMID for the Bear et al. paper is `32946773`.
   - **Claim Support**: SUPPORTED. The paper itself supports the claim, but the metadata is internally inconsistent.
   - **Required Action**: Update mapping table PMID to `32946773`.

2. **[6] Nejman et al. (2020)**
   - **Issue**: PMID `32471947` in the mapping file does not match the DOI `10.1126/science.aay9189`. The provided PMID belongs to an unrelated paper ("HTLV-1 induces T cell malignancy..."). The true PMID is `32467386`.
   - **Claim Support**: SUPPORTED.
   - **Required Action**: Update mapping table PMID to `32467386`.

3. **[7] Geller et al. (2017)**
   - **Issue**: PMID `28916614` in the mapping file does not match the DOI `10.1126/science.aah5043`. The provided PMID belongs to a paper on ALS. The true PMID is `28912244`.
   - **Claim Support**: SUPPORTED.
   - **Required Action**: Update mapping table PMID to `28912244`.

4. **[13] Boehm et al. (2015)**
   - **Issue**: This reference appears to be entirely fabricated or mismatched. The title "An armamentarium of therapeutic targets: state of the art" is not a real paper by Boehm. The DOI provided (`10.1038/nrg3915`) corresponds to a paper titled "Classifying pathogenic variation" by Bahcall. The PMID `25686524` corresponds to a paper on P2X7 receptor-mediated analgesia.
   - **Claim Support**: UNSUPPORTED.
   - **Required Action**: Remove citation [13] from the manuscript and mapping table, or replace it with a valid reference that supports target prioritisation in oncology.

### Manuscript Citations
All citation numbers `[1]` through `[14]` exist in the manuscript and map to entries in the `.bib` file.

### Required Manuscript Corrections
- The manuscript cites `[12,13]` for target prioritisation. Since `[13]` is a fabricated/mismatched reference, the manuscript text will need to be corrected in a future phase to either cite only `[12]` or replace `[13]` with a valid source. **No manuscript edits were performed during this verification phase**.

## 4. Final Decision
**MAJOR_REVISION_REQUIRED**. Four references contain inconsistent or fabricated metadata, including one completely mismatched reference (Boehm 2015) that cannot support the manuscript's claims. These metadata inconsistencies and the unsupported claim must be addressed in a subsequent phase.
