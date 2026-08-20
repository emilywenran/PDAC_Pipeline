# Phase 9B1R Correction Log

The original Phase 9B1 implementation is preserved as an audit artifact and is superseded by Phase 9B1R.

   finding_id severity
       <char>   <char>
1:    FIND_01 CRITICAL
2:    FIND_02    MAJOR
3:    FIND_03    MAJOR
4:    FIND_04    MAJOR
5:    FIND_05    MAJOR
6:    FIND_06 MODERATE
                                                                                                      finding
                                                                                                       <char>
1:                                  Omission of PurIST model intercept (-6.815) in logistic link calculation.
2:                                                   Violation of missing-gene policy for PurIST in GSE62452.
3:                               Violation of missing-gene policy (80% coverage threshold) for WGCNA modules.
4:                Use of 15-gene proxy set and rank-mean instead of ssGSEA via decoupleR on the full pathway.
5: Proxy single-gene expression used instead of VIPER TF activity, and inappropriate evidence classification.
6:                                                                         Incomplete negative control audit.
                                                                                  affected_script
                                                                                           <char>
1: 06_scripts/python/14_prepare_phase9b1_bulk_data.py; 06_scripts/R/14_phase9b1_bulk_validation.R
2: 06_scripts/python/14_prepare_phase9b1_bulk_data.py; 06_scripts/R/14_phase9b1_bulk_validation.R
3: 06_scripts/python/14_prepare_phase9b1_bulk_data.py; 06_scripts/R/14_phase9b1_bulk_validation.R
4: 06_scripts/python/14_prepare_phase9b1_bulk_data.py; 06_scripts/R/14_phase9b1_bulk_validation.R
5: 06_scripts/python/14_prepare_phase9b1_bulk_data.py; 06_scripts/R/14_phase9b1_bulk_validation.R
6: 06_scripts/python/14_prepare_phase9b1_bulk_data.py; 06_scripts/R/14_phase9b1_bulk_validation.R
                                                     affected_output
                                                              <char>
1:                                                purist_probability
2:                                                purist_probability
3: WGCNA Modules (black, blue, green, greenyellow, purple, red, tan)
4:              HALLMARK_PROTEIN_SECRETION, HALLMARK_SPERMATOGENESIS
5:                                          34 Transcription Factors
6:                                        Unrelated pathway controls
                                                                                                                                               correction_applied
                                                                                                                                                           <char>
1:                                    Include the intercept (-6.815) in the PurIST probability formula: prob = 1 / (1 + exp(-(intercept + sum(coef * (A > B))))).
2:                                                                      Modify the script to allow calculation and rescaling when coverage is >= 80% (7/8 pairs).
3:        Enforce the 80% coverage threshold. Flag WGCNA module validation in GSE71729 and GSE62452 as INSUFFICIENT_EXTERNAL_DATA and do not pool or report them.
4:                                                 Run the locked ssGSEA decoupleR method on the full MSigDB Hallmark gene sets (96 and 135 genes, respectively).
5: Rerun the VIPER algorithm via decoupleR using DoRothEA regulons (A/B/C) in R, or classify all TFs as TO_VERIFY and document that TF activity was not computed.
6:                                                                                         Compute scores and associations for the 5 unrelated Hallmark pathways.
                                outputs_recalculated
                                              <char>
1: phase9b1r_* tables, figures, and corrected report
2: phase9b1r_* tables, figures, and corrected report
3: phase9b1r_* tables, figures, and corrected report
4: phase9b1r_* tables, figures, and corrected report
5: phase9b1r_* tables, figures, and corrected report
6: phase9b1r_* tables, figures, and corrected report
                                                                                                              outputs_invalidated
                                                                                                                           <char>
1: phase9b1_* PurIST, Hallmark proxy, TF proxy, module low-coverage replication, cross-cohort synthesis, evidence classifications
2: phase9b1_* PurIST, Hallmark proxy, TF proxy, module low-coverage replication, cross-cohort synthesis, evidence classifications
3: phase9b1_* PurIST, Hallmark proxy, TF proxy, module low-coverage replication, cross-cohort synthesis, evidence classifications
4: phase9b1_* PurIST, Hallmark proxy, TF proxy, module low-coverage replication, cross-cohort synthesis, evidence classifications
5: phase9b1_* PurIST, Hallmark proxy, TF proxy, module low-coverage replication, cross-cohort synthesis, evidence classifications
6: phase9b1_* PurIST, Hallmark proxy, TF proxy, module low-coverage replication, cross-cohort synthesis, evidence classifications
   changes_scientific_conclusion
                          <char>
1:                           YES
2:                            NO
3:                           YES
4:                           YES
5:                           YES
6:                           YES

FIND_05 is now fully corrected: the executor derives TF evidence categories from the saved VIPER replication statistics and matches the locked Phase 9B1C2 audit counts (12 externally replicated, 13 partially replicated, 9 not replicated, 0 TO_VERIFY).
