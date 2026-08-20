# PDAC Subtype Reproduction Protocol

This document defines the execution protocol and validation procedures for reproducing the transcriptomic subtypes of the GSE172356 cohort. 

---

## 1. Environment Setup

### R Dependencies
Ensure R (v4.5.3 or compatible) is available and the following packages are installed:
```R
install.packages(c("readxl", "pheatmap", "openpyxl"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")
devtools::install_github("adeschen/SignatureHeatmap")
```

### Python Dependencies
Ensure Python 3.10+ is available with pandas and openpyxl:
```bash
pip install pandas openpyxl
```

---

## 2. Verified Signature Files

All signatures are stored under `02_data/reference/PDAC_subtype_signatures/`:

1. **GSE172356 Original Signature (`GSE172356_original_signatures.tsv`)**
   - Source: Chan-Seng-Yue et al., 2020 (Supplementary Table 4)
   - Size: 94 genes (6 missing from the original 100-gene panel)
   - Columns: `gene_symbol`, `original_signature`, `rank_in_signature`, `presence_in_GSE172356`, `mapped_symbol`
   - Checksum: `md5: 1fa46a3ee02166880bc58639972199c2`

2. **Moffitt 2015 Signature (`Moffitt_2015_signatures.tsv`)**
   - Source: Moffitt et al., 2015 (SignatureHeatmap)
   - Size: 50 genes (25 basal-like, 25 classical; LEMD1 is excluded in R code, leaving 49 active genes)
   - Columns: `gene_symbol`, `class`, `presence_in_GSE172356`, `mapped_symbol`
   - Checksum: `md5: fa1ec8714ff73e152014d9d564ee222d`

3. **PurIST 2020 Signature (`PurIST_signatures.tsv`)**
   - Source: Rashid et al., 2020 (fitteds_public_2019-02-12.Rdata)
   - Size: 16 genes (8 Top Scoring Pairs)
   - Columns: `pair_index`, `gene_A`, `gene_B`, `coefficient`, `direction`, `presence_in_GSE172356`, `mapped_symbol_A`, `mapped_symbol_B`
   - Checksum: `md5: 066d543aaef11b82a755d15e408edf1f`

---

## 3. Preprocessing and Mapping Rules

### Symbol Mapping for Moffitt 2015
The following genes must be mapped before hierarchical clustering:
- `CTSL2` $\rightarrow$ `CTSV`
- `ANXA8L2` $\rightarrow$ `ANXA8`
- `ATAD4` $\rightarrow$ `FLAD1`
- `LOC400573` $\rightarrow$ `TMEM238L`
- `LEMD1` must be excluded.

### Data Center/Scale
For the **GSE172356_original** and **Moffitt** methods:
1. Load `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`.
2. Extract the signature genes (handling missing genes for GSE172356_original and mapping symbols for Moffitt).
3. Compute the row median for each gene across the 62 samples.
4. Subtract the row median from the normalized counts for each cell.
5. In the clustering step, apply standard row-scaling:
   $$X_{scaled} = \frac{X - \text{mean}(X)}{\text{sd}(X)}$$
   where $X$ represents the median-centered count.

For the **PurIST** method:
- Directly use the untransformed size-factor normalized counts (no centering or scaling).

---

## 4. Subtype Assignment Rules

### Primary Method (GSE172356_original)
1. Compute the column distance matrix using Pearson correlation:
   $$d(X, Y) = 1 - \text{cor}(X, Y)$$
2. Perform average linkage hierarchical clustering.
3. Order the columns (samples) based on the column dendrogram.
4. Assign the subtypes based on the sorted order:
   - Samples 1 to 17: **Basal**
   - Samples 18 to 40: **Hybrid**
   - Samples 41 to 62: **Classical**

### Secondary Method 1 (Moffitt)
1. Compute the column distance matrix using Pearson correlation.
2. Perform average linkage hierarchical clustering.
3. Order the columns based on the column dendrogram.
4. Assign the subtypes based on the sorted order:
   - Samples 1 to 27: **Classical**
   - Samples 28 to 44: **Basal**
   - Samples 45 to 62: **Others**

### Secondary Method 2 (PurIST)
1. Evaluate the 8 top-scoring pairs for each sample:
   - $P_1 = I(\text{GPR87} > \text{REG4})$
   - $P_2 = I(\text{KRT6A} > \text{ANXA10})$
   - $P_3 = I(\text{BCAR3} > \text{GATA6})$
   - $P_4 = I(\text{PTGES} > \text{CLDN18})$
   - $P_5 = I(\text{ITGA3} > \text{LGALS4})$
   - $P_6 = I(\text{C16orf74} > \text{DDC})$
   - $P_7 = I(\text{S100A2} > \text{SLC40A1})$
   - $P_8 = I(\text{KRT5} > \text{CLRN3})$
2. Calculate the log-odds score:
   $$S = -6.815 + 1.994 \cdot P_1 + 2.031 \cdot P_2 + 1.618 \cdot P_3 + 0.922 \cdot P_4 + 1.059 \cdot P_5 + 0.929 \cdot P_6 + 2.505 \cdot P_7 + 0.485 \cdot P_8$$
3. Convert to probability:
   $$p = \frac{1}{1 + e^{-S}}$$
4. Classify:
   - If $p > 0.5$, assign **Basal-like**
   - If $p \le 0.5$, assign **Classical**

---

## 5. Verification Checks

Before locking the reproduction:
- Confirm that the sample classification of the primary method matches the column `subtype_original` in `sample_manifest.tsv` with 100% accuracy.
- Record the concordance rate (exact agreement) and Cohen's kappa between the primary method and the two secondary methods.
- Run the analysis under the predefined sensitivity analyses and log the concordance rates.
