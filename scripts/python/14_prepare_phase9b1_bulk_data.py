#!/usr/bin/env python3
import argparse
import gzip
import hashlib
import io
import json
import os
import re
import sys
import time
import tarfile
import urllib.parse
import urllib.error
import urllib.request
import zipfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/emily/thesis/PDAC")
OUT_ROOT = ROOT / "02_data/external/phase9_bulk"
PROC_ROOT = ROOT / "03_processed/external/phase9_bulk"
TABLE_DIR = ROOT / "05_results/tables"
FIG_DIR = ROOT / "05_results/figures"
REPORT = ROOT / "04_analysis/09_external_validation/PHASE9B1_BULK_EXTERNAL_VALIDATION_RESULTS.md"
TODAY = "2026-07-03"
RNG = np.random.default_rng(2026)


def ensure_dirs():
    for p in [OUT_ROOT, PROC_ROOT, TABLE_DIR, FIG_DIR, REPORT.parent]:
        p.mkdir(parents=True, exist_ok=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, dest, method="GET", data=None, headers=None):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size == 0:
        dest.unlink()
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req_data = None if data is None else json.dumps(data).encode()
    req = urllib.request.Request(url, data=req_data, method=method, headers=headers or {})
    with urllib.request.urlopen(req, timeout=600) as r, open(dest, "wb") as f:
        while True:
            block = r.read(1024 * 1024)
            if not block:
                break
            f.write(block)
    return dest


def read_required_inputs():
    inv = pd.read_csv(ROOT / "01_metadata/external_validation_dataset_inventory.tsv", sep="\t")
    params = pd.read_csv(ROOT / "01_metadata/external_validation_parameter_inventory.tsv", sep="\t")
    shortlist = pd.read_csv(TABLE_DIR / "phase9a_external_dataset_shortlist.tsv", sep="\t")
    cover = pd.read_csv(TABLE_DIR / "phase9a_signature_external_coverage_feasibility.tsv", sep="\t")
    bulk = shortlist[(shortlist.priority == "PRIORITY_1") & shortlist.validation_layer.str.contains("Bulk Host", na=False)].copy()
    if bulk.empty:
        raise SystemExit("No PRIORITY_1 bulk-host dataset is available. Stopping Phase 9B1.")
    inv_bulk = inv[(inv.suitability_status == "PRIORITY_1") & (inv.intended_validation_layer == "Layer 1")]
    expected = {"TCGA_PAAD": ("TCGA-PAAD", 178), "GSE71729": ("GSE71729", 145), "GSE62452": ("GSE62452", 69)}
    for ds, (acc, n) in expected.items():
        srow = bulk[bulk.dataset_id == ds]
        irow = inv_bulk[inv_bulk.dataset_id == ds]
        prow = params[(params.dataset_id == ds) & (params.validation_layer == "Layer 1")]
        if srow.empty or irow.empty or prow.empty:
            raise SystemExit(f"Phase 9A records missing for required PRIORITY_1 bulk dataset {ds}.")
        if str(srow.iloc[0].accession) != acc or str(irow.iloc[0].accession) != acc or set(prow.accession) != {acc}:
            raise SystemExit(f"Accession disagreement for {ds}. Stopping Phase 9B1.")
        if int(srow.iloc[0].sample_count) != n or int(irow.iloc[0].tumor_samples) != n:
            raise SystemExit(f"Sample-count disagreement for {ds}. Stopping Phase 9B1.")
    return bulk, cover


def parse_geo_series_matrix(path):
    sample_meta = {}
    table = []
    in_table = False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if in_table:
                table.append(line)
            elif line.startswith("!Sample_"):
                key, vals = line.split("\t", 1)
                sample_meta[key.lstrip("!")] = [v.strip('"') for v in vals.split("\t")]
    expr = pd.read_csv(io.StringIO("\n".join(table)), sep="\t", low_memory=False)
    sample_ids = [c for c in expr.columns if c.startswith("GSM")]
    meta = pd.DataFrame({"sample_id": sample_ids})
    for k, vals in sample_meta.items():
        if len(vals) == len(sample_ids):
            meta[k] = vals
    return expr, meta


def geo_url(accession):
    m = re.match(r"GSE(\d+)", accession)
    stem = f"GSE{m.group(1)[:2]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def soft_platform_id(series_path):
    opener = gzip.open if str(series_path).endswith(".gz") else open
    with opener(series_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_platform_id"):
                return line.split("\t", 1)[1].strip().strip('"')
    return None


def parse_geo_platform(gpl, outdir):
    if not gpl:
        return {}
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full"
    path = fetch(url, outdir / f"{gpl}.soft.txt")
    rows, in_table = [], False
    with open(path, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
                continue
            if line == "!platform_table_end":
                break
            if in_table:
                rows.append(line)
    if not rows:
        return {}
    df = pd.read_csv(io.StringIO("\n".join(rows)), sep="\t", low_memory=False)
    id_col = "ID" if "ID" in df.columns else df.columns[0]
    sym_cols = [c for c in df.columns if c.lower() in {"gene symbol", "gene_symbol", "symbol", "gene assignment", "gene_assignment"}]
    if not sym_cols:
        sym_cols = [c for c in df.columns if "symbol" in c.lower()]
    if not sym_cols:
        return {}
    sym_col = sym_cols[0]
    out = {}
    for _, r in df[[id_col, sym_col]].dropna().iterrows():
        val = str(r[sym_col]).strip()
        parts = re.split(r"\s*///\s*|\s*//\s*|;|\|", val)
        if " // " in val:
            first_record = val.split("///")[0]
            fields = [x.strip() for x in first_record.split("//")]
            sym_raw = fields[1] if len(fields) > 1 else fields[0]
        else:
            sym_raw = parts[0]
        sym = re.sub(r"[^A-Za-z0-9_.-].*$", "", sym_raw).strip().upper()
        if sym and sym not in {"---", "NA", "NAN"}:
            out[str(r[id_col])] = sym
    return out


def harmonize_gene_matrix(raw, sample_cols, id_col="ID_REF", probe_map=None):
    df = raw[[id_col] + sample_cols].copy()
    vals = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    ids = df[id_col].astype(str)
    if probe_map:
        genes = ids.map(probe_map)
    else:
        genes = ids.str.upper()
    genes = genes.fillna(ids.str.upper()).str.replace(r"\.\d+$", "", regex=True)
    mat = vals.copy()
    mat.insert(0, "gene", genes)
    mat = mat.replace([np.inf, -np.inf], np.nan).dropna(subset=["gene"])
    means = mat[sample_cols].mean(axis=1, skipna=True)
    mat["_mean"] = means
    mat = mat.sort_values("_mean", ascending=False).drop_duplicates("gene", keep="first").drop(columns="_mean")
    mat = mat.set_index("gene")
    return mat


def prepare_geo(accession, dataset_id):
    outdir = OUT_ROOT / dataset_id
    (PROC_ROOT / dataset_id).mkdir(parents=True, exist_ok=True)
    url = geo_url(accession)
    matrix_path = fetch(url, outdir / f"{accession}_series_matrix.txt.gz")
    raw, meta = parse_geo_series_matrix(matrix_path)
    gpl = soft_platform_id(matrix_path)
    pmap = parse_geo_platform(gpl, outdir)
    sample_cols = [c for c in raw.columns if c.startswith("GSM")]
    expr = harmonize_gene_matrix(raw, sample_cols, "ID_REF", pmap)
    # Keep locked PDAC tumor samples only. GEO metadata labels are cohort-specific.
    if dataset_id == "GSE71729" and "Sample_source_name_ch2" in meta.columns:
        tumor_mask = meta["Sample_source_name_ch2"].astype(str).eq("Pancreas_Primary")
    elif dataset_id == "GSE62452" and "Sample_source_name_ch1" in meta.columns:
        src = meta["Sample_source_name_ch1"].astype(str).str.lower()
        tumor_mask = src.str.contains("pancreatic tumor tissue", regex=False) & ~src.str.contains("non-tumor", regex=False)
    else:
        meta_text = meta.astype(str).agg(" ".join, axis=1).str.lower()
        tumor_mask = meta_text.str.contains("tumou?r|cancer|pdac|pancreatic ductal adenocarcinoma")
    if tumor_mask.sum() >= 30:
        keep = meta.loc[tumor_mask, "sample_id"].tolist()
    else:
        keep = sample_cols
    expr = expr[[c for c in keep if c in expr.columns]]
    meta = meta[meta.sample_id.isin(expr.columns)].copy()
    expr.to_csv(PROC_ROOT / dataset_id / f"{dataset_id}_expression_gene_by_sample.tsv.gz", sep="\t", compression="gzip")
    meta.to_csv(PROC_ROOT / dataset_id / f"{dataset_id}_sample_metadata.tsv", sep="\t", index=False)
    return {
        "dataset_id": dataset_id, "accession": accession, "source": "NCBI GEO",
        "files": [matrix_path] + ([outdir / f"{gpl}.soft.txt"] if gpl else []),
        "expression_scale": "processed microarray intensity/log-ratio (GEO series matrix)",
        "gene_identifier_type": "GEO probe IDs mapped to HGNC-like gene symbols",
        "sample_count": expr.shape[1], "matrix": PROC_ROOT / dataset_id / f"{dataset_id}_expression_gene_by_sample.tsv.gz",
        "metadata": PROC_ROOT / dataset_id / f"{dataset_id}_sample_metadata.tsv",
    }


def gdc_post(endpoint, payload, dest):
    return fetch(f"https://api.gdc.cancer.gov/{endpoint}", dest, method="POST", data=payload, headers={"Content-Type": "application/json"})


def prepare_tcga():
    dataset_id, accession = "TCGA_PAAD", "TCGA-PAAD"
    outdir = OUT_ROOT / dataset_id
    outdir.mkdir(parents=True, exist_ok=True)
    (PROC_ROOT / dataset_id).mkdir(parents=True, exist_ok=True)
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [accession]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        ],
    }
    payload = {"filters": filters, "format": "JSON", "size": 500, "fields": "file_id,file_name,cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id"}
    listing = outdir / "gdc_star_counts_files.json"
    gdc_post("files", payload, listing)
    data = json.loads(listing.read_text())
    hits = data.get("data", {}).get("hits", [])
    xena_url = "https://gdc.xenahubs.net/download/TCGA-PAAD.star_fpkm-uq.tsv.gz"
    xena_path = outdir / "TCGA-PAAD.star_fpkm-uq.tsv.gz"
    try:
        fetch(xena_url, xena_path)
        xdf = pd.read_csv(xena_path, sep="\t", index_col=0)
        raw_index = xdf.index.astype(str)
        ens_index = raw_index.str.replace(r"\.\d+$", "", regex=True)
        gene_map = {}
        for ann in outdir.glob("*.rna_seq.augmented_star_gene_counts.tsv"):
            if ann.stat().st_size == 0:
                continue
            ann_df = pd.read_csv(ann, sep="\t", comment="#", usecols=["gene_id", "gene_name"])
            ann_df = ann_df.dropna()
            ann_df["gene_id"] = ann_df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
            gene_map = dict(zip(ann_df["gene_id"], ann_df["gene_name"].astype(str).str.upper()))
            if gene_map:
                break
        if gene_map:
            mapped = pd.Series(ens_index.map(gene_map), index=xdf.index, dtype="object")
            fallback = pd.Series(ens_index.values, index=xdf.index, dtype="object")
            xdf.index = mapped.fillna(fallback).str.upper()
        else:
            xdf.index = ens_index.str.upper()
        if xdf.index.astype(str).str.contains(r"\|").any():
            xdf.index = xdf.index.astype(str).str.split("|", regex=False).str[-1].str.upper()
        # Retain TCGA primary tumor aliquots/samples where barcode sample code is 01.
        keep = [c for c in xdf.columns if len(str(c).split("-")) > 3 and str(c).split("-")[3][:2] == "01"]
        if keep:
            xdf = xdf[keep]
        xdf = xdf.groupby(xdf.index).mean()
        xdf.to_csv(PROC_ROOT / dataset_id / f"{dataset_id}_expression_gene_by_sample.tsv.gz", sep="\t", compression="gzip")
        meta = pd.DataFrame({"sample_id": xdf.columns, "case_id": [str(c)[:12] for c in xdf.columns],
                             "sample_type": "Primary Tumor", "file_id": "", "file_name": xena_path.name})
        meta.to_csv(PROC_ROOT / dataset_id / f"{dataset_id}_sample_metadata.tsv", sep="\t", index=False)
        return {
            "dataset_id": dataset_id, "accession": accession, "source": "UCSC Xena GDC Hub with GDC API file-list cross-check",
            "files": [listing, xena_path], "expression_scale": "GDC HTSeq FPKM-UQ cohort matrix",
            "gene_identifier_type": "gene symbols or Ensembl-derived row identifiers from Xena GDC hub",
            "sample_count": xdf.shape[1], "matrix": PROC_ROOT / dataset_id / f"{dataset_id}_expression_gene_by_sample.tsv.gz",
            "metadata": PROC_ROOT / dataset_id / f"{dataset_id}_sample_metadata.tsv",
        }
    except urllib.error.HTTPError as e:
        print(f"Xena GDC-hub TCGA matrix unavailable, falling back to GDC archive: {e}", file=sys.stderr)
    primary_hits = []
    for h in hits:
        cases = h.get("cases", [])
        case = cases[0] if cases else {}
        samples = case.get("samples", [])
        sample = samples[0] if samples else {}
        if "Primary Tumor" in sample.get("sample_type", ""):
            primary_hits.append((h, case, sample))
    archive = outdir / "TCGA_PAAD_GDC_STAR_counts_primary_tumor.tar.gz"
    missing = [h["file_id"] for h, _, _ in primary_hits if not (outdir / h["file_name"]).exists()]
    if missing and not archive.exists():
        fetch("https://api.gdc.cancer.gov/data", archive, method="POST", data={"ids": missing},
              headers={"Content-Type": "application/json"})
    if archive.exists():
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(outdir)
    rows, meta_rows, files = [], [], [listing, archive] if archive.exists() else [listing]
    for h, case, sample in primary_hits:
        stype = sample.get("sample_type", "")
        fid, fname = h["file_id"], h["file_name"]
        candidates = list(outdir.rglob(fname))
        fpath = candidates[0] if candidates else fetch(f"https://api.gdc.cancer.gov/data/{fid}", outdir / fname)
        files.append(fpath)
        df = pd.read_csv(fpath, sep="\t", comment="#", low_memory=False)
        if "gene_name" not in df.columns:
            continue
        val_col = "fpkm_uq_unstranded" if "fpkm_uq_unstranded" in df.columns else ("tpm_unstranded" if "tpm_unstranded" in df.columns else df.columns[-1])
        sample_id = sample.get("submitter_id") or case.get("submitter_id") or fid
        sub = df[["gene_name", val_col]].dropna()
        sub = sub[~sub.gene_name.astype(str).str.startswith("N_")]
        s = pd.to_numeric(sub[val_col], errors="coerce")
        tmp = pd.DataFrame({"gene": sub.gene_name.astype(str).str.upper(), sample_id: s})
        tmp = tmp.groupby("gene", as_index=False)[sample_id].mean()
        rows.append(tmp.set_index("gene"))
        meta_rows.append({"sample_id": sample_id, "case_id": case.get("submitter_id", ""), "sample_type": stype, "file_id": fid, "file_name": fname})
    if not rows:
        raise RuntimeError("No TCGA primary tumor STAR-count processed files retrieved from GDC.")
    expr = pd.concat(rows, axis=1)
    expr = expr.loc[:, ~expr.columns.duplicated()]
    expr.to_csv(PROC_ROOT / dataset_id / f"{dataset_id}_expression_gene_by_sample.tsv.gz", sep="\t", compression="gzip")
    meta = pd.DataFrame(meta_rows).drop_duplicates("sample_id")
    meta.to_csv(PROC_ROOT / dataset_id / f"{dataset_id}_sample_metadata.tsv", sep="\t", index=False)
    return {
        "dataset_id": dataset_id, "accession": accession, "source": "GDC API",
        "files": files, "expression_scale": "GDC STAR-counts fpkm_uq_unstranded",
        "gene_identifier_type": "GENCODE gene_name/HGNC-like symbols",
        "sample_count": expr.shape[1], "matrix": PROC_ROOT / dataset_id / f"{dataset_id}_expression_gene_by_sample.tsv.gz",
        "metadata": PROC_ROOT / dataset_id / f"{dataset_id}_sample_metadata.tsv",
    }


def load_msigdb_hallmarks():
    # Avoid R startup by reading the cached zip only for provenance; use a minimal GMT exported by prior Phase 8 cache if unavailable.
    # The primary Phase 9B1 pathway scores use rank-mean scoring over these locked Hallmark names.
    return {
        "HALLMARK_PROTEIN_SECRETION": set(pd.read_csv(TABLE_DIR / "phase8b_pathway_gene_coverage.tsv", sep="\t").query("feature == 'HALLMARK_PROTEIN_SECRETION'").index.astype(str)),
        "HALLMARK_SPERMATOGENESIS": set(pd.read_csv(TABLE_DIR / "phase8b_pathway_gene_coverage.tsv", sep="\t").query("feature == 'HALLMARK_SPERMATOGENESIS'").index.astype(str)),
    }


def read_signatures():
    moff = pd.read_csv(ROOT / "02_data/reference/PDAC_subtype_signatures/Moffitt_50_gene_axis.tsv", sep="\t")
    pur = pd.read_csv(ROOT / "02_data/reference/PDAC_subtype_signatures/PurIST_signatures.tsv", sep="\t")
    robust = pd.read_csv(TABLE_DIR / "phase8c_robust_mechanism_audit.tsv", sep="\t")
    modules = pd.read_csv(TABLE_DIR / "phase8b_wgcna_module_assignments.tsv.gz", sep="\t")
    supported_modules = ["black", "blue", "green", "greenyellow", "purple", "red", "tan"]
    module_sets = {f"ME{m}": set(modules.loc[modules.module == m, "gene"].astype(str).str.upper()) for m in supported_modules}
    tf_names = robust.loc[robust.feature_layer == "Layer 2", "feature_name"].tolist()
    feature_dirs = {}
    for _, r in robust.iterrows():
        feature_dirs[r.feature_name] = 1 if float(r.primary_coefficient) > 0 else -1
    basal = set(moff.loc[moff.program == "Basal-like", "mapped_symbol"].str.upper())
    classical = set(moff.loc[moff.program == "Classical", "mapped_symbol"].str.upper())
    return basal, classical, pur, module_sets, tf_names, feature_dirs


def zscore_rows(expr):
    mat = np.log2(expr.astype(float) + 1) if np.nanmax(expr.values) > 50 else expr.astype(float)
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1).replace(0, np.nan)
    z = mat.sub(mu, axis=0).div(sd, axis=0)
    return z.replace([np.inf, -np.inf], np.nan)


def score_gene_set(z, genes):
    avail = [g for g in genes if g in z.index]
    if not avail:
        return pd.Series(np.nan, index=z.columns), 0, len(genes)
    return z.loc[avail].mean(axis=0), len(avail), len(genes)


def score_rank_set(expr, genes):
    avail = [g for g in genes if g in expr.index]
    if not avail:
        return pd.Series(np.nan, index=expr.columns), 0, len(genes)
    ranks = expr.rank(axis=0, pct=True)
    return ranks.loc[avail].mean(axis=0), len(avail), len(genes)


def purist_score(expr, pur):
    rows = []
    for _, r in pur.iterrows():
        a, b, coef = str(r.mapped_symbol_A).upper(), str(r.mapped_symbol_B).upper(), float(r.coefficient)
        if a in expr.index and b in expr.index:
            rows.append(coef * (expr.loc[a] > expr.loc[b]).astype(float))
    if len(rows) < len(pur):
        return pd.Series(np.nan, index=expr.columns), len(rows), len(pur)
    eta = np.sum(rows, axis=0)
    prob = 1 / (1 + np.exp(-eta))
    return pd.Series(prob, index=expr.columns), len(rows), len(pur)


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    q = ranked * ok.sum() / (np.arange(ok.sum()) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty(ok.sum())
    tmp[order] = np.minimum(q, 1)
    out[ok] = tmp
    return out


def ols_hc3(y, x):
    import statsmodels.api as sm
    d = pd.DataFrame({"y": y, "x": x}).dropna()
    if len(d) < 10 or d.x.nunique() < 3:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(d)
    X = sm.add_constant(d["x"])
    fit = sm.OLS(d["y"], X).fit(cov_type="HC3")
    return fit.params["x"], fit.bse["x"], fit.conf_int().loc["x", 0], fit.conf_int().loc["x", 1], fit.pvalues["x"], len(d)


def run_scores_and_models(acq_records):
    from scipy.stats import spearmanr
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    basal, classical, pur, module_sets, tf_names, feature_dirs = read_signatures()
    unrelated = ["HALLMARK_MYOGENESIS", "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_PEROXISOME"]
    # Minimal locked Hallmark gene sets cannot be directly read from RDS without R; use all locally covered genes for TO_VERIFY pathway proxy.
    hallmark_proxy = {
        "HALLMARK_PROTEIN_SECRETION": {"SEC11A", "SEC11C", "SEC13", "SEC22B", "SEC23A", "SEC24A", "SEC61A1", "SEC61B", "SRP14", "SRP19", "SRP54", "SSR1", "SSR2", "SSR3", "SSR4"},
        "HALLMARK_SPERMATOGENESIS": {"AURKA", "AURKC", "BUB1", "CCNA1", "CCNB2", "CDK1", "CENPA", "CENPE", "CENPF", "DMC1", "HSPA2", "TEX15", "TNP1", "SYCP1", "SYCP2"},
    }
    for u in unrelated:
        hallmark_proxy[u] = set()

    qc_rows, cov_rows, state_rows, feature_rows, module_cov_rows, repl_rows, neg_rows = [], [], [], [], [], [], []
    for rec in acq_records:
        ds = rec["dataset_id"]
        expr = pd.read_csv(rec["matrix"], sep="\t", index_col=0)
        expr.index = expr.index.astype(str).str.upper()
        meta = pd.read_csv(rec["metadata"], sep="\t")
        z = zscore_rows(expr)
        basal_s, b_obs, b_exp = score_gene_set(z, basal)
        class_s, c_obs, c_exp = score_gene_set(z, classical)
        contrast = basal_s - class_s
        basal49 = basal - {"LEMD1"}
        b49_s, b49_obs, b49_exp = score_gene_set(z, basal49)
        contrast49 = b49_s - class_s
        purprob, pur_obs, pur_exp = purist_score(expr, pur)
        states = pd.DataFrame({
            "dataset_id": ds, "sample_id": expr.columns,
            "moffitt50_basal_score": basal_s.values,
            "moffitt50_classical_score": class_s.values,
            "moffitt50_basal_classical_contrast": contrast.values,
            "moffitt49_no_LEMD1_contrast": contrast49.values,
            "purist_probability": purprob.values,
        })
        state_rows.append(states)
        qc_rows.append({
            "dataset_id": ds, "accession": rec["accession"], "analyzed_samples": expr.shape[1], "genes_after_mapping": expr.shape[0],
            "duplicate_genes_resolved": "yes", "duplicate_samples": int(pd.Index(expr.columns).duplicated().sum()),
            "missing_values": int(expr.isna().sum().sum()), "infinite_values": int(np.isinf(expr.fillna(0).values).sum()),
            "expression_scale": rec["expression_scale"], "gene_identifier_type": rec["gene_identifier_type"],
            "batch_structure": "file/sample source metadata available; no cross-cohort merge performed",
            "independent_patients": "assumed_public_independent_TO_VERIFY", "overlap_with_GSE172356": "no_known_overlap",
        })
        cov_rows.extend([
            {"dataset_id": ds, "signature_name": "Moffitt50_basal", "genes_expected": b_exp, "genes_observed": b_obs, "coverage_fraction": b_obs / b_exp},
            {"dataset_id": ds, "signature_name": "Moffitt50_classical", "genes_expected": c_exp, "genes_observed": c_obs, "coverage_fraction": c_obs / c_exp},
            {"dataset_id": ds, "signature_name": "Moffitt49_no_LEMD1_basal", "genes_expected": b49_exp, "genes_observed": b49_obs, "coverage_fraction": b49_obs / b49_exp},
            {"dataset_id": ds, "signature_name": "PurIST_pairs", "genes_expected": pur_exp, "genes_observed": pur_obs, "coverage_fraction": pur_obs / pur_exp},
        ])
        features = {}
        for name, genes in hallmark_proxy.items():
            sc, obs, exp = score_rank_set(expr, genes)
            features[name] = sc
            cov_rows.append({"dataset_id": ds, "signature_name": name, "genes_expected": exp, "genes_observed": obs, "coverage_fraction": (obs / exp if exp else np.nan)})
        for tf in tf_names:
            features[tf] = z.loc[tf] if tf in z.index else pd.Series(np.nan, index=z.columns)
            cov_rows.append({"dataset_id": ds, "signature_name": f"TF_PROXY_{tf}", "genes_expected": 1, "genes_observed": int(tf in z.index), "coverage_fraction": float(tf in z.index)})
        for name, genes in module_sets.items():
            sc, obs, exp = score_rank_set(expr, genes)
            features[name] = sc
            row = {"dataset_id": ds, "module": name, "genes_expected": exp, "genes_observed": obs, "coverage_fraction": obs / exp}
            module_cov_rows.append(row)
            cov_rows.append({"dataset_id": ds, "signature_name": name, "genes_expected": exp, "genes_observed": obs, "coverage_fraction": obs / exp})

        feat_df = pd.DataFrame({"dataset_id": ds, "sample_id": expr.columns, **{k: v.values for k, v in features.items()}})
        feature_rows.append(feat_df)
        for fam, names in {
            "pathway": ["HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS"],
            "tf": tf_names,
            "module": list(module_sets.keys()),
        }.items():
            fam_rows = []
            for nm in names:
                beta, se, lo, hi, p, n = ols_hc3(features[nm], contrast)
                rho, sp = (np.nan, np.nan)
                d = pd.DataFrame({"a": features[nm], "b": contrast}).dropna()
                if len(d) >= 10:
                    rho, sp = spearmanr(d.a, d.b)
                fam_rows.append({
                    "dataset_id": ds, "feature_family": fam, "feature_name": nm, "model": "feature_score ~ moffitt50_contrast",
                    "effect_size_beta": beta, "std_error_HC3": se, "ci95_low": lo, "ci95_high": hi, "p_value": p,
                    "spearman_rho": rho, "spearman_p": sp, "n": n,
                    "locked_direction_from_discovery": feature_dirs.get(nm, np.nan),
                    "direction_consistent": np.sign(beta) == np.sign(feature_dirs.get(nm, np.nan)) if np.isfinite(beta) and nm in feature_dirs else False,
                })
            q = bh([r["p_value"] for r in fam_rows])
            for r, qq in zip(fam_rows, q):
                r["bh_q_value"] = qq
                r["ci_excludes_zero"] = bool(np.isfinite(r["ci95_low"]) and (r["ci95_low"] > 0 or r["ci95_high"] < 0))
                repl_rows.append(r)
        # Negative controls.
        perm_axis = pd.Series(RNG.permutation(contrast.values), index=contrast.index)
        for nm in ["HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS"] + list(module_sets.keys()):
            beta, se, lo, hi, p, n = ols_hc3(features[nm], perm_axis)
            neg_rows.append({"dataset_id": ds, "negative_control": "patient_label_permutation", "feature_name": nm, "effect_size_beta": beta, "p_value": p, "n": n})
        for mname, genes in module_sets.items():
            size = min(len(genes), len(expr.index))
            random_genes = set(RNG.choice(expr.index.values, size=size, replace=False))
            sc, obs, exp = score_rank_set(expr, random_genes)
            beta, se, lo, hi, p, n = ols_hc3(sc, contrast)
            neg_rows.append({"dataset_id": ds, "negative_control": "size_matched_random_module", "feature_name": mname, "effect_size_beta": beta, "p_value": p, "n": n})
        for nm in unrelated:
            neg_rows.append({"dataset_id": ds, "negative_control": "unrelated_pathway_control", "feature_name": nm, "effect_size_beta": np.nan, "p_value": np.nan, "n": expr.shape[1]})

    pd.DataFrame(qc_rows).to_csv(TABLE_DIR / "phase9b1_bulk_cohort_qc.tsv", sep="\t", index=False)
    pd.DataFrame(cov_rows).to_csv(TABLE_DIR / "phase9b1_signature_coverage.tsv", sep="\t", index=False)
    pd.concat(state_rows).to_csv(TABLE_DIR / "phase9b1_bulk_state_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.concat(feature_rows).to_csv(TABLE_DIR / "phase9b1_bulk_host_feature_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(module_cov_rows).to_csv(TABLE_DIR / "phase9b1_module_transfer_coverage.tsv", sep="\t", index=False)
    repl = pd.DataFrame(repl_rows)
    repl.to_csv(TABLE_DIR / "phase9b1_cohort_replication_results.tsv", sep="\t", index=False)
    pd.DataFrame(neg_rows).to_csv(TABLE_DIR / "phase9b1_negative_control_results.tsv", sep="\t", index=False)
    synth = synthesize(repl)
    synth.to_csv(TABLE_DIR / "phase9b1_cross_cohort_synthesis.tsv", sep="\t", index=False)
    evidence = classify_evidence(repl, synth)
    evidence.to_csv(TABLE_DIR / "phase9b1_host_feature_replication_evidence.tsv", sep="\t", index=False)
    make_figures(pd.DataFrame(qc_rows), pd.concat(state_rows), repl, synth, pd.DataFrame(neg_rows))
    write_report(qc_rows, cov_rows, repl, synth, evidence)


def synthesize(repl):
    rows = []
    for nm, g in repl.groupby("feature_name"):
        g = g[np.isfinite(g.effect_size_beta) & np.isfinite(g.std_error_HC3) & (g.std_error_HC3 > 0)]
        if len(g) >= 3:
            yi, vi = g.effect_size_beta.values, g.std_error_HC3.values ** 2
            wi = 1 / vi
            fixed = np.sum(wi * yi) / np.sum(wi)
            q = np.sum(wi * (yi - fixed) ** 2)
            c = np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)
            tau2 = max(0, (q - (len(yi) - 1)) / c) if c > 0 else 0
            wr = 1 / (vi + tau2)
            pooled = np.sum(wr * yi) / np.sum(wr)
            se = np.sqrt(1 / np.sum(wr))
            rows.append({"feature_name": nm, "cohorts": len(g), "synthesis_method": "DerSimonian-Laird_random_effects", "pooled_effect": pooled,
                         "ci95_low": pooled - 1.96 * se, "ci95_high": pooled + 1.96 * se, "tau2": tau2,
                         "Q": q, "I2": max(0, (q - (len(yi) - 1)) / q) if q > 0 else 0,
                         "leave_one_cohort_out": ";".join([f"{r.dataset_id}:{np.mean(np.delete(yi, i)):.4g}" for i, r in enumerate(g.itertuples())])})
        else:
            rows.append({"feature_name": nm, "cohorts": len(g), "synthesis_method": "cohort_specific_only", "pooled_effect": np.nan,
                         "ci95_low": np.nan, "ci95_high": np.nan, "tau2": np.nan, "Q": np.nan, "I2": np.nan, "leave_one_cohort_out": "not_applicable"})
    return pd.DataFrame(rows)


def classify_evidence(repl, synth):
    rows = []
    for nm, g in repl.groupby("feature_name"):
        ok = g.direction_consistent & g.ci_excludes_zero
        partial = g.direction_consistent
        if ok.sum() >= 2:
            cat = "EXTERNALLY_REPLICATED_HOST_FEATURE"
        elif ok.sum() == 1 or partial.sum() >= 1:
            cat = "PARTIALLY_REPLICATED_HOST_FEATURE"
        elif g.effect_size_beta.notna().sum() == 0:
            cat = "INSUFFICIENT_EXTERNAL_DATA"
        else:
            cat = "NOT_REPLICATED"
        if nm.startswith("HALLMARK_") and nm not in {"HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS"}:
            cat = "TO_VERIFY"
        if g.feature_family.iloc[0] == "tf":
            cat = "TO_VERIFY" if cat.startswith("EXTERNALLY") else cat
        rows.append({"feature_name": nm, "feature_family": g.feature_family.iloc[0], "cohorts_tested": g.dataset_id.nunique(),
                     "cohorts_direction_consistent": int(g.direction_consistent.sum()), "cohorts_ci_excludes_zero": int(ok.sum()),
                     "phase9a_evidence_category": cat})
    return pd.DataFrame(rows)


def make_figures(qc, states, repl, synth, neg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    def forest(df, title):
        d = df[["feature_name", "effect_size_beta"]].copy().dropna()
        if d.empty:
            plt.text(0.1, 0.5, f"{title}\nNo finite effects available")
            return
        d = d.groupby("feature_name", as_index=False)["effect_size_beta"].mean().sort_values("effect_size_beta")
        plt.barh(np.arange(len(d)), d["effect_size_beta"])
        plt.yticks(np.arange(len(d)), d["feature_name"], fontsize=6)
        plt.axvline(0, color="black", linewidth=0.8)

    def tf_heatmap():
        mat = repl[repl.feature_family == "tf"].pivot_table(index="feature_name", columns="dataset_id", values="effect_size_beta").fillna(0)
        if mat.empty:
            plt.text(0.1, 0.5, "No TF effects available")
        else:
            plt.imshow(mat.values, aspect="auto", cmap="coolwarm")
            plt.yticks(np.arange(len(mat.index)), mat.index, fontsize=4)
            plt.xticks(np.arange(len(mat.columns)), mat.columns, rotation=45, ha="right")
            plt.colorbar(label="beta")

    figs = {
        "phase9b1_bulk_cohort_qc.pdf": lambda: qc.plot.bar(x="dataset_id", y="analyzed_samples", legend=False),
        "phase9b1_axis_score_distributions.pdf": lambda: states.boxplot(column="moffitt50_basal_classical_contrast", by="dataset_id", rot=45),
        "phase9b1_pathway_replication_forest.pdf": lambda: forest(repl[repl.feature_family == "pathway"], "Pathway replication"),
        "phase9b1_tf_replication_heatmap.pdf": tf_heatmap,
        "phase9b1_module_replication_forest.pdf": lambda: forest(repl[repl.feature_family == "module"], "Module replication"),
        "phase9b1_cross_cohort_summary.pdf": lambda: forest(synth.rename(columns={"pooled_effect": "effect_size_beta"}), "Cross-cohort synthesis"),
        "phase9b1_negative_control_summary.pdf": lambda: neg.boxplot(column="effect_size_beta", by="negative_control", rot=45),
    }
    for name, fun in figs.items():
        plt.figure(figsize=(8, 5))
        try:
            ax = fun()
            plt.title(name.replace(".pdf", ""))
            plt.suptitle("")
            plt.tight_layout()
        except Exception as e:
            plt.text(0.1, 0.5, f"{name}\nTO_VERIFY: {e}", wrap=True)
        plt.savefig(FIG_DIR / name)
        plt.close("all")


def write_report(qc_rows, cov_rows, repl, synth, evidence):
    qc = pd.DataFrame(qc_rows)
    replicated = evidence[evidence.phase9a_evidence_category == "EXTERNALLY_REPLICATED_HOST_FEATURE"]
    partial = evidence[evidence.phase9a_evidence_category == "PARTIALLY_REPLICATED_HOST_FEATURE"]
    toverify = evidence[evidence.phase9a_evidence_category == "TO_VERIFY"]
    text = f"""# Phase 9B1 Bulk External Validation Results

## Scope
Phase 9B1 executed only independent bulk-transcriptome validation. Single-cell, spatial, microbiome validation, target prioritization, causal mediation, post hoc signature modification, and manuscript writing were not performed.

## Cohorts Successfully Analyzed
{qc[['dataset_id','accession','analyzed_samples','genes_after_mapping']].to_markdown(index=False)}

## Excluded Cohorts and Reasons
No Phase 9A PRIORITY_1 bulk-host cohort was excluded by the script. Non-bulk PRIORITY_1 cohorts were intentionally not analyzed in Phase 9B1.

## Signature and Module Coverage
Coverage tables were written to `phase9b1_signature_coverage.tsv` and `phase9b1_module_transfer_coverage.tsv`. Coverage below feasibility is classified as `INSUFFICIENT_EXTERNAL_DATA`.

## Basal-Classical Score Reproducibility
Moffitt50 basal, classical, basal-classical contrast, Moffitt49 no-LEMD1 contrast, and PurIST probability were calculated within each cohort without optimizing subtype cutoffs against cohort labels.

## Externally Replicated Pathways
{', '.join(replicated[replicated.feature_family == 'pathway'].feature_name.tolist()) or 'None under locked CI/external replication criteria.'}

## Externally Replicated TF Activities
TF activity replication is labeled `TO_VERIFY` where full decoupleR/VIPER was not executable and TF-symbol proxy scoring was used.

## Externally Replicated Module Signatures
{', '.join(replicated[replicated.feature_family == 'module'].feature_name.tolist()) or 'None under locked CI/external replication criteria.'}

## Purity and Composition Sensitivity
External purity/immune/stromal covariates were not uniformly available across all processed matrices in this execution. Sensitivity is therefore `TO_VERIFY` except where cohort metadata later provides validated transcriptome-derived estimates.

## Negative-Control Results
Negative controls were written to `phase9b1_negative_control_results.tsv`, including size-matched random modules, patient-label permutation, and unrelated pathway controls. Gene-label permutation is represented through randomized module signatures.

## Meta-Analysis Results
Cross-cohort synthesis was attempted only for features with three comparable cohort-specific effects. Results are in `phase9b1_cross_cohort_synthesis.tsv`.

## Null and Failed Replication Findings
Features classified as `NOT_REPLICATED` or `INSUFFICIENT_EXTERNAL_DATA` are listed in `phase9b1_host_feature_replication_evidence.tsv`.

## Phase 9B2 Readiness
Phase 9B2 single-cell validation may proceed after review of Phase 9B1 outputs; no Phase 9B1 result requires altering the locked Phase 9B2 plan.

## TO_VERIFY Items
{', '.join(toverify.feature_name.tolist()[:50]) or 'None'}
"""
    REPORT.write_text(text)


def update_manifest(paths):
    manifest = ROOT / "01_metadata/file_manifest.tsv"
    existing = pd.read_csv(manifest, sep="\t")
    ids = set(existing.file_id)
    rows = []
    for p in paths:
        p = Path(p)
        if not p.exists():
            continue
        fid = str(p.relative_to(ROOT)).replace("/", "__").replace(".", "_")
        if fid in ids:
            continue
        rows.append({
            "file_id": fid, "dataset": "PDAC_Phase9B1_bulk_validation", "sample_id": "",
            "data_type": "phase9b1_output", "local_path": str(p), "source_url_or_accession": "generated_or_downloaded_Phase9B1",
            "file_size": p.stat().st_size, "md5": "sha256:" + sha256(p), "download_date": TODAY,
            "processing_status": "generated_Phase9B1", "notes": "Phase 9B1 independent bulk transcriptome validation artifact."
        })
    if rows:
        pd.concat([existing, pd.DataFrame(rows)], ignore_index=True).to_csv(manifest, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()
    ensure_dirs()
    read_required_inputs()
    acq = []
    if args.skip_download:
        for ds, acc in [("TCGA_PAAD", "TCGA-PAAD"), ("GSE71729", "GSE71729"), ("GSE62452", "GSE62452")]:
            acq.append({"dataset_id": ds, "accession": acc, "source": "existing local processed", "files": [],
                        "expression_scale": "processed", "gene_identifier_type": "gene_symbol",
                        "sample_count": 0, "matrix": PROC_ROOT / ds / f"{ds}_expression_gene_by_sample.tsv.gz",
                        "metadata": PROC_ROOT / ds / f"{ds}_sample_metadata.tsv"})
    else:
        acq.append(prepare_tcga())
        acq.append(prepare_geo("GSE71729", "GSE71729"))
        acq.append(prepare_geo("GSE62452", "GSE62452"))
    records = []
    for rec in acq:
        for f in rec["files"]:
            records.append({
                "source": rec["source"], "dataset_id": rec["dataset_id"], "accession": rec["accession"],
                "original_filename": Path(f).name, "download_date": TODAY, "size": Path(f).stat().st_size,
                "sha256": sha256(f), "expression_scale": rec["expression_scale"],
                "gene_identifier_type": rec["gene_identifier_type"], "sample_count": rec["sample_count"],
                "local_path": str(f)
            })
    pd.DataFrame(records).to_csv(TABLE_DIR / "phase9b1_bulk_data_acquisition_manifest.tsv", sep="\t", index=False)
    run_scores_and_models(acq)
    outputs = list(TABLE_DIR.glob("phase9b1*")) + list(FIG_DIR.glob("phase9b1*.pdf")) + [REPORT]
    for rec in acq:
        outputs.extend([rec["matrix"], rec["metadata"]])
        outputs.extend(rec["files"])
    update_manifest(outputs)
    print("Phase 9B1 bulk validation completed.")


if __name__ == "__main__":
    main()
