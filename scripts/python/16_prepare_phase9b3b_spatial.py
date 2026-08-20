#!/usr/bin/env python3
"""Prepare Phase 9B3B spatial validation inputs from official processed files."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import tarfile
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/emily/thesis/PDAC")
RAW = ROOT / "02_data/external/phase9_spatial"
OUT = ROOT / "03_processed/external/phase9_spatial"
TABLES = ROOT / "05_results/tables"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_extract_selected(tar_path: Path, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    wanted = []
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if not (
                name.endswith(".tsv.gz")
                or name.endswith("-filtered.txt.gz")
                or name.endswith(".txt.gz")
            ):
                continue
            if "indrop" in name.lower() or "fastq" in name.lower() or name.lower().endswith(".bam"):
                continue
            target = dest / name
            if not target.exists():
                source = tar.extractfile(member)
                if source is None:
                    continue
                with target.open("wb") as out:
                    out.write(source.read())
            wanted.append(target)
    return sorted(wanted)


def read_hwang() -> None:
    src = RAW / "GSE199102"
    dst = OUT / "GSE199102"
    dst.mkdir(parents=True, exist_ok=True)
    props_path = src / "GSE199102_Broad_PDAC_WTA_AllSamples_SegmentProperties.txt.gz"
    q3_path = src / "GSE199102_hPDAC_WTA_20210222T2101_Q3Norm_TargetCountMatrix.txt.gz"
    raw_path = src / "GSE199102_hPDAC_WTA_20210222T2101_TargetCountMatrix.txt.gz"
    props = pd.read_csv(props_path, sep="\t")
    props["sample_id_matrix"] = props["Sample_ID"].str.replace("-", ".", regex=False)
    q3_cols = pd.read_csv(q3_path, sep="\t", nrows=0).columns[1:].tolist()
    raw_counts = pd.read_csv(raw_path, sep="\t", index_col=0)
    raw_counts.columns = raw_counts.columns.astype(str).str.replace("-", ".", regex=False)
    props = props[props["sample_id_matrix"].isin(q3_cols)].copy()
    props["patient_id"] = props["Patient"]
    props["section_id"] = props["Slide_name"]
    props["ROI_id"] = props["patient_id"] + "_ROI" + props["ROI_number"].astype(str)
    props["segment_id"] = props["sample_id_matrix"]
    props["compartment"] = props["Segment"].map(
        {"Epithelial": "malignant_epithelial", "CAF": "fibroblast_CAF", "Immune": "immune"}
    ).fillna(props["Segment"])
    props["cohort_id"] = props["Status"].map(
        {"Untreated": "HWANG_GSE202051_NAIVE", "Treated": "HWANG_GSE202051_TREATED"}
    )
    props["paired_segment_id"] = props.groupby(["patient_id", "ROI_id"])["segment_id"].transform(
        lambda x: ";".join(sorted(x.astype(str)))
    )
    raw_counts = raw_counts.loc[:, [c for c in raw_counts.columns if c in set(props["segment_id"])]]
    detected = (raw_counts > 0).sum(axis=0)
    library = raw_counts.sum(axis=0)
    props["detected_genes"] = props["segment_id"].map(detected)
    props["library_size"] = props["segment_id"].map(library)
    props["passes_qc"] = (props["detected_genes"] >= 1000) & (props["library_size"] >= 1000)
    # ROI composition is derived from independent morphology segment areas, not target genes.
    area = props.pivot_table(index=["patient_id", "ROI_id"], columns="Segment", values="AOI_area", aggfunc="sum")
    area = area.div(area.sum(axis=1), axis=0).fillna(0)
    area = area.rename(columns={"CAF": "CAF_fraction", "Immune": "myeloid_fraction", "Epithelial": "epithelial_fraction"})
    props = props.merge(area.reset_index(), on=["patient_id", "ROI_id"], how="left")
    props["lymphoid_fraction"] = pd.NA
    props.to_csv(dst / "GSE199102_segment_metadata_prepared.tsv", sep="\t", index=False)
    manifest = []
    for path in sorted(src.glob("*")):
        if path.is_file():
            manifest.append(
                {
                    "dataset_id": "GSE199102",
                    "path": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "source": "NCBI GEO GSE199102 supplemental",
                }
            )
    pd.DataFrame(manifest).to_csv(dst / "GSE199102_download_manifest.tsv", sep="\t", index=False)


def prepare_moncada() -> None:
    src = RAW / "GSE111672"
    dst = OUT / "GSE111672"
    dst.mkdir(parents=True, exist_ok=True)
    tar_path = src / "GSE111672_RAW.tar"
    extracted = []
    if tar_path.exists() and tar_path.stat().st_size > 700_000_000:
        extracted = safe_extract_selected(tar_path, dst / "extracted")
    manifest = []
    for path in sorted(src.glob("*")) + extracted:
        if path.is_file():
            manifest.append(
                {
                    "dataset_id": "GSE111672",
                    "path": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "source": "NCBI GEO GSE111672 supplemental",
                }
            )
    pd.DataFrame(manifest).to_csv(dst / "GSE111672_download_manifest.tsv", sep="\t", index=False)


def main() -> int:
    read_hwang()
    prepare_moncada()
    summary = {
        "prepared": ["GSE199102", "GSE111672"],
        "raw_root": str(RAW.relative_to(ROOT)),
        "processed_root": str(OUT.relative_to(ROOT)),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "phase9b3b_prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
