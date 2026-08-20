#!/usr/bin/env python3
import pandas as pd
import hashlib
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
MANIFEST_PATH = ROOT / "01_metadata/file_manifest.tsv"

files_to_hash = [
    {
        "file_id": "09_docs__planning__DECISION_LOG_md",
        "local_path": "/Users/emily/thesis/PDAC/09_docs/planning/DECISION_LOG.md"
    },
    {
        "file_id": "00_admin__PROJECT_STATUS_md",
        "local_path": "/Users/emily/thesis/PDAC/00_admin/PROJECT_STATUS.md"
    },
    {
        "file_id": "00_admin__SKILL_USAGE_LOG_tsv",
        "local_path": "/Users/emily/thesis/PDAC/00_admin/SKILL_USAGE_LOG.tsv"
    }
]

def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()

def main():
    df = pd.read_csv(MANIFEST_PATH, sep="\t")
    
    for item in files_to_hash:
        path = Path(item["local_path"])
        if not path.exists():
            print(f"Warning: {path} not found")
            continue
            
        size = path.stat().st_size
        sha = compute_sha256(path)
        
        # Update size and md5 columns
        idx = df[df["file_id"] == item["file_id"]].index
        if len(idx) > 0:
            df.loc[idx, "file_size"] = size
            df.loc[idx, "md5"] = sha
            print(f"Updated {item['file_id']}: size={size}, sha256={sha}")
        else:
            print(f"Warning: file_id {item['file_id']} not found in manifest")
            
    df.to_csv(MANIFEST_PATH, sep="\t", index=False)
    print("Saved manifest updates.")

if __name__ == "__main__":
    main()
