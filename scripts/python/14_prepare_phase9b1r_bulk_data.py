#!/usr/bin/env python3
from pathlib import Path
import sys
import pandas as pd

ROOT = Path("/Users/emily/thesis/PDAC")
PROC = ROOT / "03_processed/external/phase9_bulk"
TABLE = ROOT / "05_results/tables"

EXPECTED = {
    "TCGA_PAAD": 178,
    "GSE71729": 145,
    "GSE62452": 69,
}


def main() -> int:
    rows = []
    for cohort, expected_n in EXPECTED.items():
        expr_path = PROC / cohort / f"{cohort}_expression_gene_by_sample.tsv.gz"
        meta_path = PROC / cohort / f"{cohort}_sample_metadata.tsv"
        if not expr_path.exists() or not meta_path.exists():
            raise SystemExit(f"Missing qualified Phase 9B1 processed input for {cohort}")
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        meta = pd.read_csv(meta_path, sep="\t")
        if expr.shape[1] != expected_n:
            raise SystemExit(f"{cohort} sample count changed: observed {expr.shape[1]}, expected {expected_n}")
        rows.append({
            "cohort": cohort,
            "expression_matrix": str(expr_path.relative_to(ROOT)),
            "sample_metadata": str(meta_path.relative_to(ROOT)),
            "samples": expr.shape[1],
            "genes": expr.shape[0],
            "metadata_rows": meta.shape[0],
            "qualified_phase9a_cohort": True,
            "phase9b1r_action": "reuse_locked_qualified_processed_bulk_input",
        })
    TABLE.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TABLE / "phase9b1r_bulk_input_inventory.tsv", sep="\t", index=False)
    print("Phase 9B1R bulk inputs verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
