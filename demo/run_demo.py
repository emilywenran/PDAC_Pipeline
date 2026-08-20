#!/usr/bin/env python3
"""One-click, self-contained demo of the PurIST PDAC subtype classifier.

WHAT THIS DEMO IS
------------------
It applies the real, locked PurIST classifier coefficients used in this
project (read from reference_data/PDAC_subtype_signatures/PurIST_signatures.tsv)
to a small SYNTHETIC, seeded toy expression matrix generated in this script.

WHAT THIS DEMO IS NOT
----------------------
It is NOT a reproduction of the project's validated results. The real GSE172356
expression matrix is not redistributed in this release (see README.md and
docs/workflows/GITHUB_PUBLICATION_PACKAGE_PLAN.md, Section 5), so no real
patient data is used here. The synthetic samples below are random numbers with
no biological meaning; classifications printed for them are illustrative only.

WHY IT EXISTS
-------------
1. To verify your environment (numpy, pandas) is set up correctly.
2. To show, end-to-end and in a few seconds with no downloads, exactly how the
   locked PurIST scoring formula (intercept + coefficient-weighted gene-pair
   indicators -> logistic transform) works on this project's real coefficients.

To run the actual validated pipeline on real data, see run_pipeline.py and the
README's "Quick Start" / "Suggested Next Steps" sections.

Usage:
    python demo/run_demo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SIGNATURE_FILE = REPO_ROOT / "reference_data" / "PDAC_subtype_signatures" / "PurIST_signatures.tsv"

# Locked intercept, as recorded in docs/planning/DECISION_LOG.md and reproduced
# in analysis_reports/09_external_validation/PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md.
PURIST_INTERCEPT_BETA0 = -6.815

# Standard logistic decision boundary (probability >= 0.5 -> Basal-like).
# This is the textbook default, not independently re-verified against the
# project's internal "locked cutoff" value for this demo — treat classifications
# below as illustrative, not as a substitute for the validated pipeline.
DECISION_THRESHOLD = 0.5

N_SYNTHETIC_SAMPLES = 6
RANDOM_SEED = 42


def load_signature(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"ERROR: signature file not found at {path}", file=sys.stderr)
        sys.exit(1)
    sig = pd.read_csv(path, sep="\t")
    required_cols = {"gene_A", "gene_B", "coefficient", "direction"}
    missing = required_cols - set(sig.columns)
    if missing:
        print(f"ERROR: signature file is missing expected columns: {missing}", file=sys.stderr)
        sys.exit(1)
    return sig


def make_synthetic_expression(sig: pd.DataFrame, n_samples: int, seed: int) -> pd.DataFrame:
    """Generate a small synthetic, seeded log-expression matrix covering exactly
    the 16 genes referenced by the PurIST signature (8 pairs x 2 genes each).
    Values are draws from a log-normal distribution with no biological meaning.
    """
    rng = np.random.default_rng(seed)
    genes = pd.unique(pd.concat([sig["gene_A"], sig["gene_B"]]))
    samples = [f"synthetic_sample_{i + 1}" for i in range(n_samples)]
    values = rng.lognormal(mean=2.0, sigma=1.0, size=(len(genes), len(samples)))
    return pd.DataFrame(values, index=genes, columns=samples)


def classify_purist(sig: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    """Apply the locked PurIST formula:
        linear_predictor = beta0 + sum_i( coefficient_i * 1[expr(gene_A_i) > expr(gene_B_i)] )
        p_basal = 1 / (1 + exp(-linear_predictor))
    for every gene pair with direction 'gene_A > gene_B' (the only direction
    present in this project's locked signature file).
    """
    results = []
    for sample in expr.columns:
        linear_predictor = PURIST_INTERCEPT_BETA0
        for _, row in sig.iterrows():
            a_val = expr.loc[row["gene_A"], sample]
            b_val = expr.loc[row["gene_B"], sample]
            indicator = 1.0 if a_val > b_val else 0.0
            linear_predictor += row["coefficient"] * indicator
        p_basal = 1.0 / (1.0 + np.exp(-linear_predictor))
        label = "Basal-like" if p_basal >= DECISION_THRESHOLD else "Classical"
        results.append({"sample": sample, "p_basal_like": round(p_basal, 4), "predicted_subtype": label})
    return pd.DataFrame(results)


def main() -> None:
    print("=" * 78)
    print("PDAC PurIST classifier demo — SYNTHETIC DATA, NOT real patient data")
    print("=" * 78)
    print(f"Loading locked PurIST signature from: {SIGNATURE_FILE.relative_to(REPO_ROOT)}")
    sig = load_signature(SIGNATURE_FILE)
    print(f"Loaded {len(sig)} gene pairs. Locked intercept (beta0) = {PURIST_INTERCEPT_BETA0}")

    print(f"\nGenerating {N_SYNTHETIC_SAMPLES} synthetic samples (seed={RANDOM_SEED})...")
    expr = make_synthetic_expression(sig, N_SYNTHETIC_SAMPLES, RANDOM_SEED)
    print("Synthetic log-expression matrix (first 5 genes):")
    print(expr.head(5).round(2).to_string())

    print("\nApplying the locked PurIST scoring formula...")
    results = classify_purist(sig, expr)
    print("\nResults (illustrative only — synthetic input data):")
    print(results.to_string(index=False))

    print("\n" + "=" * 78)
    print("Demo complete. Your environment can load the project's locked signature")
    print("files and execute the classifier formula end-to-end.")
    print("For real analyses, see run_pipeline.py and README.md.")
    print("=" * 78)


if __name__ == "__main__":
    main()
