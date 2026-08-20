#!/usr/bin/env python3
"""Validate Phase 11B manuscript draft guardrails and claim traceability."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT = ROOT / "04_analysis/11_manuscript/PHASE11B_MANUSCRIPT_DRAFT.md"
TRACE = ROOT / "05_results/tables/phase11b_claim_to_text_trace.tsv"
CLAIM_MAP = ROOT / "05_results/tables/phase11a_claim_evidence_map.tsv"
PROHIBITED = ROOT / "05_results/tables/phase11a_prohibited_claims.tsv"

REQUIRED_TRACE_COLUMNS = {
    "claim_id",
    "manuscript_section",
    "exact_claim_summary",
    "supporting_phase",
    "source_report_or_table",
    "evidence_category",
    "allowed_wording",
    "prohibited_overstatement",
    "current_compliance_status",
}

REQUIRED_CLAIMS = {f"CLAIM_{idx:02d}" for idx in range(1, 17)}


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"required file missing: {path}")
    return path.read_text(encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required file missing: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"empty TSV or missing header: {path}")
        missing = REQUIRED_TRACE_COLUMNS - set(reader.fieldnames) if path == TRACE else set()
        if missing:
            raise ValueError(f"{path} missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip().lower()


def has_near(text: str, anchor_pattern: str, required_pattern: str, window: int = 240) -> bool:
    for match in re.finditer(anchor_pattern, text, flags=re.IGNORECASE):
        start = max(0, match.start() - window)
        end = min(len(text), match.end() + window)
        if re.search(required_pattern, text[start:end], flags=re.IGNORECASE | re.DOTALL):
            return True
    return False


def sentence_windows(text: str, anchor_pattern: str, radius: int = 1) -> list[str]:
    pieces = re.split(r"(?<=[.!?])\s+", text)
    windows: list[str] = []
    for idx, sentence in enumerate(pieces):
        if re.search(anchor_pattern, sentence, flags=re.IGNORECASE):
            lo = max(0, idx - radius)
            hi = min(len(pieces), idx + radius + 1)
            windows.append(" ".join(pieces[lo:hi]))
    return windows


def check_absent_near(
    errors: list[str],
    text: str,
    anchor_pattern: str,
    forbidden_pattern: str,
    label: str,
    allowed_negation_pattern: str | None = None,
    radius: int = 1,
) -> None:
    for window in sentence_windows(text, anchor_pattern, radius=radius):
        if re.search(forbidden_pattern, window, flags=re.IGNORECASE):
            if allowed_negation_pattern and re.search(allowed_negation_pattern, window, flags=re.IGNORECASE):
                continue
            errors.append(f"{label}: {window[:320]}")


def validate_trace_table(errors: list[str]) -> None:
    rows = read_tsv(TRACE)
    claim_rows = {row["claim_id"]: row for row in rows}
    missing_claims = sorted(REQUIRED_CLAIMS - set(claim_rows))
    extra_claims = sorted(set(claim_rows) - REQUIRED_CLAIMS)
    if missing_claims:
        errors.append(f"Trace table missing claim IDs: {', '.join(missing_claims)}")
    if extra_claims:
        errors.append(f"Trace table has unexpected claim IDs: {', '.join(extra_claims)}")

    valid_categories = {"Supported", "Partially Supported", "Unsupported", "Exploratory"}
    for row in rows:
        claim_id = row.get("claim_id", "")
        if row.get("evidence_category", "") not in valid_categories:
            errors.append(f"{claim_id}: invalid evidence category {row.get('evidence_category')!r}")
        if row.get("current_compliance_status") != "COMPLIANT":
            errors.append(f"{claim_id}: current_compliance_status must be COMPLIANT")
        for field in REQUIRED_TRACE_COLUMNS:
            if not row.get(field, "").strip():
                errors.append(f"{claim_id}: empty required trace field {field}")

    unsupported = {"CLAIM_09", "CLAIM_10", "CLAIM_11", "CLAIM_12", "CLAIM_13", "CLAIM_16"}
    for claim_id in unsupported:
        if claim_rows.get(claim_id, {}).get("evidence_category") != "Unsupported":
            errors.append(f"{claim_id}: unsupported claim must remain evidence_category=Unsupported")


def validate_required_disclosures(errors: list[str], manuscript: str) -> None:
    required_patterns = {
        "HALLMARK_PROTEIN_SECRETION explicitly marked PARTIAL_SPATIAL_SUPPORT": r"HALLMARK_PROTEIN_SECRETION.*?PARTIAL_SPATIAL_SUPPORT|PARTIAL_SPATIAL_SUPPORT.*?HALLMARK_PROTEIN_SECRETION",
        "HALLMARK_SPERMATOGENESIS failure disclosed": r"HALLMARK_SPERMATOGENESIS[^.]{0,220}(opposite directional effects|failed|insufficient spatial coverage|<80%)",
        "MEred and MEpurple replication failures disclosed": r"MEred\s+and\s+MEpurple[^.]{0,160}failed to replicate",
        "9 TF activity replication failures disclosed": r"9\s+(distinct\s+)?transcription factor activities[^.]{0,180}failed to replicate|9\s+activities[^.]{0,120}failed to replicate",
        "rCLR sensitivity disclosed": r"rCLR|robust centred log-ratio",
        "Herbaspirillum contamination limitation disclosed": r"Herbaspirillum[^.]{0,180}(contaminant|contamination)",
        "Moncada exploratory 1 of 6 disclosed": r"Moncada[^.]{0,160}exploratory[^.]{0,160}1\s+of\s+6|exploratory[^.]{0,160}Moncada[^.]{0,160}1\s+of\s+6",
        "CTCFL/BORIS exclusion disclosed": r"CTCFL.*?(BORIS).*?(excluded|removed|not prioritized|penalized)|CTCFL.*?(excluded|removed|not prioritized|penalized).*?(BORIS)",
    }
    for label, pattern in required_patterns.items():
        if not re.search(pattern, manuscript, flags=re.IGNORECASE | re.DOTALL):
            errors.append(f"Missing required disclosure: {label}")

    if not has_near(manuscript, r"HALLMARK_PROTEIN_SECRETION", r"PARTIAL_SPATIAL_SUPPORT", window=900):
        errors.append("HALLMARK_PROTEIN_SECRETION is not explicitly tied to PARTIAL_SPATIAL_SUPPORT nearby.")

    if not has_near(manuscript, r"Moncada", r"positive concordance in only 1 of 6 sections", window=320):
        errors.append("Moncada must be described as exploratory with positive concordance in only 1 of 6 sections.")


def validate_prohibited_language(errors: list[str], manuscript: str) -> None:
    check_absent_near(
        errors,
        manuscript,
        r"CTCFL|BORIS",
        r"\b(prioritized|prioritised|intermediate[- ]tier|therapeutic|candidate[- ]supported)\b",
        "CTCFL/BORIS prohibited support language",
        allowed_negation_pattern=r"\b(not|no|excluded|removed|penalized|penalised|blocked)\b.{0,80}\b(prioritized|prioritised|candidate|tier|therapeutic)\b|\b(prioritized|prioritised)\b.{0,80}\b(excluded|removed|not)\b",
    )
    check_absent_near(
        errors,
        manuscript,
        r"CTCFL|BORIS",
        r"cancer[- ]testis antigen|literature (support|rescue)|post[- ]hoc",
        "CTCFL/BORIS literature rescue language",
    )
    check_absent_near(
        errors,
        manuscript,
        r"microb|bacter|Ochrobactrum|Herbaspirillum|genera",
        r"\b(causes?|caused|causing|drives?|drove|induces?|mediates?|causality|causal relationship)\b",
        "Microbial causality language",
        allowed_negation_pattern=r"\b(no|not|cannot|capable of establishing|not established|not claimed|no inference|no causal relationship|directionality.*cannot)\b",
        radius=0,
    )
    check_absent_near(
        errors,
        manuscript,
        r"microb|bacter|Ochrobactrum|Herbaspirillum",
        r"\b(locali[sz](ed|ation)|spatial tissue localisation|preferentially in malignant|physical interaction|proximity)\b",
        "Microbial localization or physical interaction language",
        allowed_negation_pattern=r"\b(no|not|cannot|do not|was not determined|do not permit|not directly demonstrated|not claimed|not proposed|not supported|no inference)\b",
        radius=0,
    )
    for window in sentence_windows(manuscript, r"full spatial validation|complete spatial validation|definitive spatial validation", radius=0):
        if not re.search(r"\b(not|no|did not|does not|without|insufficient|partial)\b", window, flags=re.IGNORECASE):
            errors.append(f"Full basal-classical spatial validation language detected: {window[:320]}")


def validate_traceable_evidence(errors: list[str], manuscript: str) -> None:
    trace_rows = read_tsv(TRACE)
    claim_rows = read_tsv(CLAIM_MAP)
    read_tsv(PROHIBITED)

    trace_blob = normalize("\n".join("\t".join(row.values()) for row in trace_rows))
    claim_blob = normalize("\n".join("\t".join(row.values()) for row in claim_rows))

    major_terms = {
        "Ochrobactrum": "CLAIM_03",
        "HALLMARK_PROTEIN_SECRETION": "CLAIM_05",
        "PARTIAL_SPATIAL_SUPPORT": "CLAIM_06",
        "HALLMARK_SPERMATOGENESIS": "CLAIM_09",
        "MEred": "CLAIM_10",
        "MEpurple": "CLAIM_10",
        "Moncada": "CLAIM_14",
        "CTCFL": "CLAIM_16",
        "BORIS": "CLAIM_16",
        "Herbaspirillum": "CLAIM_02",
        "rCLR": "CLAIM_02",
    }
    for term, claim_id in major_terms.items():
        if term.lower() in manuscript.lower() and claim_id.lower() not in trace_blob:
            errors.append(f"Manuscript term {term!r} lacks trace table mapping to {claim_id}.")

    for row in trace_rows:
        summary_tokens = [token for token in re.findall(r"[A-Za-z0-9_/-]{5,}", row["exact_claim_summary"])[:3]]
        if summary_tokens and not any(token.lower() in claim_blob for token in summary_tokens):
            errors.append(f"{row['claim_id']}: trace summary does not appear connected to Phase 11A claim evidence map.")

    unsupported_promotions = [
        r"HALLMARK_SPERMATOGENESIS[^.]{0,120}\b(supported|validated|prioritized|candidate)\b",
        r"\b(MEred|MEpurple)\b[^.]{0,120}\b(supported|validated|prioritized|candidate)\b",
        r"\b(GFI1B|STAT1|ZBTB11|ZNF639|TWIST1|FOXK2|KDM5B|MAFF|TEAD4)\b[^.]{0,120}\b(supported|validated|prioritized|candidate)\b",
    ]
    for pattern in unsupported_promotions:
        if re.search(pattern, manuscript, flags=re.IGNORECASE):
            errors.append(f"Unsupported claim appears promoted without compliant evidence: {pattern}")


def main() -> int:
    errors: list[str] = []
    try:
        manuscript = read_text(MANUSCRIPT)
        validate_trace_table(errors)
        validate_required_disclosures(errors, manuscript)
        validate_prohibited_language(errors, manuscript)
        validate_traceable_evidence(errors, manuscript)
    except Exception as exc:
        errors.append(str(exc))

    if errors:
        print("PHASE 11B MANUSCRIPT VALIDATION FAILED", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print("PHASE 11B MANUSCRIPT VALIDATION PASSED")
    print("Final readiness decision: READY_FOR_PHASE11C_MANUSCRIPT_INDEPENDENT_REVIEW")
    return 0


if __name__ == "__main__":
    sys.exit(main())
