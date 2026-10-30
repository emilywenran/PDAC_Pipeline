# Workflow Archive Documentation Update

Date: 2026-10-30.

## 📦 Scope

This update adds the supplied recovered prompt archive and workflow index to the existing curated public repository. It refreshes the README and reconciles the workflow documentation with the newly available records. Analysis scripts, numerical outputs, manuscript files, software citation version and scientific conclusions are not changed by this update. Later local manuscript revisions and analysis directories have not been bulk-synchronised into this documentation update; the existing public-analysis snapshot retains its original coverage.

The local `github_release` directory is an independent Git repository configured with `origin = https://github.com/emilywenran/PDAC_Pipeline.git`. The local baseline and the remote `main` both resolved to `6f62e71ec8c8749f8b4bc0e9a2545ff6a3a6d994` during the preparation check. The enclosing analysis workspace has a separate Git history.

## 📚 Archive provenance

- [PHASE_PROMPTS.txt](../../PHASE_PROMPTS.txt): byte-for-byte copy of the supplied UTF-8 text, 418,564 bytes and 12,299 LF-delimited lines.
- SHA-256: `555d6f417989233f3b7064612d4c3efaddfb486630af3b9522b59e4428cf9dd8`.
- [WORKFLOW_INDEX.md](../../WORKFLOW_INDEX.md): 47 recovered instruction entries, from Phase 1A through Phase 11A. The publication copy aligns the source filename with `PHASE_PROMPTS.txt` and adds direct line links using GitHub LF-based line numbering. The supplied index counted 72 Unicode line separators as additional lines; the supplied source index is unchanged in the analysis workspace.
- Archive order, prompt wording and the index's execution-evidence statuses are preserved. Missing instructions are identified rather than reconstructed. Historical local development paths are deliberately retained and disclosed.

## ✅ Preparation checks

The release copy's hash and line count were checked against the supplied archive. All 47 indexed line ranges were checked for continuity and bounds. Relative file links in the current README, public workflow index and updated workflow documents were checked. The existing phase-list and describe commands and synthetic PurIST demo were smoke-checked without running the research pipeline. `git diff --check` was checked for whitespace errors.

A targeted scan of the newly supplied prompt text found no conventional private-key, OpenAI-like, GitHub-token, AWS-access-key or Google-API-key patterns. This is a limited pattern scan, not a guarantee that arbitrary sensitive content can be detected. Original local filesystem paths remain visible as documented above.

## 📝 Publication state

This file documents the prepared local update. It is not a record of a completed commit, GitHub push, software release or new analysis-validation run. The existing GitHub repository is the destination for the updated files; a successful push should be checked separately.
