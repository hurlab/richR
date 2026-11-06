# PROJECT_HANDOFF.md — richR

Updated: 2026-11-06. Artifacts are English; conversation is Korean.

## Current task — active post-release review

The operator invoked autonomous Hur Harness code review and fix, followed by
richR/richStudio compatibility verification. This is a new task after the
completed 0.2.0 release, not unfinished release work.

- Workspace: `/home/juhur/.moai/worktrees/projects/richR/richr020-review-fix`.
- Branch: `WT-review-fix-020`; source base:
  `5fdadf3cb8c75fc8b2b5682d1c9db8eab9ebaac1` (published 0.2.0).
- Development version: `0.2.0.9000`. The published `v0.2.0` tag/assets and
  released NEWS section remain immutable. README's latest-release pin remains
  `v0.2.0`. No new release, tag, stable promotion or richStudio deployment is
  authorized by this review workflow.
- Only confirmed defects receive minimal fixes. No modernization or engine
  universe-support removal is part of this task.
- Do not modify the shared primary checkout's branch or another worker's files.
  Harness and review evidence stay ignored; never stage them into the package.

## Findings and sequencing

Twelve reproduced findings are grouped into six SPECs. Reviewer R1–R7,
orchestrator P1–P3 and security SEC-001/SEC-002 are findings, not twelve SPECs.

| SPEC | Fix scope | Review snapshot |
|---|---|---|
| RICHR-010 | Literal separators and correct detail-score joins | Independent Hur reviewer APPROVED, commit `12829ba` |
| RICHR-011 | Encoding-equivalent native grouping | Independent Hur reviewer APPROVED, commit `957de6c` |
| RICHR-012 | GSEA opaque IDs and unique comparison group names | Independent Hur reviewer APPROVED, commit `c96524d` |
| RICHR-013 | Sparse-object display and exact-ID/accessor selection | Independent Hur reviewer APPROVED after revision, commit `a42a3ce`; merged |
| RICHR-014 | Empty/sparse plots, multi-result network identity, HTML-save forwarding | Independent Hur reviewer APPROVED, commit `b298d10`; merged |
| RICHR-015 | Exact manifest dimensions and integer metadata | Independent Hur reviewer APPROVED, commit `528e528`; merged |

All six approved fixes, 010 through 015, are merged into the integration
branch. The 013 revision retains the baseline complete-table console layout;
its earlier NEEDS REVISION record remains review history. These per-fix approvals
are not final integrated QA or a binding sync-auditor verdict. The sole
file-conflict dependency is **010 → 014**, both touching
`R/ggnetwork.R`. Cross-SPEC regression checks also cover accessor/plot separators,
GSEA internal-space IDs, native membership fidelity and folder interchange.
Initial plan audit required a correction to 014: preserve eligible connected
source rows without adding isolated vertices or pooling statistics. Focused
follow-up `plan-audit-014-followup.md` records PASS for that correction.

## Evidence

Evidence root: `.moai/reports/review-020/` in this worktree. These are current-run
records, not fresh tests executed by the documentation author.

| Check and retained command | Observed output |
|---|---|
| QA `timeout 900 python3 .../qa/baseline.py` | `FILES: 52 BLOCKS: 377 PASSED: 5964 FAILED: 0 ERRORS: 0 SKIPPED: 3 WARNINGS: 44` |
| QA environment-scrubbed `timeout 1200 Rscript --vanilla baseline-coverage.R` | `LINES: 5479 / 6134 = 89.32181 percent`; `EXPRESSIONS: 3697 / 4225 = 87.50296 percent` |
| QA `Rscript --vanilla baseline-parse.R` | `PARSED_FILES: 104 ERRORS: 0` |
| QA fresh installed core smoke | `CORE_SMOKE: mapping audit, ORA, tested table, interchange, themes, sensitivity, report/rerun, GSEA PASS` |
| Security `secret-scan.py`, tracked-text scope | `text_files_scanned: 301`, `candidate_count: 0` |
| Security `sast.R` | `R AST SAST: files 43 sensitive-call records 60` |
| Security OSV batch plus known-affected jsonlite control | 123 version queries, zero matches; affected jsonlite 1.8.7 returned RSEC-2023-3 |

Baseline-attribution: fresh release archive at `5fdadf3`, tree
`d180faed8bd8d0c6d7b4bacb26d4108153f378bb`, R 4.6.2 Linux,
`/tmp/hur-richr020-qa-baseline/lib/richR`. Full QA/security commands, raw logs,
limits and provenance are in `qa/baseline.md` and `security/baseline.md`.
Source reproductions are in `reviewer.md` and `orchestrator/probe.log`.
Per-fix decisions are in `reviews/SPEC-RICHR-010.md` through `015.md`.
The 014 reviewer matched installed source and independently ran six files with
`Assertions: 652 failed: 0 errors: 0 warnings: 0 skips: 0`. This is scoped
review evidence, not a package-wide integrated test result.

## richStudio compatibility handoff

Use `.moai/reports/review-020/consumer/richStudio-handoff.md` for the reusable
owner-facing handoff, with references to three H1–H7 reports and the independent
background-capability probe. The sibling source was read-only. Its source pin
still references richR 0.1.6; no deployed dependency version is inferred.
The operator will copy the handoff to the richStudio Claude CLI owner. Do not
edit that sibling, change its pin or deploy it from this task.

## Gaps and residual risk

Final integrated installation/suite, global coverage comparison, package check,
post-fix security probes and binding sync audit remain pending. Do not close
SPECs or report this task completed until those records exist. Baseline gaps:
three opt-in provider tests skipped; 44 alias warnings retained; native
Windows/macOS, lintr and hosted CI unobserved. Bounded OSV/secret/AST results do
not certify all dependencies or execution paths. Filesystem interchange remains
nontransactional and is not a filesystem sandbox.

The consumer checks used source-level functions and controlled fixtures, not a
running deployed Shiny session. The original biological demo analysis cannot be
reconstructed from its saved result TSV alone. Legacy result provenance may be
unknown and must not be invented.

## Prior completed release

0.2.0 release completion remains historical. Its authoritative prior handoff is
preserved verbatim at `.moai/reports/review-020/prior-PROJECT_HANDOFF.md`;
the newer historical log is preserved below the new entry in PROJECT_LOG.md.
No release verification or release-asset download was rerun by this documentation
step. Retain all worktrees until parent integration and remote landing are verified.
