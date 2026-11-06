# PROJECT_HANDOFF.md — richR 0.2.1 release preparation

Updated: 2026-11-06. Artifacts are English; conversation is Korean.

## Current task

The owner authorized a new 0.2.1 release, promotion of the stable branch, and an
upstream contribution to guokai8/richR. This supersedes the preceding review's
restriction against release publication; the old 0.2.0 tag and assets remain
immutable.

- Workspace: `/home/juhur/.moai/worktrees/projects/richR/richr-release-021`.
- Branch: `WT-release-021`; source base:
  `fa383355357ee21b21707c1e35abe79520b4b21e`.
- Release metadata is being prepared for version `0.2.1`. GitHub publication,
  assets and branch promotion are pending fresh release verification.
- The repository's existing default branch is `stable`. Keep that designation;
  advance `stable` and `master` to the validated release commit.
- Upstream work is separately scoped to avoid duplicating the existing open
  contributions. No upstream PR publication is claimed by this document.
- Distribution policy: GitHub now; any future package-repository submission
  targets Bioconductor. CRAN submission is not a release gate. No Bioconductor
  submission-readiness assessment is claimed.

## Release scope

The patch packages twelve confirmed fixes across six reviewed areas: literal
gene separators/detail joins, encoding-equivalent native keys, GSEA identifiers
and group names, result access/display, plot boundaries/network exports, and
integer validation during folder import. Package-title and README metadata
corrections are included. The 0.2.0 significance defaults and engine custom
background support remain unchanged.

## Verification attribution

The preceding development runtime at
`7df78b07aabee622625192db4a245734e35f03da` passed 6,369 assertions with zero
failures/errors; line coverage was 5,561/6,183 = 89.9401585%. Its independent
review passed all 26 acceptance criteria. The metadata follow-up at
`30e55792951451ac38bd8c2018b423a4aee0053c` passed standard package check with
0 ERROR / 0 WARNING / 0 NOTE and 175 scoped fresh-built native/ORA assertions.
These are **previous development-baseline results**, not fresh 0.2.1 results.

Authoritative previous reports remain in the prior isolated worktree:

- [Final review report](/home/juhur/.moai/worktrees/projects/richR/richr020-review-fix/.moai/reports/review-020/final-report.md)
- [Prior full runtime QA](/home/juhur/.moai/worktrees/projects/richR/richr020-review-fix/.moai/reports/review-020/qa/final-qa.md)
- [Prior standard package check](/home/juhur/.moai/worktrees/projects/richR/richr020-review-fix/.moai/reports/review-020/qa/standard-check-30e557929514/standard-check.md)

Fresh 0.2.1 QA, asset hashes, publication receipts, remote branch refs and the
upstream PR URL will be recorded under `.moai/reports/release-021/` as observed.
Release-note draft: `.moai/reports/release-021/release-notes.md`.

## Historical review handoff

The snapshot below describes the completed preceding task and its former
workspace. Relative evidence paths in that snapshot refer to the prior review
worktree, not this release worktree.

---

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
| RICHR-014 | Empty/sparse plots, multi-result network identity, HTML-save forwarding | Independent Hur reviewer APPROVED at `b298d10`, then `stats::sd` revision `0ca51d3`; merged |
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

## Integration attempt history

The first freeze `6086335719a0f5465201d0b88d32b7b614a2d416` passed its fresh
installed suite (6,369 passes, 0 failures/errors) and scoped security probes.
Its completed package check nevertheless reported `Status: 1 WARNING, 4 NOTEs`,
including `no visible global function definition for ‘sd’` in ggheatmap. It is
**not the final accepted baseline**. Preserve its QA/security directories as
attempt history; do not attribute a clean code-usage check to that attempt.

The author changed only both `sd(x)` calls to `stats::sd(x)` in revision
`0ca51d33c96d8c1f56bb3a60a7adb08068fb10d5`. Independent review APPROVED the
revision, with 69 assertions and 0 failures/errors/warnings/skips, 0 detached-
stats sd visibility diagnostics and `Package code-usage issue count: 0`.
That scoped usage check is not a replacement for the corrected full package check.
The fresh corrected freeze and its current results follow.

## Corrected frozen integration — final gates still active

Frozen source: `7df78b07aabee622625192db4a245734e35f03da`, Git tree
`4bb38b4c2232b4d5c442e5905b6a438ce97cfc5f`, development version 0.2.0.9000.
The fresh archive installation is
`/tmp/hur-richr020-qa-final-7df78b07aabe/lib/richR`. QA provenance, installed
runner, suite output and raw results are in `qa/final-7df78b07aabe/`.
The retained `run-installed.py` exports the exact archive, scrubs library/test
environment variables, installs privately, asserts the loaded package path and
version, then runs the installed-package suite. `suite.log` records:

```text
FILES: 59 BLOCKS: 411 PASSED: 6369 FAILED: 0 ERRORS: 0 SKIPPED: 3 WARNINGS: 84
```

Warning classification (`classify-warnings.R`, `warning-classification.log`)
mechanically compares baseline and final raw results:

```text
WARNING_TOTALS: baseline 44 final 84 same_baseline_files 44 new_files 40
NEW_MESSAGE_TYPES: 0
EXISTING_FILE_WARNING_SIGNATURES_IDENTICAL: TRUE
BASELINE_ALIAS_GENERIC_WARNINGS: 24
```

The baseline's 44 warnings comprise 24 alias-generic and 20 igraph deprecations;
they are not all alias warnings. Final warnings comprise 24 alias-generic and
60 igraph deprecations. New `test-review-network-exports.R` exercises existing
igraph deprecated calls and contributes the extra 40. Earlier historical release
notes in the session log calling all 44 warnings alias warnings are imprecise;
this current-run classification is the correction, not a rewritten old baseline.

All suite and scoped security results in this section were rerun against the
corrected archive/private installation; they are not carried forward from 608.
The completed scoped security report is
`security/postfix-7df78b07aabe/final-security.md`. Current-run commands include
tracked-text `secret-scan.py`, AST `sast.R`, public manifest/typed-integer/HTML
probes, dependency inventory and live `osv-query.py` with a known-positive control.
Observed: 308 text files, 0 secret candidates; 43 R files/61 sensitive-call records;
123 advisory queries/0 matches, with affected jsonlite 1.8.7 returning RSEC-2023-3.
Dependency versions are unchanged. Public malformed-integer cases are rejected,
valid round-trip control is TRUE, and richNetwork/richNetmap produce nonempty
requested HTML. Runtime probes assert the same fresh final namespace. Final
source recheck reports `SOURCE_HASH_MISMATCHES: []` for 51 recorded files.

No new Critical/High finding was established within that limited security scope.
This is not a vulnerability-free assertion.

Final QA raw evidence in `qa/final-7df78b07aabe/` now records:

| Gate / command evidence | Observed result |
|---|---|
| Full-package `coverage.R` and `validate-results.R` | Lines 5561/6183 = 89.9401585%; baseline 5479/6134 = 89.3218128%; +0.6183457 percentage points. Expressions 3775/4270 = 88.40749%. Instrumented suite matches uninstrumented suite. |
| `parse.R`, `parse.log` | `PARSED_FILES: 111 ERRORS: 0` |
| `joint-smoke.R`, `joint-smoke.log` | 34 regression blocks, 405 passes, 0 failures/errors/skips, 40 known igraph warnings; 23 further successes; explicit-universe independent phyper controls and all six integrated contracts PASS |
| `smoke.R`, `smoke.log` | Eight successes; mapping/ORA/tested/interchange/themes/sensitivity/report-rerun/GSEA core smoke PASS |
| Browser artifact smoke | Exit 0, 14 unique IDs, 4 resolving internal links, 2 embedded PNGs; QA visually inspected the generated report screenshot |
| `run-packagecheck.py`, build/check logs and copied `00check.log` | Build/check driver exit 0; `* DONE`, `Status: 1 WARNING, 3 NOTEs`; R code possible-problems check OK |

The package check explicitly uses `R CMD check --as-cran --no-tests`,
`NOT_CRAN=true`, `_R_CHECK_CRAN_INCOMING_=false`; its test step is SKIPPED.
The separate exact-source installed suite above supplies full-suite evidence.
The warning is missing host `inconsolata.sty`; notes are a compiler flag,
missing HTML Tidy and manual-build detritus. They remain reported diagnostics,
not a warning-free package check. The corrected check has no undefined-sd NOTE.

The final read-only binding sync audit records **PASS, 97.75/100, 26/26 ACs**
at `.moai/reports/review-020/final-audit.md`, against exact source `7df78b0`.
The auditor independently recounted source hashes, line-coverage CSVs and warning
CSVs, and inspected saved independent runtime evidence; its read-only sandbox
could not rerun R. QA's consolidated report is `qa/final-qa.md`.

All six scoped code SPECs are completed after that audit. E4 records exact-source
acceptance and `sync_commit_sha: pending-backfill` until parent commit. E2/E3
remain author history. Plan/acceptance/progress files without frontmatter retain
their existing structure. No SPEC body requirement or acceptance criterion changed.

## Check-profile follow-up — standard check passed

The original 7df package check used `NOT_CRAN=true` and
`_R_CHECK_CRAN_INCOMING_=false` as well as `--as-cran --no-tests`. It did not
observe incoming checks. The canonical `qa/final-qa.md` now states that caveat.
The frozen report copy was restored byte-identical; its separate
`qa/final-7df78b07aabe/report-clarification-receipt.json` preserves the addendum
history rather than rewriting the frozen evidence manifest.

The isolated default-incoming check at the same 7df source completed with
**1 ERROR, 0 WARNINGs, 0 NOTEs**. Compilation flags, real PDF manual, HTML
validation and detritus checks each report OK. It is not a CRAN-ready result.
The structured failure is incoming feasibility: archived `RichR` name collision,
large development-version component and GO.db availability in consulted software
repositories, plus the then-existing URL diagnostics. Official GO.db annotation
page is reachable; no package rename, version change or annotation dependency
removal is authorized by those policy/repository findings.

`qa/clean-check-7df78b07aabe/clean-check.md` records real Inconsolata Type1 font
embedding, scoped tidy 5.6 provisioning and a temporary Makevars removal of only
`-mno-omit-leaf-frame-pointer`, retaining other hardening. No check suppression
or incoming override was used. The check-built library subsequently passed 175
native/ORA assertions. Original collector export failed on a list-valued CSV
column after saving raw results; a separate reader verified those raw counts and
exported scalar columns without rerunning or hiding that exit.

The independently approved metadata-only task was committed as
`30e55792951451ac38bd8c2018b423a4aee0053c`: DESCRIPTION Title case, two README
HTTPS URLs and broken disabled-Actions badge removal. Version 0.2.0.9000,
Package richR, GO.db and the README v0.2.0 release pin remain. Actions remains
disabled; no CI pass is asserted. QA's
`qa/standard-check-30e557929514/runtime-equivalence.json` confirms 305
runtime/native/test/help/data/inst/vignette/NAMESPACE files byte-identical to 7df;
the committed archive differs only in DESCRIPTION and README.

Standard `R CMD check --no-tests` on metadata commit `30e5579` now completed with
**exit 0, Status: OK (0 ERRORs, 0 WARNINGs, 0 NOTEs)**. Its profile uses neither
`--as-cran` nor an incoming override. Evidence in
`qa/standard-check-30e557929514/`: `commands.json`, `outcome.json`,
`provenance.json`, `00check.log` and `postcompile-regression.log`; consolidated
report: `qa/standard-check-30e557929514/standard-check.md`. The fresh
check-built private library passes `TOTAL 175 FAIL 0 ERROR 0 WARN 0 SKIP 0`.
The real PDF manual check is OK. This standard profile does not run a separate
HTML-manual validation step; that step passed in the strict 7df check with
byte-identical Rd content. No HTML step is inferred from the standard result.

**Code and environment review are complete; final commit and remote receipt
remain pending.** The parent will write the actual commit/push observations to
`.moai/reports/review-020/git-receipt.json`; that path is a planned receipt,
not evidence of an already-created commit or remote update. The standard PASS
does not replace the preserved strict incoming ERROR. No CRAN-ready claim,
new release/tag or stable promotion is made.

## richStudio compatibility handoff

Use `.moai/reports/review-020/consumer/richStudio-handoff.md` for the reusable
owner-facing handoff, with references to three H1–H7 reports and the independent
background-capability probe. The sibling source was read-only. Its source pin
still references richR 0.1.6; no deployed dependency version is inferred.
The operator will copy the handoff to the richStudio Claude CLI owner. Do not
edit that sibling, change its pin or deploy it from this task.

## Gaps and residual risk

Binding code audit and isolated standard check passed; final commit and remote
integration receipt remain pending. Strict incoming feasibility remains an observed ERROR. Three opt-in provider tests remain skipped; observed
warning classifications are above. Native
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
