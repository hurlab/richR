# PROJECT_HANDOFF.md — richR

## 1. Project Overview

richR is an R package for functional enrichment analysis (GO, KEGG, Reactome, MSigDB, GSEA) with visualization. Originally developed by Kai Guo (`guokai8/richR`), now maintained at `hurlab/richR`.

- **Last updated:** 2026-04-16 CDT
- **Fork version:** 0.1.4 (upstream @ `fd96d4f`, with PR #9 richUpset fix cherry-picked)
- **Last coding CLI used:** Claude Code CLI (Claude Opus 4.6)

## 2. Current State

| Feature / Milestone | Status | Notes |
|---|---|---|
| v0.1.0 — Codebase review, bug fixes, documentation (PR #6) | MERGED upstream 2026-04-11 | Squashed into upstream commit `086f655` |
| v0.1.2 — show methods, validation, GMT import, tests (PR #8) | MERGED upstream 2026-04-11 | Squashed into upstream commit `086f655` |
| v0.1.1 — Color-enhanced UpSet plot (PR #7) | Superseded by PR #9 | PR #7 proposed a new `ggupset()` function; PR #9 instead mainlines the fix directly into `richUpset()` |
| PR #9 — richUpset alignment fix + multi-set color blending | OPEN, pending review | https://github.com/guokai8/richR/pull/9, targets `master` |
| Upstream sync → `fd96d4f` | Completed 2026-04-15 | Branch `sync/upstream-v0.1.2` created as review staging |
| Fork release 0.1.4 | Completed 2026-04-16 | `hurlab/0.1.4` branch: upstream + PR #9 patch + fork README/figures/docs |
| ggUpSet retired | Completed 2026-04-16 | No separate `ggUpSet()` or `ggupset()` function; improvements are in `richUpset()` directly; `ggupset` remains as a simple alias from upstream's `zzz_ggplots_aliases.R` |
| README: 4 figures embedded in Visualization sections | Completed 2026-04-16 | Bar, Dot, GSEA, Compare-dot, UpSet figures restored from safety tag |
| README: citation updated to hurlab/richR + 2026 + v0.1.4 | Completed 2026-04-16 | DOI removed pending hurlab-owned Zenodo release |
| R CMD check clean build | Not started | Blocked by missing Bioconductor + CRAN dependencies |
| Integration test coverage | Not started | Need tests for core enrichment pipeline |

## 3. Execution Plan Status

| Phase | Status | Last Updated | Notes |
|---|---|---|---|
| Phase 0: Upstream sync | Completed | 2026-04-15 | Safety tag `pre-upstream-sync-backup`; sync branch on origin |
| Phase 0.5: Upstream PR for richUpset | Completed (awaiting merge) | 2026-04-16 | PR #9 opened |
| Phase 0.6: Hurlab fork 0.1.4 release | Completed (awaiting master promotion) | 2026-04-16 | `hurlab/0.1.4` branch built |
| Phase 1: Environment & Build Health | Not started | 2026-04-16 | Bioconductor (`GO.db`, `fgsea`, `S4Vectors`, `DOSE`) + CRAN (`igraph`, `visNetwork`, `reshape2`, `intergraph`, `ggrepel`, `GGally`, `sna`, `testthat`, `devtools`) still missing |
| Phase 2: Documentation consistency | Partially complete | 2026-04-16 | CLAUDE.md still references removed `exportPattern`; man pages not regenerated |
| Phase 3: Test coverage expansion | Not started | 2026-04-16 | Upstream has test-helpers/validation/s4-classes from PR #8; integration tests still missing |
| Phase 4: Additional figures | Not started | 2026-04-16 | Only 4 example figures available; volcano/scatter/termSim/network/etc. need regeneration once deps are installed |

## 4. Outstanding Work

| Item | Status | Last Updated |
|---|---|---|
| Promote `hurlab/0.1.4` → `master` | Pending user approval | 2026-04-16 |
| Force-push `master` to origin | Pending | 2026-04-16 |
| Delete stale merged `feature/*` branches on origin | Pending | 2026-04-16 |
| Monitor PR #9 for Kai's review | External (blocked) | 2026-04-16 |
| Install missing R dependencies | Not started | 2026-04-16 |
| Run `R CMD check` | Not started | 2026-04-16 |
| Update CLAUDE.md (remove stale `exportPattern` reference) | Not started | 2026-04-16 |
| Regenerate additional example figures (volcano, scatter, etc.) | Not started | 2026-04-16 |
| Write integration tests | Not started | 2026-04-16 |
| Close PR #7 once PR #9 merges (optional) | Pending | 2026-04-16 |

## 5. Risks, Open Questions, and Assumptions

| Item | Status | Notes |
|---|---|---|
| `0.1.3-hurlab.1` would fail R CMD check (non-numeric `hurlab` token) | Resolved | Using plain `0.1.4` — fork leads upstream's 0.1.3 release line |
| Missing deps prevent R CMD check | Open | Upstream added `sna` + dropped `tidyr`/`methods` from Imports |
| Safety tag `pre-upstream-sync-backup` retained | Active | Local-only; points to pre-sync master tip |
| Target audience (CRAN vs GitHub-only) | Open (default: GitHub-only) | |

## 6. Verification Status

| Item | Method | Result | Date |
|---|---|---|---|
| PR #6 MERGED upstream | `gh pr view` | mergedAt 2026-04-11 | 2026-04-15 |
| PR #8 MERGED upstream | `gh pr view` | mergedAt 2026-04-11 | 2026-04-15 |
| PR #9 OPEN | `gh pr create` returned pull/9 URL | State OPEN | 2026-04-16 |
| `R/ggupset.R` parses cleanly with richUpset fix + `.blend_colors` | `Rscript -e parse()` | 3 top-level expressions OK | 2026-04-16 |
| Figures restored from safety tag | `git checkout pre-upstream-sync-backup -- man/figures/*` | 5 PNGs present in man/figures/ | 2026-04-16 |
| R CMD check | Not verified | — | Blocked by missing dependencies |

## 7. Restart Instructions

- **Starting point:** branch `hurlab/0.1.4` at upstream `fd96d4f` + fork commits.
- **Safety tag:** `pre-upstream-sync-backup` for full rollback.
- **Open upstream PR:** https://github.com/guokai8/richR/pull/9
- **Recommended next actions:**
  1. If user approves, promote: `git checkout master && git reset --hard hurlab/0.1.4 && git push origin master --force-with-lease`
  2. Delete stale branches on origin: `feature/codebase-review-docs-v0.1.0`, `feature/enhancements-v0.1.2`
  3. Keep `feature/upset-plot-v0.1.1` until PR #7 is closed (superseded by PR #9)
  4. Install deps + run R CMD check + fix any surfaced issues
  5. Monitor PR #9 for Kai's review; on merge, bump hurlab version to 0.1.3 (drop the `.1` patch segment) and sync
- **Key branches:**
  - `hurlab/0.1.4` — staged fork release, built on upstream + PR #9 + fork-specific additions
  - `upstream-pr/richupset-align-and-blend` — clean PR #9 branch (upstream + single commit)
  - `sync/upstream-v0.1.2` — earlier sync-review staging (retained until master is promoted)
  - `master` — unchanged from pre-sync state until promotion
- **Last updated:** 2026-04-16 CDT
