# PROJECT_LOG.md — richR

## Session 2026-10-17 CDT — Fork sync and harness cleanup reconciled

- Primary development line confirmed as `hurlab/richR` remote `master`, version 0.1.4.
- Prior Claude Code session merged upstream `590efdb` and `ef0fb86` into fork
  `19f38fb`, and removed MoAI commit `23f977c` from upstream PR #11.
  This session verified remote state instead of repeating those mutations.
- Created isolated `.moai/worktrees/fork-wrapup`, branch `WT-fork-wrapup`,
  from remote `origin/master` at `19f38fb`. Primary checkout stayed on
  `upstream-pr/richr-v0.1.3` at `eb4c9da`.
- Resolved the previously absent `origin/HEAD` from the remote's advertised
  default branch (`master`); baseline divergence check returned `0 0`.
- Added `/.moai-backups/` to local `.git/info/exclude`. Primary checkout
  `git status --short` then returned no paths. Harness files remain local.
- Both `R/misc.R` and `R/species.R` idconvert definitions returned TP53 and BRCA1
  for Entrez IDs 7157 and 672 with `AnnotationDbi attached=FALSE`.
  Used R 4.6.1 with its usual libraries first and existing R 4.5 annotation
  packages as fallback; no dependency installation or full local package check.
- Inspected CI run 39079967110 for baseline `19f38fb`: Ubuntu, macOS, Windows
  package-check jobs succeeded, each with `Status: 4 NOTEs`.
- Replaced stale master-promotion and PR #9 restart instructions in the handoff.
  PR #11 is OPEN at `eb4c9da`, 3 commits, 24 files, with no harness paths in
  the inspected families; PR #7/#9 are CLOSED and unmerged.
- No production code changed. No upstream reminder was sent and no PR was modified.
- Evidence commands, runtime-check expression, gaps, and residual CI NOTES are
  recorded in PROJECT_HANDOFF.md. Older log entries remain historical.
- Remaining: upstream action on PR #11; optional future dependency/environment
  work and NOTE cleanup require their own scope. The user authorized committing
  the two wrap-up documents and pushing them to hurlab master.

---

## Session 2026-04-16 CDT — PR #9 opened, ggUpSet retired, hurlab 0.1.4 staged

- **Coding CLI used:** Claude Code CLI (Claude Opus 4.6, 1M context)

### Phase(s) worked on

- Refreshed `gh` auth with `workflow` scope via `gh auth refresh -s workflow -h github.com`
  (stored PAT credential helper lacked the scope needed to push workflow files).
- Built a clean, PR-scoped branch `upstream-pr/richupset-align-and-blend`
  that mainlines the v0.1.2.1 UpSet improvements directly into
  `richUpset()` (replacing the previous plan to keep them as a
  separate `ggUpSet()` function on the fork). Opened PR #9 on
  `guokai8/richR`.
- Built the hurlab fork's 0.1.4 release as branch
  `hurlab/0.1.4`: upstream + PR #9 commit + README
  restructure + 4 restored example figures + fork bookkeeping.

### Concrete changes implemented

1. **Upstream PR branch (`upstream-pr/richupset-align-and-blend`)**:
   - `R/ggupset.R` — inlined the bar-color blending and the
     `align_plots + ggdraw/draw_plot` panel assembly into
     `richUpset()`. Added internal `.blend_colors()` helper.
     Updated roxygen `@param main.bar.color` to describe the new
     default behavior. Added `@importFrom cowplot align_plots ggdraw
     draw_plot` and `@importFrom grDevices col2rgb rgb`.
   - `NAMESPACE` — added `importFrom(cowplot, align_plots|ggdraw|draw_plot)`
     and `importFrom(grDevices, col2rgb|rgb)` at the correct
     alphabetical positions.
   - `DESCRIPTION` — version `0.1.2` → `0.1.3`.
   - `UPDATES.md` — new neutral v0.1.3 entry describing the bug fix
     and enhancement without fork-specific language.
   - Opened [PR #9](https://github.com/guokai8/richR/pull/9) on
     `guokai8/richR`.

2. **Hurlab fork branch (`hurlab/0.1.4`)**:
   - Cherry-picked the PR #9 commit.
   - Version bumped to `0.1.4` in `DESCRIPTION`. (The literal
     `0.1.3-hurlab.1` string the user originally proposed is invalid
     per R's version parser, which requires numeric components; after
     briefly staging `0.1.3.1`, the user opted for a clean `0.1.4`
     bump that signals the hurlab fork leads upstream's 0.1.3 release
     line until PR #9 merges.)
   - Restored 5 example figures from `pre-upstream-sync-backup` tag
     into `man/figures/`: `example_barplot.png`,
     `example_comparedot.png`, `example_dotplot.png`,
     `example_gsea.png`, `example_upset.png`.
   - README rewrite:
     - Version badge `0.1.2` → `0.1.4`, badge URL `guokai8` → `hurlab`
     - Install URL `guokai8/richR` → `hurlab/richR`
     - Key Features: new bullet for publication-quality UpSet plots
     - Visualization sections: embedded 4 figures (bar+dot side-by-side,
       GSEA, compare-dot, UpSet)
     - UpSet section: rewritten around `richUpset()` with dedicated
       feature bullets (alignment, blended colors, flexible input);
       4-set example with `set.seed(123)`
     - Citation block: year 2020 → 2026, version 0.1.2 → 0.1.4,
       URL `guokai8` → `hurlab`, DOI line removed pending hurlab
       Zenodo release
     - Contact: points to hurlab/richR issues (upstream still referenced
       for upstream-affecting bugs)
   - `UPDATES.md`: new `0.1.4` fork-packaging entry above the
     shared `0.1.3` entry.
   - `.gitignore`: added tooling dirs (`.claude/`, `.moai/`, `.agency/`,
     `.mcp.json`).
   - `.Rbuildignore`: added `PROJECT_HANDOFF.md`, `PROJECT_LOG.md`,
     tooling dirs, safety tag.
   - `PROJECT_HANDOFF.md`, `PROJECT_LOG.md`: rewritten for this
     session (this file).
   - Did NOT create `R/ggUpSet_hurlab.R` (retired concept).
   - Did NOT edit `R/zzz_ggplots_aliases.R` (upstream's
     `ggupset <- richUpset` simple alias remains intact — users
     calling `ggupset()` get `richUpset()` with the new behavior).

### Files/modules/functions touched

On `upstream-pr/richupset-align-and-blend` (for PR #9):
- `R/ggupset.R` — `richUpset()`, `.blend_colors()` (new), roxygen `@importFrom`
- `NAMESPACE` — 5 new importFrom entries
- `DESCRIPTION` — version
- `UPDATES.md` — v0.1.3 entry

On `hurlab/0.1.4` (additional over PR #9):
- `DESCRIPTION` — version `0.1.3` → `0.1.4`
- `README.md` — badges, install URL, features bullet, figures embed,
  UpSet rewrite, citation, contact
- `UPDATES.md` — v0.1.4 fork entry
- `.gitignore`, `.Rbuildignore` — tooling ignores
- `PROJECT_HANDOFF.md`, `PROJECT_LOG.md` — updated
- `man/figures/*.png` — 5 restored

### Key technical decisions

- **Mainline into `richUpset()` vs maintain `ggUpSet()`**: Earlier
  session preserved the fork version as `ggUpSet()`; user requested
  the simpler path of proposing changes upstream and retiring
  `ggUpSet`. Result: one function, two patches (alignment + blend),
  no fork-specific API surface.
- **Version `0.1.4` (fork leads)**: R's version parser requires
  numeric-only components, ruling out the originally proposed
  `0.1.3-hurlab.1`. After briefly considering `0.1.3.1`, the user
  chose a clean `0.1.4` bump. The hurlab fork therefore leads
  upstream's 0.1.3 release line; PR #9 bumps upstream to 0.1.3 when
  merged, and the fork stays ahead.
- **Minimal PR scope (R + DESCRIPTION + NAMESPACE + UPDATES only)**:
  No README changes, no fork citation, no PROJECT docs. Maximizes
  merge likelihood.
- **Neutral language in upstream UPDATES.md entry**: Per user's
  Q4 preference.
- **`ggupset` remains as alias via upstream's `R/zzz_ggplots_aliases.R`**:
  Since `richUpset()` now has the fixed behavior, the upstream alias
  `ggupset <- richUpset` continues to work — users calling either
  get the same improved function.

### Verification performed

- `Rscript -e parse()` on `R/ggupset.R` — 3 expressions OK
- `git diff --stat upstream/master` on PR branch — 4 files, +82/-17
- `gh pr create` returned `https://github.com/guokai8/richR/pull/9`
- `grep guokai8|ggupset|hurlab README.md` — all references
  intentional and correctly placed
- 5 figures present in `man/figures/`

### Git state at end of session

- Branch `upstream-pr/richupset-align-and-blend` at commit `9d46437`, pushed
- Branch `hurlab/0.1.4` at commit on top of upstream + PR #9,
  NOT yet pushed (pending next commit + master promotion)
- Branch `sync/upstream-v0.1.2` retained for now (can be deleted after
  master promotion)
- Branch `master` unchanged from pre-sync state
- Tag `pre-upstream-sync-backup` retained
- PR #9 open on `guokai8/richR`

### Items NOT yet performed (next actions)

- Commit README/UPDATES/figures/etc. on `hurlab/0.1.4`
- Promote master to hurlab branch + force-push
- Delete stale branches on origin
- Install R dependencies and run `R CMD check`
- Monitor PR #9 for Kai's review

---

## Session 2026-04-15 CDT — Upstream sync after PR #6/#8 merges

- Created safety tag `pre-upstream-sync-backup` on old master.
- Created branch `sync/upstream-v0.1.2` at upstream `fd96d4f`.
- First draft of fork preservation: added `R/ggupset_hurlab.R` (later
  retired this session) and removed upstream's `ggupset <- richUpset`
  alias (also undone this session).
- Established pattern of non-destructive branch-based review before
  promotion.

---

## Session 2026-03-28 16:20 CDT

- **Coding CLI used:** Claude Code CLI (Claude Opus 4.6, 1M context)

### Phase(s) worked on
- UpSet plot alignment and color fix
- README consistency
- PR updates to upstream (guokai8/richR)

### Concrete changes implemented

1. **UpSet bar-to-dot alignment fix** (`R/ggupset.R`):
   - Root cause: nested `plot_grid()` assembled bar chart and dot matrix without vertical axis alignment; different y-axis label widths caused progressive x-axis drift.
   - Fix: replaced nested `plot_grid()` with `align_plots()` (for panel width matching) + `ggdraw()`/`draw_plot()` (for pixel-accurate positioning). Both panels now start at the same x-coordinate with the same width.
   - For `set.size.show = FALSE` case: added `align = "v", axis = "lr"` to `plot_grid()`.

2. **Multi-set intersection bar color blending** (`R/ggupset.R`):
   - Added `.blend_colors()` helper: averages RGB values of participating set colors.
   - Single-set bars keep the set color; multi-set bars get a blended color.
   - `main.bar.color` parameter still overrides all bars to a uniform color.

3. **NAMESPACE updates** (`NAMESPACE`):
   - Added: `importFrom(cowplot, align_plots)`, `importFrom(cowplot, draw_plot)`, `importFrom(cowplot, ggdraw)`, `importFrom(grDevices, col2rgb)`, `importFrom(grDevices, rgb)`.

4. **Example figure regenerated** (`man/figures/example_upset.png`):
   - Used `set.seed(123)`, 250-gene universe, 150/130/110/100 sample sizes.
   - Shows proper alignment and blended bar colors.

5. **README updates** (`README.md`):
   - UpSet example: replaced `hsako$GeneID`-based example with self-contained reproducible dataset matching the figure.
   - Fixed undefined `res1`/`res2` → `resko1`/`resko2`.
   - Added `main.bar.color` override example.
   - Updated "Key features" to mention color blending.
   - Fixed citation: year 2025 → 2026, version 0.1.0 → 0.1.2.

6. **UPDATES.md**: Added v0.1.2.1 entry for alignment fix and color blending.

### Items completed in this session
- UpSet plot bar-to-dot alignment fix — Completed
- UpSet multi-set bar color blending — Completed
- Example figure regeneration — Completed
- README consistency (UpSet example, citation) — Completed
- PR #7 (`feature/upset-plot-v0.1.1`) updated with alignment fix — Completed
- PR #8 (`feature/enhancements-v0.1.2`) updated with all fixes — Completed
- PROJECT_HANDOFF.md and PROJECT_LOG.md created — Completed
