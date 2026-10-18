# PROJECT_HANDOFF.md — richR

## 1. Project overview and baseline

richR is an R package for functional enrichment analysis and visualization.
Kai Guo's `guokai8/richR` is upstream; `hurlab/richR` is our primary development line.

- Updated: 2026-10-17 CDT.
- Fork version: 0.1.4 (`DESCRIPTION` at `origin/master`).
- Verified code baseline: `19f38fb997ea8ae5236e1d52b6b3a0cbefea5493`.
- Wrap-up workspace: `.moai/worktrees/fork-wrapup`, branch `WT-fork-wrapup`, based on `origin/master`.
- Primary checkout remains on `upstream-pr/richr-v0.1.3` at `eb4c9da`; it is a PR preparation branch, not the current fork baseline.
- The earlier handoff's master-promotion and PR #9 instructions are superseded by this update. Do not repeat its reset/force-push plan.

## 2. Current state

| Item | Verified state |
|---|---|
| Fork release | Version 0.1.4 is already on remote `master` |
| Upstream sync | `19f38fb` merged upstream commits `590efdb` and `ef0fb86`; both `idconvert()` definitions use `AnnotationDbi::mapIds()` |
| Upstream-only commits | `git rev-list --count HEAD..upstream/master` returned `0` on the code baseline |
| Fork improvements | richUpset alignment/color blending, sna dependency, regenerated Rd documentation, fork README/figures remain in the fork |
| Upstream contribution | [PR #11](https://github.com/guokai8/richR/pull/11) is OPEN at `eb4c9da`, with 3 commits and 24 files |
| Earlier PRs | PR #7 and PR #9 are CLOSED and unmerged; follow PR #11 instead |
| Harness cleanup | PR #11 and the fork baseline have no tracked MoAI/agent harness paths in the inspected path families |
| Local harness files | Retained locally and excluded through `.git/info/exclude`; `.moai-backups/` was added to that local exclusion list |
| CI | The baseline passed R-CMD-check on Ubuntu, macOS and Windows; each reported `Status: 4 NOTEs` |

MoAI/agent harness files must never be committed to richR branches. This includes
`.moai/`, `.moai-backups/`, `.claude/`, `.agents/`, `.agency/`, `.codex/`,
`.git_hooks/`, `.mcp.json`, `.worktreeinclude`, `AGENTS.md`, and `CLAUDE.md`.
Keep machine-local exclusion rules in `.git/info/exclude` rather than committing
new harness wiring. Existing R-package ignore rules were left unchanged.

## 3. Verification evidence

### Claim

Both changed source definitions convert human Entrez IDs without attaching
AnnotationDbi. The upstream sync has landed in the fork, and its CI completed.

### Evidence

- `git show -s --format='%h %p %s' 19f38fb`:
  `19f38fb 3d4a496 ef0fb86 Merge upstream guokai8/richR master: use AnnotationDbi::mapIds in idconvert()`.
- `git rev-list --count HEAD..upstream/master`: `0`.
- `gh pr view 11 -R guokai8/richR --json state,headRefOid,commits,files`:
  `OPEN`, head `eb4c9daa39fa7002842d4e54cec620f5ac8dd3c6`, 3 commits, 24 files.
  Filtering its file paths for harness directories and root harness files returned `[]`.
- `git ls-files` for the harness path families listed above returned no paths.
- `git check-ignore -v .moai-backups/` in the primary checkout:
  `.git/info/exclude:30:/.moai-backups/ .moai-backups/` (fields separated by whitespace).
- `gh run view 39079967110 -R hurlab/richR --json jobs`:
  `macOS-latest (release)`, `ubuntu-latest (release)`, and `windows-latest (release)`
  each returned `conclusion: success`, including the package-check step.
- `gh run view 39079967110 -R hurlab/richR --log`:
  all three platforms reported `Status: 4 NOTEs`.
  [CI run](https://github.com/hurlab/richR/actions/runs/39079967110).

The direct runtime check used `timeout 120 Rscript --vanilla -e` with the following
R expression, from the wrap-up worktree:

```r
.libPaths(c(.libPaths(), "/home/juhur/R/x86_64-pc-linux-gnu-library/4.5"))
cat("R:", as.character(getRversion()), "\n")
stopifnot(requireNamespace("org.Hs.eg.db", quietly = TRUE))
e <- new.env(parent = baseenv())
sys.source("R/species.R", e)
stopifnot(!"package:AnnotationDbi" %in% search())
for (f in c("R/misc.R", "R/species.R")) {
  for (x in parse(f)) {
    if (is.call(x) && identical(x[[1]], as.name("<-")) &&
        identical(x[[2]], as.name("idconvert"))) eval(x, e)
  }
  actual <- e$idconvert("human", c("7157", "672"), "ENTREZID", "SYMBOL")
  stopifnot(identical(unname(actual), c("TP53", "BRCA1")))
  cat(f, ": 7157=", actual[[1]], " 672=", actual[[2]],
      "; AnnotationDbi attached=", "package:AnnotationDbi" %in% search(), "\n", sep = "")
}
cat("PASS: both source definitions convert IDs without attaching AnnotationDbi\n")
```

Observed output (excluding the annotation database's mapping messages):

```text
R: 4.6.1
R/misc.R: 7157=TP53 672=BRCA1; AnnotationDbi attached=FALSE
R/species.R: 7157=TP53 672=BRCA1; AnnotationDbi attached=FALSE
PASS: both source definitions convert IDs without attaching AnnotationDbi
```

### Baseline-attribution

Local runtime checks ran against the files at `19f38fb` in the isolated wrap-up
worktree during this session. GitHub's recorded CI run targets exactly the same
commit; it was inspected in this session, not rerun locally.
The runtime used R 4.6.1, prioritizing its normal libraries and adding the existing
R 4.5 user library last for annotation packages. No packages were installed.

### Gaps

- Full local package installation and R CMD check were not run. The default R 4.6
  libraries could not load AnnotationDbi, org.Hs.eg.db or sna. The direct check
  succeeded only after adding the existing R 4.5 library as a fallback.
- The direct check evaluates each changed function separately with species helpers;
  it does not test loading the complete installed package or every identifier type.
- No new enrichment integration tests or additional example figures were produced.

### Residual-risk

CI still reports four NOTE categories: hidden files, top-level packaging files,
unused declared Imports, and code-analysis findings about unqualified names.
These are observed CI findings, not a claim of an entirely clean CRAN check.
Cross-version library fallback is a limitation of the direct local check.

## 4. Remaining work and restart

The requested fork sync, upstream PR cleanup, local backup exclusion, and state
reconciliation are complete. PR #11 still awaits upstream action. No automatic
review reminder, branch deletion, dependency installation, or new feature work is
scheduled by this handoff.

For subsequent fork work, create a fresh isolated worktree from remote
`origin/master`. Do not resume from the older upstream PR branch by accident.
The user authorized committing these two wrap-up documents and pushing them to
`hurlab/richR` master. Retain the worktree until the remote update is verified.
Historical entries in PROJECT_LOG.md describe their original session state.
