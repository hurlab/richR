# richR Updates

## Version 0.2.0 (unreleased)

### Significance cutoffs (changes reported results)

- Every function that reports or plots enrichment results now applies one
  significance rule: a term is kept when its raw p-value is at most `pvalue`
  **and**, unless `padj` is `NULL`, its adjusted p-value is at most `padj`.
  Both comparisons use `<=`, and a row with a missing p-value (or a missing
  adjusted p-value while `padj` is set) is dropped. In 0.1.x the comparisons
  were strict (`<`) and a non-`NULL` `padj` replaced the raw cutoff instead
  of adding to it. The rule covers the ORA functions (`enrich()`,
  `richGO()`, `richKEGG()`, `richLevel()`, `richDAVID()`), `richGSEA()`,
  `compareResult()`, `compareGSEA()`, `richCompareDot()`, `getPathways()`,
  `filterPathways()`, `richGSEAcurve()` and the plot functions.
- The defaults are now `pvalue = 0.05` and `padj = 0.05` everywhere,
  including the plot functions; in 0.1.x `padj` defaulted to `NULL`. A
  default call therefore reports fewer terms than in 0.1.5, and on random or
  weak input it can report none. To get the old behaviour back:
  - `padj = NULL` filters on the raw p-value only, as the 0.1.x defaults did;
  - `pvalue = 1` filters on the adjusted p-value alone, as a 0.1.x call with
    `padj` set did;
  - `pvalue = 1, padj = 1` disables significance filtering (size filters
    still apply).
- Cutoffs are validated with one rule everywhere: a single finite number in
  (0, 1]. `richGSEA()` now rejects `pvalue = 0` and `padj = 0`, which it
  accepted in 0.1.x.
- `compareResult()` and `compareGSEA()` no longer apply the extra
  `Padj < 1` / `Pvalue < 1` guards of 0.1.x; terms with a p-value or an
  adjusted p-value of exactly 1 are kept when the cutoffs allow them.
- `getPathways()` and `filterPathways()` follow the same rule: both
  `padj_cutoff` and `pval_cutoff` apply together (in 0.1.x `pval_cutoff`
  was used only when `padj_cutoff` was `NULL`). `padj_cutoff = NULL` filters
  on the raw p-value only, and `pval_cutoff = NULL` applies no raw cutoff.
- `richCompareDot(include.all = TRUE)` draws every input row without
  applying `pvalue` or `padj`; filter in `compareResult()` instead.
- `richVolcano()` does not filter rows; with the new default `padj = 0.05`
  its dashed reference line is drawn at the adjusted cutoff, where in 0.1.x
  it was drawn at `pvalue`. Pass `padj = NULL` to draw it at `pvalue`.
- Filtering at analysis time is permanent: rows removed by `richGO()`,
  `richKEGG()`, `richGSEA()` and the other analysis functions are not stored
  in the result, so a later plot or `compareResult()` call with cutoffs of 1
  cannot bring them back. Run the analysis with `pvalue = 1, padj = 1` when
  you want to choose cutoffs later.
- `richGSEAcurve()` applies `pvalue` and `padj` only when it picks the
  pathways itself; its explicit `pathways` and `pathway_pattern` selections
  are not significance-filtered, by design, and choose among the rows the
  `GSEAResult` holds.
- Printing a `richResult` or `GSEAResult` shows the stored cutoffs with the
  inclusive rule, for example `(pvalue <= 0.05, padj <= 0.05)`; in 0.1.x a
  `richResult` printed `(pvalue < 0.05 )` and never the adjusted cutoff.
  `summary(<GSEAResult>)` counts the rows passing the stored cutoffs instead
  of a fixed `padj < 0.05`.
- `richDAVID()` is deprecated and warns on every call: it depends on
  RDAVIDWebService, which is no longer available from Bioconductor, and it
  will be removed in the next minor release. Use `enrich()`, `richGO()` or
  `richKEGG()` with a locally built annotation instead. Its documentation now
  states that Benjamini, DAVID's adaptive FDR and Bonferroni are distinct
  quantities and that DAVID applies `count = minSize` on the server, so
  `keepRich` cannot restore smaller terms.

### ORA tested family (changes adjusted p-values)

- ORA now adjusts p-values over the terms whose annotated size lies in
  [`minGSSize`, `maxGSSize`] (inclusive) and that share at least one gene
  with the query, the convention of clusterProfiler and DOSE. In 0.1.x the
  size bounds were applied after adjustment, so terms outside them still
  entered the adjustment; adjusted p-values can therefore be smaller than in
  0.1.5 when the annotation holds terms outside the size bounds. `minSize`,
  `maxSize`, `keepRich` and the cutoffs remain reporting filters and never
  change an adjusted value.
- New argument `padj.family` on `enrich()`, `richGO()`, `richKEGG()` and
  `richLevel()` (and through `...` on `batchEnrich()`). The default
  `"overlap"` is the family above; `"all"` adjusts over every size-eligible
  term, zero-overlap terms contributing p = 1 (they are never reported).
  The default family is a convention, not an FDR guarantee: with a universe
  of 100 single-gene terms and one selected gene under the null
  (`minGSSize = 1, minSize = 1`), the only overlapping term has p = 0.01,
  BH over that one term gives 0.01, and a null term is reported every time;
  BH over all 100 terms gives 1. Use `padj.family = "all"` when that matters.
  `richGO()` adjusts within each ontology and `richLevel()` within each
  level. `enrichmentInfo()` records the family in `testing$family_policy`
  and the argument in `options$padj.family`.
- `keepRich = TRUE` now exempts fully rich terms (`RichFactor == 1`) from
  the `minSize` filter only. It no longer admits terms below `minGSSize`,
  which are outside the tested family. With the default sizes
  (`minGSSize = 10`, `minSize = 2`) every tested term already has at least
  10 annotated genes, so `keepRich` has no effect unless `minGSSize` is
  lowered below `minSize`.

### Annotation identity (changes ORA results and labels)

- One identifier rule for every annotation builder, ORA function and
  `richGSEA()`: column 1 is the gene; the identifier is the first of `GOALL`,
  `PATH`, `Module` or `Term`, otherwise column 2; the label is `Annot`.
  `Term` (like `GOALL`, `PATH` and `Module`) is the identifier column in an
  annotation frame, while in an ORA result `Annot` is the identifier and
  `Term` the label (in a GSEA result, `pathway` and `Term`). The help pages of
  `enrich()`, `richGSEA()`, `readGMT()` and `buildAnnotFromList()` state this
  convention.
- `enrich()` results now show labels: `Term` holds the label taken from the
  annotation frame's `Annot` column (the identifier where no label exists)
  and `Annot` the identifier. In 0.1.x `enrich()` reported the identifier in
  both. The result and detail columns are unchanged.
- Migration: in 0.2.0 ORA, as in `richGSEA()`, `Term` is the identifier
  column of an annotation frame and `Annot` the label, whatever the column
  order. In 0.1.x, ORA keyed a three-column frame by column 2. Two layouts
  with column 2 = `Annot` therefore change keying:
  - `(GeneID, Annot = label, Term = ID)`: 0.1.x keyed it by the label, so
    identifiers sharing a label were merged; 0.2.0 keys it by `Term`, the
    identifier, and tests them separately. Gene sets, counts and result rows
    can change; no action is needed.
  - `(GeneID, Annot = ID, Term = label)` (named like a result): 0.1.x keyed it
    by the identifier in `Annot`; 0.2.0 keys it by the label in `Term`, so
    identifiers sharing a label are merged into one gene set, with a
    conflicting-labels warning. For example, with
    `Annot = c("T1", "T1", "T2", "T2")` and `Term = c("same", "same", "same", "same")`,
    0.1.x tested two gene sets, `T1` and `T2`; 0.2.0 tests one, `same`, and
    reports `Annot = "same"`, `Term = "T1"`. Rename the columns so that `Term`
    holds the identifier and `Annot` the label:
    `names(frame)[2:3] <- c("Term", "Annot")`.
- An identifier that carries different labels in one annotation takes the
  first label in row order, and ORA and GSEA warn once, naming the gene sets
  (`conflicting labels (Annot) for gene set(s) ...`) and the column
  convention.
- Builder messages count gene sets by the identifier column rather than by
  column position, and the set-name builders (`readGMT()`, `buildAnnotFromList()`, set-name
  `buildMSIGDB()` collections) write the set name to both `Term` and `Annot`.
- `buildAnnot(anntype = "Reactome")` now keeps the pathway identifier: the
  annotation has the columns `GeneID`, `Term` (the `R-HSA-` identifier) and
  `Annot` (the pathway name without the species prefix). Reactome pathways
  that share a name are tested as separate gene sets; in 0.1.x they were
  merged by name (with reactome.db 1.96.0, 2845 human pathways with genes
  became 2836 gene sets, because 18 pathways share 9 names). ORA results
  report the identifier in `Annot` and the name in `Term`; GSEA results in
  `pathway` and `Term`. No pathway is dropped: every human pathway with at
  least one gene is kept.
- `buildAnnot(anntype = "Reactome")` no longer needs reactome.db to be
  attached; in 0.1.x it failed with `object 'reactomePATHNAME2ID' not found`
  unless `library(reactome.db)` had been called first.
- ORA plots draw one category per identifier. Where two terms share a label,
  the identifier is appended (`"label [ID]"`), in `richBar()`, `richDot()`,
  `richNetplot()`, `richNetwork()`, `richHeatmap()`, `richCompareDot()`,
  `richClusterDot()` and the other ORA plot functions; in 0.1.x such terms
  were merged into one bar or point, or the plot failed while naming network
  nodes. The labels are computed after `short = TRUE` shortening, so labels
  that become equal once shortened stay apart (also for GSEA plots).
- `richHeatmap()`, `richCompareDot()` and `richClusterDot()` match rows
  across results by identifier (`Annot`) instead of by label; a legacy
  result without `Annot` is still matched by `Term`.

### KEGG labels and provenance

- `richKEGG(builtin = FALSE)` and `buildAnnot(anntype = "KEGG", builtin =
  FALSE)` read pathway labels from live KEGG. When that listing fails, or
  lacks the reference pathway `04148`, richR now warns
  (`KEGG pathway labels: <reason>; using the bundled KEGG snapshot.`) and uses
  the bundled labels. In 0.1.x a failed listing stopped with an error and a
  listing without `04148` switched to the bundled labels silently. There is
  no option to silence the warning.
- KEGG and KEGG module annotations, and `richKEGG()` results, record the
  label source actually used, separately from the membership provider:
  `enrichmentInfo(<Annot>)$provenance$annotation$label_source` and
  `enrichmentInfo(<richResult>)$provenance$labels` hold `source`
  (`"KEGGREST"` or `"bundled"`), `release` and `retrieved`. A fallback is
  never recorded as a live retrieval.
- For live KEGG labels the KEGG version is recorded in `label_source$release`
  and in `source_version`: the release KEGG states, verbatim (for example
  `Release 116.0+/10-28`); otherwise `last update YYYY-MM-DD`, built from the
  update date KEGG states (for example `last update 2026-10-28` from
  `Last update 2026/10/28`); otherwise `NA`. KEGG stopped stating a release
  number in its database information; richR never infers one. For bundled
  labels, a fallback or KEGG modules the KEGG version is `NA`.
- `source_version` is recorded only where the source states a version (the
  msigdbr `db_version`, the KEGG version above) and is `NA` otherwise; the
  provider package and its version stay in their own fields.

### Gene-set fingerprints

- Every annotation built by richR (`buildAnnot()`, `buildOwn()`,
  `buildMSIGDB()`, `readGMT()`, `buildAnnotFromList()`) and every ORA and
  GSEA result records two gene-set fingerprints in
  `enrichmentInfo(x)$provenance$annotation$fingerprint`: `membership`, the
  sha256 of the identifier-gene pairs, and `full`, the sha256 of the pairs
  plus each identifier's label, with `algorithm = "sha256"`,
  `canonicalization = "richR-geneset-v1"` and the numbers of pairs, gene sets
  and genes. Equal `membership` means the same gene sets were tested; equal
  `full` means the same labels would be printed too. Both ignore row order,
  duplicated rows, row names, extra columns and factor versus character
  columns. Results fingerprint the annotation actually passed to the
  analysis, so an annotation edited after it was built is fingerprinted as
  used. `digest` is a new import.

### GSEA

- `richGSEA()` keys gene sets by annotation identifier. Distinct gene sets
  that share a label were merged into one set in 0.1.x; they are now tested
  separately. `@result$pathway` holds the identifier (for example a KEGG or
  GO ID), and a new `Term` column, directly after `leadingEdge`, holds the
  label. Code that compares `pathway` with pathway names must switch to
  `Term`, for example `subset(result(res), Term == "Cell cycle")`.
- `richGSEAcurve(pathways = )`, `compareGSEA(pathway = )`,
  `richGSEAplot(term = )` and `getGenes(<GSEAResult>, term = )` accept an
  identifier, a KEGG identifier with or without its organism prefix
  (`"hsa04110"` or `"04110"`), or a unique `Term` label. A label shared by
  several pathways, or a value that matches nothing, is an error.
  `richGSEAplot()` also still draws a `GSEAResult` saved from 0.1.x, whose
  `pathway` column holds labels.
- Plots of a `GSEAResult` show the `Term` label; where two pathways share a
  label, the identifier is appended (`"label [ID]"`) so they are drawn as two
  categories with their own values.
- Removed the unused `minGSSize` and `maxGSSize` arguments of the internal
  GSEA routine: they were validated and recorded but never passed to fgsea,
  and the exported `richGSEA()` methods never exposed them. The GSEA family
  is bounded by `minSize` / `maxSize` on the gene-set members found in the
  ranked list.

### Bundled KEGG data (changes KEGG names, levels and counts)

- The bundled `kegg.db`, `module` and `path` were rebuilt from live KEGG on
  2026-11-05. KEGG stated no release number that day, only its last update,
  so the KEGG version is recorded as `last update 2026-10-28`. The rebuild
  depended on KEGG on that date and is not byte-reproducible. Each object
  records its retrieval date, KEGG version and source URL in
  `attr(<object>, "richR_source")`, and its help page states the same values
  and the row count. A KEGG label fallback (`richKEGG(builtin = FALSE)` or
  `buildAnnot(anntype = "KEGG", builtin = FALSE)` when the live listing fails)
  now reports the snapshot's retrieval date in `enrichmentInfo()`.
- `kegg.db` grew from 567 to 587 pathways (21 added, `03269` removed, 12
  renamed, for example `04145` "Phagosome" is now "Phagocytosis"), `module`
  from 473 to 583 modules (112 added, `M00166` and `M00167` removed, 32
  renamed) and `path` from 580 rows to 587 pathways (8 added, 11 `Level3` and
  1 `Level2` values changed). `richKEGG()`, `richLevel()` and `richGSEA()`
  results on KEGG annotations show the new names and levels.
- `path` now has exactly one row per pathway. The earlier snapshot listed
  `00541` twice ("O-Antigen nucleotide sugar biosynthesis" under "Glycan
  biosynthesis and metabolism", and "Biosynthesis of various nucleotide
  sugars" under "Carbohydrate metabolism"), so `richLevel()` counted the
  genes of `00541` under both rows, while `richKEGG()` and `richGSEA()` used
  only the first. `richLevel()` now counts each pathway under one hierarchy
  row; `00541` is the only pathway affected. Where KEGG lists a pathway under
  several hierarchy rows, `data-raw/path.R` keeps the first in KEGG's order.

### Tested terms, comparing and storing results

`@result` and `compareResult()` are unchanged: every entry point reports the
same rows, columns, values and order as before, and `compareResult()`,
`richCompareDot()`, `batchEnrich()` and the plots do not read the new tested
table.

- `richResult` and `GSEAResult` objects have a new slot, `tested`, that
  records what the analysis tested, not only what it reported. For ORA it
  holds one row for every term of the tested family that shares at least one
  gene with the query; for GSEA one row for every pathway fgsea returned. Its
  columns are those of `@result`, with the values the analysis computed, plus
  a logical column `reported` that is `TRUE` exactly for the rows `@result`
  holds. Terms removed by `minSize`, `maxSize`, `keepRich` or the
  significance cutoffs stay in the table with `reported = FALSE`, so running
  with the default cutoffs no longer loses the statistics of the terms that
  did not pass.
- New `testedTerms(x)` returns the tested table as a data frame (zero rows
  when nothing was tested).
- `enrichmentInfo(x)$testing$zero_overlap_count` counts the tested terms with
  zero overlap, which the tested table does not list. It is non-zero only for
  ORA with `padj.family = "all"` and `0` otherwise; the tested table's rows
  plus `zero_overlap_count` equal `testing$tested_count`. The metadata of
  `richResult` and `GSEAResult` objects now has `schema_version = 2`;
  `Annot` metadata and the unavailable record keep `schema_version = 1`.
- Objects without a recorded tested table (saved by an earlier richR, built
  with `as(<data.frame>, "richResult")`, by `richDAVID()` or by `new()`) are
  reported as unavailable, never reconstructed: `testedTerms()` returns
  `NULL` with the warning "tested table unavailable", and no function treats
  `@result` as the tested set or fills in zeros. Objects saved before 0.2.0
  fail `validObject()` ("slots in class definition but not in object:
  \"tested\""), but richR functions do not call `validObject()` on their
  inputs and keep accepting them.
- `filter()`, `arrange()`, `mutate()`, `slice()`, `select()`, `rename()` and
  `group_by()` on a result object change only `@result`; the tested table and
  the `enrichmentInfo()` record stay those of the original analysis, so a
  filtered object is a view of that analysis, not a new one.
- New `bindRich()` binds a list of `richResult` objects, or of `GSEAResult`
  objects, into one data frame with one row per group and identifier, taken
  from the union of every group's `@result` and tested table. Its columns are,
  in order: `group` (a factor with the groups in input order), the identifier
  (`Annot` for ORA, `pathway` for GSEA), `Term`, `state`, `execution`, then
  the other columns of the groups' results and tested tables. `execution` is
  `"success"` (the group's `@result` has rows), `"empty"` (it has none) or
  `"failed"` (a `NULL` entry, as `batchEnrich()` returns for a failed group).
  `state` is `NA` for a failed group and otherwise `"reported"` (in the
  group's `@result`), `"tested_not_reported"` (in its tested table only),
  `"not_tested"` (absent from a complete tested table),
  `"untested_or_zero_overlap"` (absent, and the group counted zero-overlap
  terms it does not list) or `"unavailable"` (absent from a group whose tested
  table or metadata is not recorded, so what it tested is unknown). Each row's
  `Term` and statistics come from that group's own record, and are `NA` in
  every other state: `bindRich()` applies no cutoff, re-adjusts no p-value,
  computes no statistic and never writes `0` for a missing value. Groups are
  named as in `compareResult()`; ORA and GSEA results cannot be mixed.
- New `writeRich(x, path)` stores a `richResult` or `GSEAResult` as a folder:
  `manifest.json` plus UTF-8 tab-separated tables `result.tsv`, `tested.tsv`
  (when a tested table is recorded), `detail.tsv` (`richResult`) or
  `input.tsv` (`GSEAResult`, the ranked genes and scores in stored order).
  Each table has a header row and one line per row. In text, factor and list
  cells a backslash is written `\\`, TAB `\t`, newline `\n` and carriage
  return `\r`; a missing value of any class is written `\N`, so `NA`, `""`
  and the text `"NA"` stay distinct; identifiers keep their text (`"00010"`
  is not read back as `10`). Doubles are written with 17 significant digits,
  which on Linux read back bit-identical (the default 15-digit text form
  changes most of them); `NaN`, `Inf` and `-Inf` are spelled out, and list
  cells such as a GSEA `leadingEdge` are written as JSON arrays. A table with
  no columns is recorded in the manifest only. The manifest records the
  format name `richR-result`, format version `1`, the class, every other
  slot, the `enrichmentInfo()` record and, per table, its file, column names
  and classes, row names and row count. `writeRich()` checks everything before
  it writes: it refuses an existing non-empty folder unless
  `overwrite = TRUE`, and even then any folder it did not write itself; a
  grouped `@result` (convert it with `x@result <- as.data.frame(x@result)`);
  an `@result` without its identifier or `Term`; and columns of other
  classes.
- New `readRich(path)` rebuilds the object `writeRich()` stored, every value
  exact: the tested table is unavailable exactly when it was in the stored
  object, and `enrichmentInfo()` returns the stored record unchanged, never
  upgraded or re-dated. Version policy: the folder format version must be `1`
  and, when the folder records metadata, its `schema_version` `1` or `2`;
  `schema_version` 1 is read as it is, without conversion, and a missing,
  non-integer, older or newer version is an error naming the found and
  supported versions. A folder written from an object without metadata is
  read back without metadata.
  `readRich()` also stops, naming the folder, file, column and first offending
  row, on a missing manifest, a missing or unlisted table file, a header or
  row count that differs from the manifest, a line with more or fewer fields
  than the manifest's columns, an unknown escape sequence or a cell that does
  not parse to its column's class (a double only in the form `writeRich()`
  writes it); it never returns a partially read or silently coerced object.
  The metadata and slots in the manifest are read as plain data only (`NULL`,
  logical, integer, double, character and lists, with names): a manifest
  holding a function, call, environment, namespace, raw or serialized bytes,
  an S4 object or a classed value is refused, naming the field, so reading a
  folder never rebuilds code or loads a package.

### Identifier mapping, background size and sensitivity grids

Every existing entry point called without the new `universe_policy`
argument reports the same rows, statistics and metadata fields as before;
`idconvert()`, `buildAnnot()`, `richLevel()`, `compareResult()` and the plots
are unchanged.

- New `mapGenes(x, mapping, from, to, multiple = c("drop", "first"))`
  converts a gene vector, or a ranking (a named numeric vector), through a
  mapping table (a data frame whose `from` and `to` columns are character or
  factor; a factor is read as its labels, any other column type is an error
  naming the column), an `OrgDb` object, or a species name listed by
  `showData()` (its OrgDb package must be installed). It never selects a
  mapping silently: every input entry gets exactly one decision, decided in
  this order — `invalid` (missing or blank identifier; for a ranking also a
  missing or non-finite score), `unmapped` (no candidate), `ambiguous`
  (several candidates under `multiple = "drop"`, the default), then, among the
  entries still kept, `mapped` for the first entry of each target and
  `merged` for later entries of the same target. With `multiple = "first"`
  the first candidate in source order is kept and the `reason` names the
  candidate count. A source pair whose target is missing or blank is not a
  candidate, and a repeated pair counts once. For a ranking, entries sharing a
  target are kept once when their scores are identical and are all dropped as
  `conflict` when they differ: `mapGenes()` never chooses among entries by
  score value or magnitude. The result holds each kept target once, in the
  order of its first kept entry (for a ranking with the input scores
  unchanged, so it can go straight to `richGSEA()`); under `"drop"` it does
  not depend on the row order of the source. When every entry is dropped the
  result is empty and a warning states the count per decision. Malformed
  calls are errors naming the argument and the offending value.
- The conversion is recorded in the attribute `"richR_mapping"` of the
  result, as plain data (named lists and atomic vectors, no `class` or
  `row.names`): `decisions`, a column list with one element per input entry
  (`position`, `input_id`, `candidates`, `target`, `decision`, `reason`, and
  `score` for a ranking), and `source`: `kind` (`"table"` or `"OrgDb"`),
  `from`, `to`, `multiple`, the richR version and, for an OrgDb, its package,
  installed version and `metadata = list(name, value)` from
  `AnnotationDbi::metadata()`; for a table, its row count and a fingerprint
  `list(algorithm = "sha256", canonicalization = "richR-idmap-v1", hash,
  n_rows)`. The hash is the SHA-256 of the bytes of the text made of the line
  `richR-idmap-v1` and one line `<from>\t<to>` per table row in row order, each
  ended by LF, fields in UTF-8 with backslash, TAB, LF and CR written `\\`,
  `\t`, `\n`, `\r` (the `richR-geneset-v1` rule) and a missing value `\N`.
  Values a source does not state are `NA`.
- `enrich()`, `richGO()`, `richKEGG()` (and `batchEnrich()`) and `richGSEA()`
  keep a matching mapping record as `enrichmentInfo(x)$input$mapping`. A
  record whose targets (and scores) no longer equal the submitted vector —
  for example after `v[1] <- "X"` — is stored as
  `list(available = FALSE, reason = "mapping record does not match the
  submitted input")` with a warning; an input without a record (also after
  `unique()`, `c()` or `[`, which drop attributes) records none. ORA results
  also record `input$occurrences`, the number of times each identifier of
  `input$ids` occurred in the valid input; GSEA results record 1 per ranked
  gene. The record survives `writeRich()` / `readRich()` unchanged. The new
  fields are additive within `schema_version = 2`.
- New `mappingAudit(x)` lists one decision per identifier, read only from the
  recorded metadata and mapping record — never reconstructed from `@result`,
  `@gene` or `@detail`. For a `mapGenes()` vector it returns the decision
  table as a data frame. For a `richResult` it returns the columns `role`,
  `input_id`, `mapping`, `id`, `occurrences`, `decision`, `reason`: in the
  plain layout one `query` row per unique valid query identifier (with its
  recorded count) and one `invalid` row per missing or blank entry; in the
  joined layout (the input carried a mapping record) one `query` row per
  mapping entry, in input order, joining the original identifier and its
  mapping decision to the engine decision of its target (`not_submitted` for
  dropped entries; `occurrences` on the target's `mapped` row only). Query
  decisions are `effective`, `not_annotated`, `outside_background` and
  `invalid`; `background` rows are `effective` or `not_annotated`. The counts
  reconcile with the engine: the distinct `effective` query ids number
  `testing$effective_n`, the `effective` background rows number
  `testing$effective_N`, and the recorded occurrences sum to
  `input$original_count`. For a `GSEAResult` all rows are `ranked`
  (`ranked_annotated`, `ranked_not_annotated` — the gene stays in the ranking
  but belongs to no gene set — or `not_submitted`), with the gene's `score` in
  the result's `input` slot and `tied` when another ranked gene has the same
  score. Results recorded before this version give the plain layout with
  `occurrences` `NA` and are read as `universe_policy = "annotated"`; objects
  without metadata, with the unavailable record, or vectors without a mapping
  record (or with a `richR_mapping` attribute that is not a `mapGenes()`
  record) give `NULL` with a warning "mapping audit unavailable".
- New argument `universe_policy = c("annotated", "measured")` on the `enrich()`,
  `richGO()` and `richKEGG()` methods (and through `batchEnrich()`'s `...`).
  The default, `"annotated"`, keeps `N` = the unique valid universe genes that
  are annotated in the annotation the analysis tests (for `richGO()`, the
  selected ontology) and `n` = the query genes among them. `"measured"`
  requires an explicit `universe` and also counts its tested-but-unannotated
  genes: `N` = the unique valid universe genes, `n` = the unique valid query
  genes in the universe. Term sizes, overlaps and the tested family are the
  same under both; `FoldEnrichment` and `zscore` follow the counted `N` and
  `n`. Under `"measured"`, `background$effective_ids`, `input$effective_ids`
  and `testing$effective_N` / `effective_n` hold the counted sets, while
  `mapped_ids`, `unmapped_ids` and `excluded_by_universe_ids` keep their
  annotation-membership meaning. The choice is recorded as
  `options$universe_policy`. `"measured"` without a `universe`, or any other
  value, is an error; a universe or query with no annotated gene keeps the
  existing errors. `richLevel()` takes no `universe` and gets no
  `universe_policy`; GSEA has no universe.
- New `richSensitivity(x, object, backgrounds, maxGSSize = c(500, 1000), ...)`
  runs one independent `enrich()` per cell of the named `backgrounds` levels
  (`NULL` for the annotation-derived background or a character vector used as
  `universe`) crossed with the `maxGSSize` bounds, every other argument
  identical; a cell is never derived by filtering another cell (a filtered run
  has different adjusted p-values). Cells are named `"<background>/<bound>"`.
  It returns `cells` (execution `success`, `empty` or `failed` with the error,
  the cell's warnings — each distinct message is signalled once —,
  `background_size`, `effective_N`, `effective_n`, `tested_count`),
  `results` (the per-cell results, `NULL` when failed), `states`
  (`bindRich()` of the results with `background` and `maxGSSize` columns),
  `summary` and `design`. The summary has one row per term that at least one
  cell tested (a term no cell tests does not appear): `cells`, `ran`
  (success or empty), `tested` (reported or tested-not-reported), `reported`,
  `not_tested`, `undetermined` (untested-or-zero-overlap or unavailable),
  `support = reported / tested` (`NA` when untested), the ranges of the tested
  cells' `Padj` and `FoldEnrichment`, and the best and worst `Pvalue` rank
  among each tested cell's tested terms; values no cell records are `NA`,
  never `0`. `design` records the backgrounds as passed, their sizes, the
  bounds, the engine, the `...` arguments (including `universe_policy`, which
  applies to every cell), the criterion (state `reported` under the cells'
  cutoffs), the definitions of `ran`, `tested` and `support` and the richR
  version, and reproduces the grid. The summary is descriptive stability
  evidence across prespecified settings — not a p-value, a confidence
  interval, an independent replication or a combined test; no statistic is
  combined across cells. Invalid backgrounds or bounds, `universe` or
  `maxGSSize` in `...`, `universe_policy = "measured"` with a `NULL` level, a
  numeric ranking or unnamed numeric identifiers (pass them as character), or
  a grid whose every cell fails are errors. The grid is ORA only.
- Known residual: `idconvert()` and the conversions inside `buildAnnot()` and
  `buildMSIGDB()` keep their first-candidate behaviour
  (`AnnotationDbi::mapIds(multiVals = "first")`); use `mapGenes()` for a
  recorded conversion.

### Plots show what they draw (changes default plot output)

A plot function now only displays an enrichment result: it reports the rows
it leaves out and no longer draws invented or mislabelled values. Plot
defaults, argument names, aliases and layouts are unchanged; every other plot
output is as before. See `?richR_plot_selection`.

- **Selection message.** The table plots (`richBar()`, `richDot()`,
  `richNetplot()`, `richNetwork()`, `richNetmap()`, `richLollipop()`,
  `richCircle()`, `richNES()`, `richScatter()`, `richECDF()`,
  `richTermSim()`, `richGeneDot()`, `richGeneHeat()`, `richGeneBar()`,
  `richVolcano()`, `richHeatmap()`, `richCompareDot()` and their aliases)
  apply their cutoffs and `top` exactly as before, but previously did so
  silently: a result holding 6 terms could be drawn with 3, and a
  `compareResult()` table built with `pvalue = 1` was re-filtered by
  `richCompareDot()`'s defaults. When the selection drops rows, a call now
  signals one message of class `richR_plot_selection`, for example
  `richBar: drawing 3 of 6 rows (Pvalue <= 0.05 and Padj <= 0.05; first 50 in
  input order)`; none is signalled when every row is drawn. `richNES()` and
  `richECDF()` name their order rule (`by NES, descending` / `ascending`),
  `richHeatmap()` counts identifiers over all results, and `richNetmap()`
  signals one message counted against the results' own rows. The drawn rows
  are unchanged, so the earlier silent output is
  `suppressMessages(<same call>)`. Each plot's help page gives the call that
  draws every row: `pvalue = 1, padj = NULL, top = Inf` for most plots,
  `top = Inf` for `richVolcano()` and `richHeatmap()`, `include.all = TRUE`
  for `richCompareDot()`, and a finite `top` of at least the total row count
  for `richNetmap()` (`top = Inf` is not supported there).
- **GSEA lollipop, gene-bar and scatter plots draw NES.** For a `GSEAResult`,
  `richLollipop()` and `richGeneBar()` previously drew a constant
  `RichFactor = 1` for every pathway (axis title `RichFactor`), and
  `richScatter()` drew every pathway at `log2` fold enrichment 0 under the
  title `log2(Fold Enrichment)` against `log10(Pathway Size)`. They now draw
  each pathway's `NES` (axis title `NES`); `richScatter()` draws `NES`
  against `log10(size)`, titled `log10(Gene-set size)`. The plot data no
  longer carries the invented `RichFactor` and `Annotated` columns. ORA input
  is unchanged.
- **Effect titles without `FoldEnrichment`.** `richVolcano()` and
  `richScatter()` on ORA input without a `FoldEnrichment` column (results
  saved before richR 0.1.2) draw the same values as before, but previously
  titled them `log2(Fold Enrichment)` (volcano colour `log2(FE)`). The
  titles now say what is drawn: `log2(Rich Factor = Significant/Annotated)`
  (volcano colour `log2(Rich Factor)`); a frame without `Significant` and
  `Annotated` as well is drawn by `richVolcano()` with every point at 0,
  titled `No effect column: all points at 0` (colour `No effect column`).
  With `FoldEnrichment` present the titles are unchanged.
- **Missing fold changes in gene plots.** When `fc` holds no value for a
  member gene of a drawn term, `richGeneDot()` and `richGeneHeat()`
  previously recorded and drew it as fold change 0, white like a measured
  zero. The gene's `fc_value` is now `NA` and it is drawn in the fixed colour
  `grey55` (`#8C8C8C`): a smallest-size point in `richGeneDot()`, a tile in
  `richGeneHeat()`, whose cells for genes not in a term keep `na.fill`
  (`grey90` by default). A caption names the encoding. Output is unchanged
  when every member gene has a fold change or `fc` is `NULL`. To draw the
  missing genes as 0 again, supply them: `fc[setdiff(genes, names(fc))] <- 0`.
- **Absent terms in `richHeatmap()`.** A term absent from one of the
  results was filled with `p = 1` and drawn as a measured non-significant
  value; its cell is now `NA` in the drawn data and in the
  `returnData = TRUE` table, and is drawn in the missing-value colour. Code
  reading the `returnData = TRUE` table must expect `NA` there;
  `tab[is.na(tab)] <- 1` restores the earlier table. The fill legend
  previously named the other p-value column (`-log10(Padj)` over p-values
  with the default `usePadj = FALSE`, and the reverse); it now names the
  column drawn, `-log10(Pvalue)` or `-log10(Padj)`, and the value column of
  the plot data carries that name. The row order is unchanged (an absent
  cell still counts as `p = 1` for ordering).
- Known limits, unchanged: for a `GSEAResult`, `richGeneDot()` and
  `richGeneHeat()` stop with an error when a selected pathway has an empty
  leading edge, and the network plots never draw such a pathway;
  `richNetwork()` and `richNetmap()` omit a selected term that shares no gene
  with another selected term.

### Requirements and packaging

- richR now requires R 4.6 or later (`Depends: R (>= 4.6)`). The package
  metadata follow Bioconductor conventions: non-empty `biocViews`, plus `URL`
  and `BugReports` fields.
- The `LazyData` field was removed. `data(path, package = "richR")` and the
  other `data()` calls work as before, and richR's own functions load the
  bundled data with `utils::data()`, so analysis results are unchanged.
  Direct access such as `richR::path`, which worked in 0.1.5, now fails with
  `'path' is not an exported object from 'namespace:richR'`; load the dataset
  with `data()` first. The same applies to `godata`, `kegg.db`, `module` and
  `pathway`. Four of the five load under their own name (`data(path)`,
  `data(godata)`, `data(module)`, `data(pathway)`); the KEGG table `kegg.db` is
  stored in `data/kegg.rdata`, so it is loaded with `data(kegg)`, which creates
  the object `kegg.db`. `data(kegg.db)` loads nothing and warns that the data
  set was not found.
- Added `data-raw/` scripts (`kegg.R`, `module.R`, `path.R`, `pathway.R`) that
  rebuild the bundled KEGG data and record the retrieval date and KEGG
  version. `kegg.R`, `module.R` and `path.R` were run for this release (see
  "Bundled KEGG data" below); `pathway.R` was not, and `godata` and `pathway`
  are unchanged, with unknown retrieval dates.
- `jsonlite` is a new import, used by `writeRich()` and `readRich()` for the
  folder manifest.
- Corrected the documentation of the bundled data: `godata` is a named list of
  GO terms mapped to themselves and their offspring, not a data frame of gene
  annotations.

### Bug fixes

- `compareGSEA()` with a `padj` cutoff no longer stops with an error: in
  0.1.x it referred to an undefined local `Padj` cutoff. The caller's
  `pvalue` and `padj` are now applied as described under "Significance
  cutoffs".
- `compareGSEA()` builds its NES, `padj` and `pval` tables keyed by group and
  pathway. In 0.1.x, groups whose retained pathways differed in set or order
  received values from the wrong pathway, and pathways retained only by a
  later group were dropped.
- `richGSEAcurve()` and `compareGSEA()` stop with "no pathway passes the
  cutoffs" when the cutoffs leave nothing to draw, instead of failing with an
  unrelated data-frame error; the stricter 0.2.0 defaults make this case more
  common. `richGSEAcurve()` and `richVolcano()` (with its aliases
  `ggvolcano()` and `ggVolcano()`) now validate `pvalue` and `padj` even when
  they do not filter rows.
- `detail(<GSEAResult>)` now returns `as.data.frame(x@result)` with a message
  that the class has no `detail` slot. In 0.1.5 the call failed with
  `no slot of name "detail"` because an outdated duplicate method read a slot
  that `GSEAResult` does not have.
- `richDot()` now honours `short = TRUE` and shortens long term labels like
  the other dot and bar plots; in 0.1.x both of its methods accepted `short`
  but ignored it.
- `group_by()` on `richResult`, `GSEAResult` and `Annot` objects works again.
  With dplyr 1.2 it failed with "The `add` argument of `group_by()` was
  deprecated in dplyr 1.0.0 and is now defunct"; the methods now pass `.add`
  to dplyr and keep their own `add` argument.
- `enrich()` on an `Annot` object now honours `padj.family` (new in this
  release). The `Annot` method's arguments were out of the generic's order,
  so S4 dispatch passed `padj.family` into the method's `builtin` argument
  and the adjustment always used the default `"overlap"` family. Only
  `Annot` calls that passed `padj.family = "all"` change (directly, through
  `batchEnrich()` or in a `richSensitivity()` grid): their
  adjusted p-values, tested count and recorded `padj.family` now follow the
  `"all"` family, as the same call on an annotation data frame already did.
  `richGO()` and `richKEGG()` on `Annot` objects were not affected.

### Internal changes

These changes do not alter results or the set of exported functions (the
four new exports of 0.2.0 are described under "Tested terms, comparing and
storing results"), and the only characterization snapshot they changed is
the `detail()` bug fix above. (The GSEA, `compareGSEA()`, `enrich()` and
`compareResult()` snapshots were later re-pinned for the behaviour changes
described in the sections above.)

- Removed `R/misc.R`, which held outdated duplicate definitions of 32 functions
  and 37 S4 methods or coercions that silently overrode, or were overridden by,
  those in `R/accessors.R`, `R/kappa.R`, `R/utils.R` and `R/species.R`. Each
  function and method is now defined once; the four definitions found only in
  `R/misc.R` (`reverseList_bk()`, `.clean.char()`, `.calGSEA()` and the
  `data.frame` to `richResult` coercion) moved to `R/utils.R`, `R/richGSEA.R`
  and `R/accessors.R`.
- `DESCRIPTION` now lists an explicit `Collate` field, with the class and
  generic definitions first.
- Added characterization tests (testthat edition 3 snapshots on small
  fixtures) that pin the ORA, GSEA, `richLevel`, `compareResult`,
  `compareGSEA`, `richCluster`, accessor and coercion behaviour of 0.1.5.

## Version 0.1.6 (2026-11-05)

Bug-fix release on top of 0.1.5; no changes to statistics, cutoffs or result columns.

- `richDot(short = TRUE)` now shortens long term labels; the option was
  previously ignored.
- With `builtin = FALSE`, `richKEGG()` no longer stops when the live KEGG
  pathway listing cannot be reached. It uses the bundled KEGG labels and warns,
  and it now also warns when it falls back because the live listing is
  incomplete.
- ORA plots draw two terms that share a label as two categories,
  `label [ID]`. Before, `richBar()` and `richDot()` with `order = TRUE`,
  `richNetwork()`, `richNetplot()` and `richHeatmap()` stopped with an error,
  and the other plots merged the two terms into one. `richHeatmap()` now
  matches rows across results by term ID.
- `compareGSEA()` with a `padj` cutoff no longer stops with
  "object 'Padj' not found", and it labels each group's NES, adjusted and raw
  p-values from that group's own pathway instead of by position.

## Version 0.1.5 (2026-10-20)

- Bundled KEGG and pathway data load from the richR namespace into local
  environments, including in namespace-only use.
- Kappa clustering preserves named term groups for equal-sized clusters and
  keeps representatives, memberships, and scores aligned.

## Version 0.1.4 (2026-10-20)

### Core enrichment reliability

- ORA deduplicates query genes, supports an explicit effective background,
  normalizes invalid or duplicate annotation memberships, and records analysis
  counts, options, and provider provenance in enrichment diagnostics.
- ORA adjusts p-values over the positive-overlap tested family before report
  filters; `keepRich` uses the documented reporting predicate.
- GSEA adjusts over the full tested family before significance filtering and
  exposes backend controls including `eps`, `scoreType`, `nproc`,
  `BPPARAM`, and `seed`. Seeded calls restore the caller's RNG state.

### Batch enrichment and comparisons

- `compareResult()` validates group identities, preserves real result rows and
  extra columns, and can retain empty group labels without adding placeholder
  terms. Its strict significance cutoffs are unchanged.
- `richCompareDot(..., include.all = TRUE)` keeps all requested group labels
  on the x-axis, including when pathway filtering leaves no points.
- `batchEnrich()` combines the rows returned by each group without applying a
  second cutoff and reports each group as `success`, `empty`, or `failed` in a
  typed `status` data frame.

### Packaging and namespace

- Bundled KEGG data load into local environments, preserving caller's
  `.GlobalEnv` bindings.
- Function imports are explicit; unused `reshape2` and `tidyr`
  dependencies were removed while required `sna` remains.
- Built source archives include `NEWS.md` and exclude repository-only
  files.

## Historical Version 0.1.4 release-note draft (2026-04-16; publication unverified)

The following prior release-note text is retained for context; it does not
verify that a public 0.1.4 release was published.

First hurlab-led release, shipping the `richUpset()` alignment fix and
multi-set color blending described in the 0.1.3 section below, plus
fork-specific README/figures/packaging. Proposed upstream as
[guokai8/richR PR #9](https://github.com/guokai8/richR/pull/9)
(pending review). The fork leads upstream at 0.1.4 until PR #9 merges
and upstream cuts an official release.

Functional changes over upstream:

- All of the `richUpset()` bug fix and enhancement in 0.1.3 below.

Packaging changes (fork-only):

- Restored example figures in `man/figures/` (bar, dot, GSEA,
  compare-dot, UpSet) and embedded them in the README Visualization
  section so each plot family has a rendered example.
- Rewrote the README UpSet section around the updated `richUpset()`
  behavior (pixel-accurate alignment + blended multi-set bar colors).
- README install URL, version badge, and citation updated to
  `hurlab/richR`.
- Declared `sna` as a package Import. `GGally::ggnet2()` (used by
  `R/ggrich.R`) requires `sna` at runtime, which caused the upstream
  CI vignette build to fail. This is also the fix bundled in PR #9.

---

## Version 0.1.3 (2026-04-16)

### Bug Fixes

- **`richUpset()` bar-to-dot alignment**: The intersection bar chart
  and the dot matrix previously used a nested `plot_grid()` layout
  where each row laid out the left gutter and the main plot
  independently. Differences in y-axis label widths (bar-chart tick
  labels vs. set-name labels) produced distinct panel start positions,
  causing the bars to drift out of alignment with the dots beneath
  them. Panels are now pre-aligned with `cowplot::align_plots()` and
  placed at explicit coordinates with `ggdraw()` + `draw_plot()` so
  the bar chart and matrix always share identical x-ranges regardless
  of label width. The `set.size.show = FALSE` branch additionally
  gained `align = "v", axis = "lr"` to `plot_grid()` for the same
  reason.

### Enhancements

- **`richUpset()` multi-set bar colors**: Previously, bars representing
  multi-set intersections inherited the color of the largest
  participating set. Bars for different intersections could therefore
  share a color, making it hard to tell visually which sets were
  combined. Multi-set bars now use an RGB-blend of the participating
  set colors (computed by the new internal helper `.blend_colors()`).
  Single-set bars still use the set color, and setting
  `main.bar.color` still overrides every bar to a uniform color.

### Internal

- New internal helper `.blend_colors()` in `R/ggupset.R` that averages
  RGB channels of a vector of color specifications.
- NAMESPACE: added `importFrom(cowplot, align_plots)`,
  `importFrom(cowplot, ggdraw)`, `importFrom(cowplot, draw_plot)`,
  `importFrom(grDevices, col2rgb)`, `importFrom(grDevices, rgb)`.

---

## Version 0.1.2 (2026-03-28)

### New Features
- **`show()` methods** for richResult, GSEAResult, and Annot classes: typing an
  object at the console now prints a concise summary (organism, ontology, gene
  count, top terms) instead of raw S4 slot dump.
- **Input validation**: all enrichment functions now validate gene input (empty
  lists, NAs, duplicates), check parameter ranges (minGSSize < maxGSSize), and
  report gene-annotation mapping rates. GSEA validates for named numeric input
  and removes NA/Inf values.
- **Empty result handling**: enrichment functions now return informative messages
  when no significant terms are found, instead of cryptic errors.
- **`readGMT()`**: import gene sets from GMT files (the community standard format
  used by MSigDB, Enrichr, GSEA) into Annot objects.
- **`buildAnnotFromList()`**: convert named lists of gene vectors into Annot
  objects for custom enrichment analysis.
- **`getGenes()`**: extract genes from a specific enrichment term, properly
  parsing the comma-separated GeneID field.
- **`theme_richR()`**: consistent publication-ready ggplot2 theme for all richR
  visualizations.
- **`batchEnrich()`**: run enrichment on multiple gene lists in one call, with
  progress reporting and combined output ready for `comparedot()`.
- **Unit test infrastructure**: added testthat tests for S4 class construction,
  show methods, input validation, accessor methods, and edge cases.

### Bug Fixes
- **`data()` anti-pattern**: replaced bare `data(kegg)`, `data(path)`,
  `data(module)` calls inside functions with `data(..., envir)` to avoid
  polluting the global environment.

### Internal
- New functions (`readGMT`, `buildAnnotFromList`, `getGenes`, `theme_richR`,
  `batchEnrich`) are explicitly exported.
- Added verbose `message()` calls to `buildAnnot()` reporting progress and
  database statistics.

---

## Version 0.1.1 (2026-03-28)

### New Features
- **`ggupset()`**: Color-enhanced UpSet plot for visualizing set intersections
  across multiple gene lists. Built entirely with ggplot2 (no dependency on the
  UpSetR package). Inspired by the colored UpSet implementation from the
  katelynhur/covid19VAE project.
  - Per-set coloring of bars, matrix dots, and set labels via `mycol` parameter
  - Customizable intersection line color via `line.color` 
  - Toggle per-set vs. uniform dot coloring via `color.by.set` 
  - Order intersections by frequency or degree
  - Accepts named lists of genes, richResult objects, or GSEAResult objects
  - Bar labels colored to match their respective bars

### Documentation
- Added UpSet plot section to README.md with examples and visualization image
- Version bump: 0.1.0 -> 0.1.1

---

## Version 0.1.0 (2026-03-28)

This release includes a comprehensive codebase review, bug fixes, documentation
improvements, and major refactoring to improve maintainability.
