# names/colnames/rownames/row.names/dim are pinned for all three classes

    list(richResult = list(names = c("Annot", "Term", "Annotated", 
    "Significant", "RichFactor", "FoldEnrichment", "zscore", "Pvalue", 
    "Padj", "GeneID"), colnames = c("Annot", "Term", "Annotated", 
    "Significant", "RichFactor", "FoldEnrichment", "zscore", "Pvalue", 
    "Padj", "GeneID"), rownames = c("GO:0001", "GO:0002", "GO:0003"
    ), row.names = c("GO:0001", "GO:0002", "GO:0003"), dim = c(3L, 
    10L)), GSEAResult = list(names = c("pathway", "pval", "padj", 
    "log2err", "ES", "NES", "size", "leadingEdge"), colnames = c("pathway", 
    "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge"
    ), rownames = c("1", "2", "3"), row.names = c("1", "2", "3"), 
        dim = c(3L, 8L)), Annot = list(names = c("GeneID", "Term", 
    "Annot"), colnames = c("GeneID", "Term", "Annot"), rownames = c("1", 
    "2", "3", "4"), row.names = c("1", "2", "3", "4"), dim = 4:3))

# head() and tail() print their banner and return the leading/trailing rows

    Code
      cat("class:", cls, "\n")
    Output
      class: richResult 
    Code
      print(head(x, 2))
    Output
      === Total significant terms is: 3 10 ===
                Annot       Term Annotated Significant RichFactor FoldEnrichment
      GO:0001 GO:0001  apoptosis       100          10        0.1            2.5
      GO:0002 GO:0002 cell cycle       200          20        0.1            2.0
              zscore Pvalue Padj  GeneID
      GO:0001    3.5  0.001 0.01   A,B,C
      GO:0002    3.0  0.010 0.05 D,E,F,G
    Code
      print(tail(x, 2))
    Output
      === Total significant terms is: 3 10 ===
                Annot       Term Annotated Significant RichFactor FoldEnrichment
      GO:0002 GO:0002 cell cycle       200          20        0.1            2.0
      GO:0003 GO:0003  signaling        50           5        0.1            1.5
              zscore Pvalue Padj  GeneID
      GO:0002    3.0   0.01 0.05 D,E,F,G
      GO:0003    2.5   0.04 0.08     A,H

---

    Code
      cat("class:", cls, "\n")
    Output
      class: GSEAResult 
    Code
      print(head(x, 2))
    Output
      === Total significant terms is: 3 8 ===
        pathway  pval  padj log2err   ES  NES size leadingEdge
      1      p1 0.001 0.003     0.5  0.8  2.1    5         a,b
      2      p2 0.010 0.015     0.4 -0.7 -1.9    8           c
    Code
      print(tail(x, 2))
    Output
      === Total significant terms is: 3 8 ===
        pathway pval  padj log2err   ES  NES size leadingEdge
      2      p2 0.01 0.015     0.4 -0.7 -1.9    8           c
      3      p3 0.20 0.200     0.1  0.3  0.9    7       d,e,f

---

    Code
      cat("class:", cls, "\n")
    Output
      class: Annot 
    Code
      print(head(x, 2))
    Output
      === species is: human and Annotation is GO keytype is SYMBOL ===
        GeneID    Term     Annot
      1      A GO:0001 apoptosis
      2      B GO:0001 apoptosis
    Code
      print(tail(x, 2))
    Output
      === species is: human and Annotation is GO keytype is SYMBOL ===
        GeneID    Term      Annot
      3      C GO:0002 cell cycle
      4      D GO:0002 cell cycle

# [ and $ subset the underlying data.frame for all three classes

    list(richResult = list(rows_cols = structure(list(Annot = c("GO:0001", 
    "GO:0002"), Term = c("apoptosis", "cell cycle")), row.names = c("GO:0001", 
    "GO:0002"), class = "data.frame"), first_col = c("GO:0001", "GO:0002", 
    "GO:0003"), single_cell = "apoptosis", one_row_drop_false = structure(list(
        Annot = "GO:0001", Term = "apoptosis"), row.names = "GO:0001", class = "data.frame"), 
        dollar = c("GO:0001", "GO:0002", "GO:0003")), GSEAResult = list(
        rows_cols = structure(list(pathway = c("p1", "p2"), pval = c(0.001, 
        0.01)), row.names = c("1", "2"), class = "data.frame"), first_col = c("p1", 
        "p2", "p3"), single_cell = 0.001, one_row_drop_false = structure(list(
            pathway = "p1", pval = 0.001), row.names = "1", class = "data.frame"), 
        dollar = NULL), Annot = list(rows_cols = structure(list(GeneID = c("A", 
    "B"), Term = c("GO:0001", "GO:0001")), row.names = 1:2, class = "data.frame"), 
        first_col = c("A", "B", "C", "D"), single_cell = "GO:0001", 
        one_row_drop_false = structure(list(GeneID = "A", Term = "GO:0001"), row.names = 1L, class = "data.frame"), 
        dollar = c("apoptosis", "apoptosis", "cell cycle", "cell cycle"
        )))

# result() and summary() are pinned for richResult and GSEAResult

    list(rich = structure(list(Annot = c("GO:0001", "GO:0002", "GO:0003"
    ), Term = c("apoptosis", "cell cycle", "signaling"), Annotated = c(100, 
    200, 50), Significant = c(10, 20, 5), RichFactor = c(0.1, 0.1, 
    0.1), FoldEnrichment = c(2.5, 2, 1.5), zscore = c(3.5, 3, 2.5
    ), Pvalue = c(0.001, 0.01, 0.04), Padj = c(0.01, 0.05, 0.08), 
        GeneID = c("A,B,C", "D,E,F,G", "A,H")), class = "data.frame", row.names = c("GO:0001", 
    "GO:0002", "GO:0003")), gsea = structure(list(pathway = c("p1", 
    "p2", "p3"), pval = c(0.001, 0.01, 0.2), padj = c(0.003, 0.015, 
    0.2), log2err = c(0.5, 0.4, 0.1), ES = c(0.8, -0.7, 0.3), NES = c(2.1, 
    -1.9, 0.9), size = c(5L, 8L, 7L), leadingEdge = c("a,b", "c", 
    "d,e,f")), class = "data.frame", row.names = c("1", "2", "3")))

---

    Code
      summary(.char_richResult())
    Output
      Total input genes is: 10 and significant biological term is: 3 
    Code
      summary(.char_gsea_result())
    Output
      Total significant biological term is: 2 

# detail() is pinned for richResult and for GSEAResult (current failure mode)

    structure(list(TERM = c("GO:0001", "GO:0001", "GO:0002"), Annot = c("apoptosis", 
    "apoptosis", "cell cycle"), GeneID = c("A", "B", "D"), Pvalue = c(0.001, 
    0.001, 0.01), Padj = c(0.01, 0.01, 0.05)), class = "data.frame", row.names = c(NA, 
    -3L))

---

    Code
      detail(.char_gsea_result())
    Condition
      Error in `as.data.frame()`:
      ! no slot of name "detail" for this object of class "GSEAResult"

# as.data.frame() and as(<class>, 'data.frame') pin the three coercions

    list(richResult = list(as.data.frame = structure(list(Annot = c("GO:0001", 
    "GO:0002", "GO:0003"), Term = c("apoptosis", "cell cycle", "signaling"
    ), Annotated = c(100, 200, 50), Significant = c(10, 20, 5), RichFactor = c(0.1, 
    0.1, 0.1), FoldEnrichment = c(2.5, 2, 1.5), zscore = c(3.5, 3, 
    2.5), Pvalue = c(0.001, 0.01, 0.04), Padj = c(0.01, 0.05, 0.08
    ), GeneID = c("A,B,C", "D,E,F,G", "A,H")), class = "data.frame", row.names = c("GO:0001", 
    "GO:0002", "GO:0003")), as = structure(list(Annot = c("GO:0001", 
    "GO:0002", "GO:0003"), Term = c("apoptosis", "cell cycle", "signaling"
    ), Annotated = c(100, 200, 50), Significant = c(10, 20, 5), RichFactor = c(0.1, 
    0.1, 0.1), FoldEnrichment = c(2.5, 2, 1.5), zscore = c(3.5, 3, 
    2.5), Pvalue = c(0.001, 0.01, 0.04), Padj = c(0.01, 0.05, 0.08
    ), GeneID = c("A,B,C", "D,E,F,G", "A,H")), class = "data.frame", row.names = c("GO:0001", 
    "GO:0002", "GO:0003"))), GSEAResult = list(as.data.frame = structure(list(
        pathway = c("p1", "p2", "p3"), pval = c(0.001, 0.01, 0.2), 
        padj = c(0.003, 0.015, 0.2), log2err = c(0.5, 0.4, 0.1), 
        ES = c(0.8, -0.7, 0.3), NES = c(2.1, -1.9, 0.9), size = c(5L, 
        8L, 7L), leadingEdge = c("a,b", "c", "d,e,f")), class = "data.frame", row.names = c("1", 
    "2", "3")), as = structure(list(pathway = c("p1", "p2", "p3"), 
        pval = c(0.001, 0.01, 0.2), padj = c(0.003, 0.015, 0.2), 
        log2err = c(0.5, 0.4, 0.1), ES = c(0.8, -0.7, 0.3), NES = c(2.1, 
        -1.9, 0.9), size = c(5L, 8L, 7L), leadingEdge = c("a,b", 
        "c", "d,e,f")), class = "data.frame", row.names = c("1", 
    "2", "3"))), Annot = list(as.data.frame = structure(list(GeneID = c("A", 
    "B", "C", "D"), Term = c("GO:0001", "GO:0001", "GO:0002", "GO:0002"
    ), Annot = c("apoptosis", "apoptosis", "cell cycle", "cell cycle"
    )), class = "data.frame", row.names = c(NA, -4L)), as = structure(list(
        GeneID = c("A", "B", "C", "D"), Term = c("GO:0001", "GO:0001", 
        "GO:0002", "GO:0002"), Annot = c("apoptosis", "apoptosis", 
        "cell cycle", "cell cycle")), class = "data.frame", row.names = c(NA, 
    -4L))))

# as(<data.frame>, 'Annot') and as(<data.frame>, 'richResult') are pinned

    list(annot = new("Annot", species = character(0), anntype = character(0), 
        keytype = character(0), annot = structure(list(GeneID = c("A", 
        "B", "C"), Term = c("T1", "T1", "T2"), Annot = c("first", 
        "first", "second")), class = "data.frame", row.names = c(NA, 
        -3L))), richResult = new("richResult", result = structure(list(
        Annot = c("GO:0001", "GO:0002", "GO:0003"), Term = c("apoptosis", 
        "cell cycle", "signaling"), Annotated = c(100, 200, 50), 
        Significant = c(10, 20, 5), Pvalue = c(0.001, 0.01, 0.04), 
        Padj = c(0.01, 0.05, 0.08), GeneID = c("A,B,C", "D,E,F,G", 
        "A,H")), class = "data.frame", row.names = c("GO:0001", "GO:0002", 
    "GO:0003")), detail = structure(list(), names = character(0), row.names = integer(0), class = "data.frame"), 
        pvalueCutoff = numeric(0), pAdjustMethod = character(0), 
        padjCutoff = numeric(0), genenumber = 8L, organism = character(0), 
        ontology = character(0), gene = c("A", "B", "C", "D", "E", 
        "F", "G", "H"), keytype = character(0), sep = ","))

# rbind() on richResult and GSEAResult objects is pinned

    list(rich_distinct = structure(list(Annot = c("GO:0001", "GO:0002", 
    "GO:0003", "GO:0010", "GO:0011"), Term = c("apoptosis", "cell cycle", 
    "signaling", "metabolism", "transport"), Annotated = c(100, 200, 
    50, 80, 60), Significant = c(10, 20, 5, 8, 6), RichFactor = c(0.1, 
    0.1, 0.1, 0.1, 0.1), FoldEnrichment = c(2.5, 2, 1.5, 2.2, 1.8
    ), zscore = c(3.5, 3, 2.5, 3.2, 2.8), Pvalue = c(0.001, 0.01, 
    0.04, 0.002, 0.02), Padj = c(0.01, 0.05, 0.08, 0.02, 0.06), GeneID = c("A,B,C", 
    "D,E,F,G", "A,H", "B,I", "C,J")), row.names = c("GO:0001", "GO:0002", 
    "GO:0003", "GO:0010", "GO:0011"), class = "data.frame"), rich_duplicate_rownames = structure(list(
        Annot = c("GO:0001", "GO:0002", "GO:0003", "GO:0001", "GO:0002", 
        "GO:0003"), Term = c("apoptosis", "cell cycle", "signaling", 
        "apoptosis", "cell cycle", "signaling"), Annotated = c(100, 
        200, 50, 100, 200, 50), Significant = c(10, 20, 5, 10, 20, 
        5), RichFactor = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1), FoldEnrichment = c(2.5, 
        2, 1.5, 2.5, 2, 1.5), zscore = c(3.5, 3, 2.5, 3.5, 3, 2.5
        ), Pvalue = c(0.001, 0.01, 0.04, 0.001, 0.01, 0.04), Padj = c(0.01, 
        0.05, 0.08, 0.01, 0.05, 0.08), GeneID = c("A,B,C", "D,E,F,G", 
        "A,H", "A,B,C", "D,E,F,G", "A,H")), row.names = c("GO:0001", 
    "GO:0002", "GO:0003", "GO:00011", "GO:00021", "GO:00031"), class = "data.frame"), 
        gsea_distinct = structure(list(pathway = c("p1", "p2", "p3", 
        "p4", "p5"), pval = c(0.001, 0.01, 0.2, 0.002, 0.03), padj = c(0.003, 
        0.015, 0.2, 0.004, 0.03), log2err = c(0.5, 0.4, 0.1, 0.5, 
        0.3), ES = c(0.8, -0.7, 0.3, 0.6, -0.5), NES = c(2.1, -1.9, 
        0.9, 1.8, -1.6), size = c(5L, 8L, 7L, 6L, 9L), leadingEdge = c("a,b", 
        "c", "d,e,f", "a,c", "e,f")), row.names = c("1", "2", "3", 
        "11", "21"), class = "data.frame"), gsea_duplicate_rownames = structure(list(
            pathway = c("p1", "p2", "p3", "p1", "p2", "p3"), pval = c(0.001, 
            0.01, 0.2, 0.001, 0.01, 0.2), padj = c(0.003, 0.015, 
            0.2, 0.003, 0.015, 0.2), log2err = c(0.5, 0.4, 0.1, 0.5, 
            0.4, 0.1), ES = c(0.8, -0.7, 0.3, 0.8, -0.7, 0.3), NES = c(2.1, 
            -1.9, 0.9, 2.1, -1.9, 0.9), size = c(5L, 8L, 7L, 5L, 
            8L, 7L), leadingEdge = c("a,b", "c", "d,e,f", "a,b", 
            "c", "d,e,f")), row.names = c("1", "2", "3", "11", "21", 
        "31"), class = "data.frame"))

