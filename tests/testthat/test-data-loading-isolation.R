package_data_reference <- function(name) {
  reference <- new.env(parent = emptyenv())
  utils::data(list = name, package = "richR", envir = reference)
  reference
}

with_global_data_state <- function(values = list(), absent = character(), code) {
  names_to_restore <- unique(c(names(values), absent))
  previous <- lapply(names_to_restore, function(name) {
    present <- exists(name, envir = .GlobalEnv, inherits = FALSE)
    list(present = present,
         value = if (present) get(name, envir = .GlobalEnv, inherits = FALSE) else NULL)
  })
  names(previous) <- names_to_restore
  on.exit(for (name in names_to_restore) {
    state <- previous[[name]]
    if (state$present) {
      assign(name, state$value, envir = .GlobalEnv)
    } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }, add = TRUE)

  for (name in names(values)) assign(name, values[[name]], envir = .GlobalEnv)
  for (name in absent) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }
  force(code)
}

run_kegg_gsea_with_bundled_path <- function() {
  reference <- package_data_reference("path")
  usable <- which(!is.na(reference$path$Level3) & nzchar(reference$path$Level3))
  term <- as.character(reference$path$Level3[[usable[[1]]]])
  annotation <- data.frame(GeneID = paste0("g", 1:6), Term = rep(term, 6))
  stats <- setNames(c(-10:-1, 1:10), paste0("g", 1:20))
  result <- suppressWarnings(richR::richGSEA(
    stats, annotation, KEGG = TRUE, pvalue = 1, minSize = 1, maxSize = 20,
    nproc = 1, seed = 42
  ))
  list(result = result, reference = reference$path)
}

test_that("KEGG GSEA does not alter present or absent global path bindings", {
  marker <- data.frame(sentinel = 123L)
  with_global_data_state(values = list(path = marker), code = {
    actual <- run_kegg_gsea_with_bundled_path()
    expect_identical(get("path", envir = .GlobalEnv), marker)
    expected <- actual$reference
    rownames(expected) <- expected$Level3
    expected <- expected[actual$result@result$pathway, , drop = FALSE]
    expect_s4_class(actual$result, "GSEAResult")
    expect_identical(actual$result@result$Level3, expected$Level3)
    expect_identical(actual$result@result$Level2, expected$Level2)
  })

  with_global_data_state(absent = "path", code = {
    actual <- run_kegg_gsea_with_bundled_path()
    expect_false(exists("path", envir = .GlobalEnv, inherits = FALSE))
    expected <- actual$reference
    rownames(expected) <- expected$Level3
    expected <- expected[actual$result@result$pathway, , drop = FALSE]
    expect_identical(actual$result@result$Level3, expected$Level3)
    expect_identical(actual$result@result$Level2, expected$Level2)
  })
})

test_that("bundled KEGG and module loaders isolate existing and absent globals", {
  kegg_reference <- package_data_reference("kegg")$kegg.db
  module_reference <- package_data_reference("module")$module
  marker <- list(sentinel = TRUE)
  active_kegg <- getFromNamespace(".get_kg_dat", "richR")
  active_module <- getFromNamespace(".get_kgm.data", "richR")
  loaders <- list(kegg = active_kegg, module = active_module)

  for (globals in list(list(kegg.db = marker, module = marker), character())) {
    if (length(globals)) {
      with_global_data_state(values = globals, code = {
        for (loader_name in names(loaders)) {
          actual <- loaders[[loader_name]]()
          if (identical(loader_name, "kegg")) {
            expect_identical(actual, kegg_reference)
          } else {
            expect_identical(actual, module_reference)
          }
        }
        expect_identical(get("kegg.db", envir = .GlobalEnv), marker)
        expect_identical(get("module", envir = .GlobalEnv), marker)
      })
    } else {
      with_global_data_state(absent = c("kegg.db", "module"), code = {
        for (loader_name in names(loaders)) {
          actual <- loaders[[loader_name]]()
          if (identical(loader_name, "kegg")) {
            expect_identical(actual, kegg_reference)
          } else {
            expect_identical(actual, module_reference)
          }
        }
        expect_false(exists("kegg.db", envir = .GlobalEnv, inherits = FALSE))
        expect_false(exists("module", envir = .GlobalEnv, inherits = FALSE))
      })
    }
  }

  with_global_data_state(values = list(kegg.db = marker), code = {
    fallback_kegg <- testthat::with_mocked_bindings(
      keggList = function(...) setNames("fallback trigger", "path:hsa00001"),
      active_kegg(builtin = FALSE),
      .package = "richR"
    )
    expect_identical(fallback_kegg, kegg_reference)
    expect_identical(get("kegg.db", envir = .GlobalEnv), marker)
  })
})
