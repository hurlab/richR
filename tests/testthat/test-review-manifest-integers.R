# SPEC-RICHR-015: reject lossy integer metadata through the public folder API.
.mi_folder <- function(x = .rp_fq()$q1) {
  p <- tempfile("manifest-integers-")
  writeRich(x, p)
  p
}
.mi_edit <- function(p, edit) {
  f <- file.path(p, "manifest.json")
  m <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  writeLines(jsonlite::toJSON(edit(m), auto_unbox = TRUE, null = "null",
                              digits = I(17)), f)
}
.mi_refuse <- function(p, field) {
  warnings <- character()
  e <- withCallingHandlers(tryCatch(readRich(p), error = identity),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(warnings, 0L)
  expect_s3_class(e, "error")
  if (inherits(e, "error")) {
    expect_match(conditionMessage(e), "manifest.json", fixed = TRUE)
    expect_match(conditionMessage(e), p, fixed = TRUE)
    expect_match(conditionMessage(e), field)
  }
}

test_that("manifest dimensions refuse fractions, overflow and nonfinite values", {
  for (field in c("nrow", "ncol")) {
    for (v in c(-1, 0.5, 1e200)) {
      p <- .mi_folder()
      .mi_edit(p, function(m) {
        i <- which(vapply(m$tables, function(e) e$slot == "detail", logical(1)))
        m$tables[[i]][[field]] <- v
        if (field == "nrow") m$tables[[i]]$row_names <- list(type = "automatic", n = -v)
        m
      })
      .mi_refuse(p, "detail")
      unlink(p, recursive = TRUE)
    }
  }
  p <- .mi_folder()
  on.exit(unlink(p, recursive = TRUE))
  .mi_edit(p, function(m) {
    i <- which(vapply(m$tables, function(e) e$slot == "detail", logical(1)))
    m$tables[[i]]$nrow <- 1e200
    m$tables[[i]]$row_names <- list(type = "automatic", n = -1e200)
    m
  })
  f <- file.path(p, "manifest.json")
  writeLines(gsub("1e+200", "1e999", readLines(f), fixed = TRUE), f)
  .mi_refuse(p, "detail")
})

test_that("integer and automatic row names refuse lossy values", {
  for (rn in list(list(type = "integer", values = c(1, 2.5, 3)),
                  list(type = "integer", values = c(1, 1e200, 3)),
                  list(type = "automatic", n = -3.5))) {
    x <- .rp_fq()$q1
    x@detail <- data.frame(x = 1:3)
    p <- .mi_folder(x)
    .mi_edit(p, function(m) {
      i <- which(vapply(m$tables, function(e) e$slot == "detail", logical(1)))
      m$tables[[i]]$row_names <- rn
      m
    })
    .mi_refuse(p, "row names.*detail")
    unlink(p, recursive = TRUE)
  }
})

test_that("typed integer metadata is checked before reconstruction", {
  for (v in list(2.5, 1e200, "2.5", "garbage")) {
    p <- .mi_folder()
    .mi_edit(p, function(m) {
      i <- match("schema_version", unlist(m$metadata$attributes$names$value))
      m$metadata$value[[i]] <- list(type = "integer", value = list(v))
      m
    })
    .mi_refuse(p, "metadata.*schema_version")
    unlink(p, recursive = TRUE)
  }
})

test_that("legitimate dimensions, row names and missing integers round trip", {
  for (rn in list(NULL, c(3L, 1L, 2L), c("x", "x", "z"))) {
    for (nc in 0:1) {
      x <- .rp_fq()$q1
      tab <- if (nc) data.frame(x = 1:3) else data.frame(row.names = 1:3)
      if (!is.null(rn)) attr(tab, "row.names") <- rn
      x@detail <- tab
      meta <- attr(x, "richR_metadata")
      meta$integer_control <- c(NA_integer_, 0L, -.Machine$integer.max, .Machine$integer.max)
      attr(x, "richR_metadata") <- meta
      p <- .mi_folder(x)
      y <- readRich(p)
      expect_identical(y@detail, x@detail)
      expect_identical(attr(y, "richR_metadata"), attr(x, "richR_metadata"))
      unlink(p, recursive = TRUE)
    }
  }
  x <- .rp_fq()$q1
  x@detail <- data.frame()
  p <- .mi_folder(x)
  on.exit(unlink(p, recursive = TRUE))
  expect_identical(readRich(p)@detail, x@detail)
})
