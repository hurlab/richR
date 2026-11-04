# richDAVID() deprecation (SPEC-RICHR-002 REQ-025 / AC-025). No network is used:
# the call stops at the RDAVIDWebService check, after the deprecation warning.

test_that("richDAVID() warns that it is deprecated before the package-missing error", {
  skip_if(requireNamespace("RDAVIDWebService", quietly = TRUE),
          "RDAVIDWebService is installed; the package-missing path cannot be reached")
  expect_warning(
    expect_error(richDAVID(c("1", "2")), "RDAVIDWebService"),
    "deprecated"
  )
  expect_warning(
    expect_error(richDAVID(c("1", "2")), "RDAVIDWebService"),
    "no longer available from Bioconductor.*next minor release"
  )
})

test_that("the deprecation warning comes before any other check", {
  expect_warning(
    expect_error(richDAVID(c("1", "2"), pvalue = 0), "must be between 0 and 1"),
    "deprecated"
  )
})
