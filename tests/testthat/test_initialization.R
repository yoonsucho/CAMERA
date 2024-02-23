context("CAMERA")
library(CAMERA)

test_that("initialization", {
  x <- CAMERA$new(
    exposure_ids=c("ieu-a-2", "bbj-a-1"),
    outcome_ids=c("ieu-a-7", "bbj-a-109"),
    pops = c("EUR", "EAS")
  )
  eids=c("ieu-a-2", "bbj-a-1")
  oids=c("ieu-a-7", "bbj-a-109")
  expect_identical(length(eids), length(x$exposure_ids))
  expect_identical(length(oids), length(x$outcome_ids))
})

test_that("extract_instrument", {
  x <- CAMERA$new(
    exposure_ids=c("ieu-a-2", "bbj-a-1"),
    outcome_ids=c("ieu-a-7", "bbj-a-109"),
    pops = c("EUR", "EAS")
  )
  x$extract_instruments()
  expect_true(nrow(x$instrument_raw) > 2)
})


test_that("import", {
  example_file <- system.file(file.path("extdata", "example-CAMERA.rds"), package = "CAMeRa")
  xold <- readRDS(example_file)
  x <- CAMERA$new()
  x$import(xold)
  rm(xold)
  expect_true(inherits(x, "CAMERA"))
})
