context("pleiotropy")

example_file <- system.file(file.path("extdata", "example-CAMERA.rds"), package = "CAMeRa")
xold <- readRDS(example_file)
x <- CAMERA$new()
x$import(xold)
rm(xold)
test_that("pleiotropy", {
  x$harmonise()
  x$cross_estimate()
  x$pleiotropy()
  x$pleiotropy_outliers
  x$pleiotropy_Q_outliers
  x$pleiotropy_agreement
  x$plot_pleiotropy_heterogeneity(pthresh=1)
  expect_true(!is.null(x$pleiotropy_outliers))
})
