test_that("Test for bspline function", {
  x <- c(1, 2, 3, 4)
  xl <- 1
  xr <- 4
  nseg <- 3
  bdeg <- 2

  result <- bspline(x, xl, xr, nseg, bdeg)

  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(is.matrix(result[[1]]))
})
