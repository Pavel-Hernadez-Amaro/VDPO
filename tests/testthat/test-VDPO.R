test_that("the 'VDPO' function with one 'ffvd' term works as expected", {
  data <- VDPO_example_vd
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))
  res <- VDPO(formula = formula, data = data)

  expect_equal(res$theta_ffvd[33], 0.22216830638841419954)
  expect_equal(res$theta_ffvd[7],  0.35787233214991615027)
})


test_that("the 'VDPO' giving the formula as a character vector works as expected", {
  data <- VDPO_example_vd
  formula <- "y ~ ffvd(X_se, nbasis = c(10, 10, 10))"
  res <- VDPO(formula = formula, data = data)

  expect_equal(res$theta_ffvd[33], 0.22216830638841419954)
  expect_equal(res$theta_ffvd[7],  0.35787233214991615027)
})

test_that("the 'VDPO' function with one 'ffpo' term works as expected", {
  data <- VDPO_example_po
  formula <- y ~ ffpo(X = x, grid = grid)
  res <- VDPO(formula = formula, data = data, family = stats::binomial())

  expect_equal(res$theta_ffpo[5], 110.89449899223717466157)
  expect_equal(res$theta_ffpo[7], 85.41612803392936825730)
})


test_that("the 'VDPO' function throws an error when a non-dataframe object is used as data", {
  data <- stats::rnorm(100)
  formula <- y ~ ffpo(X = doesnt_matter)
  expect_error(
    VDPO(formula = formula, data = data)
  )
})

test_that("the 'VDPO' function throws an error when the formula doesn't contains
          any 'ffvd', 'ffpo' or 'ffpo_2d' special", {
  formula <- y ~ SOP::f(X_se, nseg = 5)
  expect_error(
    VDPO(formula = formula, data = VDPO_example_po)
  )
})

test_that("the 'summary' method for 'VDPO' objects is identical to the 'summary'
          method for the fitted model inside the 'VDPO'n object", {
  data <- VDPO::VDPO_example_vd
  formula <- y ~ ffvd(X_se, nbasis = c(5, 5, 5))
  res <- VDPO(formula = formula, data = data)

  expect_identical(summary(res), summary(res$fit))
})




