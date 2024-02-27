test_that("the 'VDPO' function with one 'ffvd' term works as expected", {
  data <- VDPO::VDPO_example_vd
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))
  res <- VDPO(formula = formula, data = data)

  expect_equal(res$theta_ffvd[33], 0.22216830638841419954)
  expect_equal(res$theta_ffvd[7],  0.35787233214991615027)
})


test_that("the 'VDPO' function with one 'ffpo' term works as expected", {
  data <- VDPO::VDPO_example_po
  formula <- y ~ ffpo(X = x, grid = grid)
  res <- VDPO(formula = formula, data = data, family = stats::binomial())

  expect_equal(res$theta_ffpo[5], 110.89449899223717466157)
  expect_equal(res$theta_ffpo[7], 85.41612803392936825730)
})

