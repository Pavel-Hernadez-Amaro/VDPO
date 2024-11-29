test_that("po_fit handles basic case correctly", {
  data_1d <- data_generator_po_1d(noise_sd = 0.1)

  formula <- response ~ ffpo(X = data_1d$noisy_curves_miss, grid = data_1d$grid, nbasis = c(10, 10))

  res_po <- po_fit(
    formula = formula,
    data = data_1d
  )

  expect_s3_class(res_po, "po_fit")
  expect_type(res_po, "list")
  expect_named(res_po, c("fit", "theta", "ffpo_evals"))

  expect_true(attr(res_po, "N") > 0)
})

test_that("po_fit handles offset", {
  data_1d <- data_generator_po_1d(noise_sd = 0.1)

  formula <- response ~ ffpo(X = data_1d$noisy_curves_miss, grid = data_1d$grid, nbasis = c(10, 10))

  offset <- rep(0.1, length(data_1d$response))

  res_po <- po_fit(
    formula = formula,
    data = data_1d,
    offset = offset
  )

  expect_s3_class(res_po, "po_fit")
})

# test_that("po_fit handles different family distributions", {
#   data_1d <- data_generator_po_1d(noise_sd = 0.1)
#
#   formula <- response ~ ffpo(X = data_1d$noisy_curves_miss, grid = data_1d$grid, nbasis = c(10, 10))
#
#   res_po_poisson <- po_fit(
#     formula = formula,
#     data = data_1d,
#     family = stats::poisson()
#   )
#
#   expect_s3_class(res_po_poisson, "po_fit")
# })


test_that("po_fit throws appropriate errors", {
  # Generate test data
  data_1d <- data_generator_po_1d(noise_sd = 0.1)

  # Invalid formula without ffpo
  expect_error(
    po_fit(response ~ noisy_curves_miss, data = data_1d),
    regexp = "this function should be used with at least one 'ffpo' term"
  )

  # Invalid data type
  expect_error(
    po_fit(response ~ ffpo(X = data_1d$noisy_curves_miss, grid = data_1d$grid),
           data = c()),
    regexp = "the 'data' argument should be a list"
  )

  # Insufficient data
  expect_error(
    po_fit(response ~ ffpo(X = data_1d$noisy_curves_miss, grid = data_1d$grid),
           data = list(response = c(1))),
    regexp = "'data' should have at least a response and a covariate"
  )
})
