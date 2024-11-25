test_that("po_2d_fit handles basic case correctly", {
  # data_2d_internal <- generate_2d_po_functional_data(n = 10, noise_sd = 0.1)
  formula <- response ~ ffpo_2d(
    X = data_2d_internal$noisy_surfaces_miss,
    miss_points = data_2d_internal$miss_points,
    missing_points = data_2d_internal$missing_points
  )

  res_po_2d <- po_2d_fit(
    formula = formula,
    data = data_2d_internal
  )

  expect_s3_class(res_po_2d, "po_2d_fit")
  expect_type(res_po_2d, "list")
  expect_named(res_po_2d, c("fit", "theta", "ffpo_2d_evals"))

  expect_true(attr(res_po_2d, "N") > 0)
})

test_that("po_2d_fit handles offset correctly", {
  formula <- response ~ ffpo_2d(
    X = data_2d_internal$noisy_surfaces_miss,
    miss_points = data_2d_internal$miss_points,
    missing_points = data_2d_internal$missing_points
  )

  offset <- rep(0.1, length(data_2d_internal$response))

  res_po_2d <- po_2d_fit(
    formula = formula,
    data = data_2d_internal,
    offset = offset
  )

  expect_s3_class(res_po_2d, "po_2d_fit")
})

test_that("po_2d_fit throws appropriate errors for invalid inputs", {
  expect_error(
    po_2d_fit(response ~ noisy_surfaces_miss, data = data_2d_internal),
    regexp = "this function should be used with at least one 'ffpo_2d' term"
  )

  expect_error(
    po_2d_fit(response ~ ffpo_2d(
      X = data_2d_internal$noisy_surfaces_miss,
      miss_points = data_2d_internal$miss_points,
      missing_points = data_2d_internal$missing_points),
      data = c()),
    regexp = "the 'data' argument should be a list"
  )

  expect_error(
    po_2d_fit(response ~ ffpo_2d(
      X = data_2d_internal$noisy_surfaces_miss,
      miss_points = data_2d_internal$miss_points,
      missing_points = data_2d_internal$missing_points),
      data = list(response = c(1))),
    regexp = "'data' should have at least a response and a covariate"
  )
})

# test_that("po_2d_fit handles different coefficient surfaces", {
#   data_2d_saddle <- generate_2d_po_functional_data(beta_type = "saddle", noise_sd = 0.1)
#   data_2d_exp <- generate_2d_po_functional_data(beta_type = "exp", noise_sd = 0.1)
#   formula <- response ~ ffpo_2d(
#     X = data_2d_saddle$noisy_surfaces_miss,
#     miss_points = data_2d_saddle$miss_points,
#     missing_points = data_2d_saddle$missing_points
#   )
#
#   res_po_2d_saddle <- po_2d_fit(
#     formula = formula,
#     data = data_2d_saddle
#   )
#
#   res_po_2d_exp <- po_2d_fit(
#     formula = formula,
#     data = data_2d_exp
#   )
#
#   expect_s3_class(res_po_2d_saddle, "po_2d_fit")
#   expect_s3_class(res_po_2d_exp, "po_2d_fit")
#
#   expect_false(all.equal(res_po_2d_saddle$theta, res_po_2d_exp$theta))
# })
