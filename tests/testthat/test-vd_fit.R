test_that("vd_fit handles basic case correctly", {
  data <- data_generator_vd(beta_index = 1)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))

  result <- vd_fit(formula = formula, data = data)

  expect_s3_class(result, "vd_fit")
  expect_type(result, "list")
  expect_named(result, c("fit", "Beta", "theta", "covar_theta", "M",
                         "ffvd_evals", "theta_no_functional", "theta_f"))

  expect_true(length(result$Beta) == 1)
  expect_true(all(sapply(result$Beta, is.matrix)))
  expect_true(is.matrix(result$covar_theta))
  expect_true(is.list(result$M))
})

test_that("vd_fit handles multiple ffvd terms", {
  data <- data_generator_vd(beta_index = 1)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) +
    ffvd(Y_se, nbasis = c(10, 10, 10))

  result <- vd_fit(formula = formula, data = data)

  expect_length(result$Beta, 2)
  expect_length(result$M, 2)
  expect_true(all(sapply(result$Beta, is.matrix)))
})

test_that("vd_fit handles non-aligned data", {
  data_not_aligned <- data_generator_vd(aligned = FALSE, beta_index = 2)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))

  result <- vd_fit(formula = formula, data = data_not_aligned)

  expect_s3_class(result, "vd_fit")
  expect_true(is.list(result$M))
  expect_true(all(sapply(result$M, is.matrix)))
})

test_that("vd_fit handles additional covariates", {
  data <- data_generator_vd(beta_index = 2, use_x = TRUE)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) + x2

  result <- vd_fit(formula = formula, data = data)

  expect_false(is.null(result$theta_no_functional))
  expect_length(result$theta_no_functional, 1)
})

test_that("vd_fit handles f() terms", {
  data <- data_generator_vd(beta_index = 2, use_x = TRUE, use_f = TRUE)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) +
    f(x1, nseg = 20, pord = 2, degree = 3) + x2

  result <- vd_fit(formula = formula, data = data)

  expect_false(is.null(result$theta_f))
  expect_false(is.null(result$theta_no_functional))
})

test_that("vd_fit handles offset", {
  data <- data_generator_vd(beta_index = 1)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))
  offset <- rep(0.1, length(data[["y"]]))

  result <- vd_fit(formula = formula, data = data, offset = offset)

  expect_s3_class(result, "vd_fit")
})

test_that("process_ffvd_evals processes evaluations correctly", {
  data <- data_generator_vd(beta_index = 1)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))
  result <- vd_fit(formula = formula, data = data)

  expect_named(result$ffvd_evals, "ffvd_X_se")
  expect_type(result$ffvd_evals, "list")
})

test_that("calculate_beta_ffvd computes Beta correctly", {
  data <- data_generator_vd(beta_index = 1)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))
  result <- vd_fit(formula = formula, data = data)

  expect_type(result$Beta, "list")
  expect_true(all(sapply(result$Beta, is.matrix)))
})

test_that("calculate_theta_aux handles different cases", {
  # Only non-special indices
  data <- data_generator_vd(beta_index = 2, use_x = TRUE)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) + x2
  result1 <- vd_fit(formula = formula, data = data)

  # Only f terms
  data <- data_generator_vd(beta_index = 2, use_x = TRUE, use_f = TRUE)
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) +
    f(x1, nseg = 30, pord = 2, degree = 3)
  expect_error(result2 <- vd_fit(formula = formula, data = data))

  # Both f terms and non-special indices
  formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) +
    f(x1, nseg = 25, pord = 2, degree = 3) + x2
  result3 <- vd_fit(formula = formula, data = data)

  expect_false(is.null(result1$theta_no_functional))
  expect_false(is.null(result3$theta_f))
  expect_false(is.null(result3$theta_no_functional))
})

# Test error cases
test_that("vd_fit throws appropriate errors", {
  data <- data_generator_vd(beta_index = 1)

  # Test invalid formula
  expect_error(vd_fit(y ~ x, data = data),
               "this function should be used with at least one 'ffvd' term")

  # Test invalid data type
  expect_error(vd_fit(y ~ ffvd(X_se), data = list()),
               "'data' should have at least a response and a covariate")
})

