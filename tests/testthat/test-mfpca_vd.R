make_mfpca_data <- function(seed_shift, N, maxcols) {
  set.seed(100 + seed_shift)
  D  <- matrix(NA_real_, N, maxcols)
  Ti <- matrix(NA_real_, N, maxcols)
  for (i in seq_len(N)) {
    ni <- sample(10:maxcols, 1)
    tt <- sort(stats::runif(ni, 0, ni))
    D[i, seq_len(ni)]  <- stats::rnorm(1) * sin(tt) + stats::rnorm(1) * cos(tt / 2) +
      stats::rnorm(ni, 0, 0.1)
    Ti[i, seq_len(ni)] <- tt
  }
  list(D = D, Ti = Ti)
}

test_that("mfpca_vd handles the Times input correctly", {
  N <- 30
  v1 <- make_mfpca_data(1, N, 20)
  v2 <- make_mfpca_data(2, N, 20)

  res <- mfpca_vd(
    Data  = list(v1$D, v2$D),
    Times = list(v1$Ti, v2$Ti),
    m_npcs = 2, u_npcs = 3, k_m = 6
  )

  expect_s3_class(res, "mfpca_vd")
  expect_named(res, c("scores_m", "efunctions_m", "efunctions_u", "scores_u",
                      "evalues_u", "evalues_m", "var_u", "mean_model",
                      "M_grid", "argvals_u"))
  expect_equal(nrow(res$scores_m), N)
  expect_length(res$M_grid, N)
  expect_length(res$efunctions_u, 2)
  expect_length(res$mean_model, 2)
})

test_that("mfpca_vd requires exactly one of Times or M_grid", {
  D <- list(matrix(stats::rnorm(20), 5, 4), matrix(stats::rnorm(20), 5, 4))
  expect_error(mfpca_vd(Data = D))
  expect_error(mfpca_vd(Data = D, Times = D, M_grid = 1:5))
})
