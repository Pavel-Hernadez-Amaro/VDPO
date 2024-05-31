test_that("the number of intervals should be an even number", {
  args <- partial_inprod_arguments_generator()

  # expect lack of error
  expect_error(do.call(partial_inprod, args), NA)

  # expect error
  args$n_intervals <- args$n_intervals - 1
  expect_error(do.call(partial_inprod, args))
})


test_that("'rng' should be a vector with two elements", {
  args <- partial_inprod_arguments_generator()

  # expect lack of error
  expect_error(do.call(partial_inprod, args), NA)

  # expect error
  args$rng <- args$rng[1]
  expect_error(do.call(partial_inprod, args))
})

test_that("'bdeg' should be a vector with two elements", {
  args <- partial_inprod_arguments_generator()

  # expect lack of error
  expect_error(do.call(partial_inprod, args), NA)

  # expect error
  args$bdeg<- args$bdeg[1]
  expect_error(do.call(partial_inprod, args))
})

test_that("'partial_inprod' with an odd number of inervals", {
  # Setup
  n_intervals <- 3
  knots1 <- c(1, 2)
  knots2 <- c(3, 4)
  bdeg <- c(1, 1)
  spline_domain <- matrix(1:4, nrow = 2)
  rng <- c(0, 1)

  # Call the function and expect an error
  expect_error(partial_inprod(n_intervals, knots1, knots2, bdeg, spline_domain, rng),
               "the 'n_intervals' parameter should be an even number")
})
