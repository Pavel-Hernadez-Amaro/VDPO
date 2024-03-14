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

# the main functionality of the function is already tested in 'VDPO'
