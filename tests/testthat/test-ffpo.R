test_that("the 'X' argument should be a matrix", {
  data <- VDPO_example_po
  x <- c(data$x)

  expect_error(ffpo(x, data$grid))
})

test_that("the add_grid function works", {
  data <- VDPO_example_po
  grid <- c(data$grid)# transform it into an atomic vector
  data_without_grid <- subset(data, select = -grid)

  expect_identical(add_grid(data_without_grid, grid), data)
})
