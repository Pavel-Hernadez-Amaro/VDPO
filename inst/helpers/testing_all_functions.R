data <- data_generator_vd(beta_index = 2)
formula <- y ~ ffvd(X_se)
res <- VDPO(formula = formula, data = data)

# this raises an error
data <- data_generator_po()
formula <- y ~ ffpo(X = x, grid = grid)
res <- VDPO(formula = formula, data = data, family = stats::binomial())

data <- data_generator_po()
formula <- y ~ ffpo(X = data$x, grid = data$grid)
res <- VDPO(formula = formula, data = data, family = stats::binomial())

data <- data_generator_po2d(N = 5, px = 20, py = 21)
# res11 <- ffpo_2d(X_miss = data$X_miss, miss_points = data$miss_points, missing_points = data$missing_points)

# sapply(data, length)

data2 <- data.frame(y = data$y)
data2[["nu"]] <-  data$nu
data2[["X_true"]] <-  data$X_true
data2[["X_real"]] <-  data$X_real
data2[["X_miss"]] <-  data$X_miss
data2[["miss_points"]] <-  data$miss_points
data2[["missing_points"]] <-  data$missing_points

formula <- y ~ ffpo_2d(X_miss = data2$X_miss, miss_points = data2$miss_points, missing_points = data2$missing_points)
res <- VDPO(formula = formula, data = data2, family = stats::gaussian())
