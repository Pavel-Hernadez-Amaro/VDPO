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
data2[["nu"]] <- data$nu
data2[["X_true"]] <- data$X_true
data2[["X_real"]] <- data$X_real
data2[["X_miss"]] <- data$X_miss
data2[["miss_points"]] <-  data$miss_points
data2[["missing_points"]] <-  data$missing_points

formula <- y ~ ffpo_2d(X_miss = data2$X_miss, miss_points = data2$miss_points, missing_points = data2$missing_points)
res <- VDPO(formula = formula, data = data2, family = stats::gaussian())


# 30092024
data <- data_generator_vd(beta_index = 2)
formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) + ffvd(Y_se, nbasis = c(11, 10, 10))
res <- vd_fit(formula = formula, data = data)

formula2 <- y ~ ffvd_old(X_se, nbasis = c(30, 50, 30))
t0 <- proc.time()
res2 <- VDPO(formula = formula2, data = data)
t1 <- proc.time()

time2 <- t1 - t0


# data <- data_generator_vd(beta_index = 1)
formula <- y ~ ffvd(X_se, nbasis = c(30, 50, 30))
t0 <- proc.time()
res3 <- vd_fit(formula = formula, data = data)
t1 <- proc.time()

time3 <- t1 - t0

data <- data_generator_vd(beta_index = 1)
formula <- y ~ ffvd(X_se, nbasis = c(30, 50, 30)) + ffvd(Y_se, nbasis = c(30, 50, 30))
res <- vd_fit(formula = formula, data = data)

data_not_aligned <- data_generator_vd(aligned = FALSE, beta_index = 2)
formula <- y ~ ffvd(X_se, nbasis = c(30, 50, 30))
t0 <- proc.time()
res_not_aligned <- vd_fit(formula = formula, data = data_not_aligned)
t1 <- proc.time()

time_not_aligned <- t1 - t0


# 20241007
data = data_generator_vd(beta_index = 2, use_x = TRUE)
data = data_generator_vd(beta_index = 2, use_x = TRUE, use_f = TRUE)
formula = y ~ ffvd(X_se, nbasis = c(10, 10, 10)) + f(x1, nseg = 30, pord = 2, degree = 3) + x2 # + f(x2, nseg = 30, pord = 2, degree = 3)
family = stats::gaussian()
offset = NULL

res <- vd_fit(formula, data)




