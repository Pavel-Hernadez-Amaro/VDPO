px <- 20
py <- 20
N <- 100
Rsq <- 0.95
x <- seq(from = 0, to = 1, length.out = px)
y <- seq(from = 0, to = 1, length.out = py)
# e_1 = 0.2;e_2 = 0.2;e_data = 0.02
#
# a_data <- matrix(0, nrow = length(x), ncol = length(y))
# for (i in seq_along(x)) {
#   for (j in seq_along(y)) {
#     a_data[i, j] = rnorm(1, 0, e_data)
#   }
# }
#
# a1 <- rnorm(1,0,e_1)
# a2 <- rnorm(1,0,e_2)


X <- lapply(1:N, function(a)
  list(
    DATA_T = Data_H(x, y)$DATA_T,
    DATA_N = Data_H(x, y)$DATA_N
  )
)


X_true <- lapply(seq_along(X), function (i) X[[i]]$DATA_T)
X_real <- lapply(seq_along(X), function (i) X[[i]]$DATA_N)



all.equal(X[[1]]$DATA_T, X[[2]]$DATA_T)
all.equal(X[[1]]$DATA_N, X[[2]]$DATA_N)

nu <- sapply(1:5, function(a) response_int_H(Data_H, f_Beta = Beta_H_saddle, x = x, y = y))


var_e <- (1/Rsq - 1) * var(nu) #(1-Rsq)*var(nu[ind,])

response=nu+rnorm(N,sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL #rnorm(nu[ind,])
