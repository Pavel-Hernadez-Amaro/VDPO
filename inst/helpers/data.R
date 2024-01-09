library(refund)
library(fda)
library(mgcv)
library(SOP)

#' data generator
#'
#' @param N Number of subjects.
#' @param J Number of maximum observations per subject.
#' @param R Number of simulations per the simulation study.
#' @param case Total number of scenarios simulated.
#' @param aligned If the data is aligned or not.
#' @param multivariate If TRUE, the data is generated with 2 variables.
#'
#' @return
dg <- function(N = 100, J = 100, R = 1, case = 1, Rsq = 0.95, aligned = TRUE, multivariate = FALSE) {
  for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE
    for (iter in 1:R) {
      set.seed(42+iter)

      M <-  round(runif(N, 10, J), digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS

      if (max(M) > J)
        M[which(M>J)] <- J

      if (min(M) <= 10)
        M[which(M <= 10)] <- 10

      T <- max(M)

      t <- 1:T

      M <- sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY

      ############ HERE WE GENERATE THE FUNCTIONAL DATA

      X_s  <- matrix(NA, N, T) # NOT NOISY
      X_se <- matrix(NA, N, T) # NOISY
      Y_s  <- matrix(NA, N, T) # NOT NOISY
      Y_se <- matrix(NA, N, T) # NOISY


      for (i in 1:N) {
        u1 <- rnorm(1)

        temp <- matrix(NA, 10, T)

        for (k in 1:10) {

          v_i1 <- rnorm(1, 0, 4 / k^2)
          v_i2 <- rnorm(1, 0, 4 / k^2)

          if (aligned) {
            temp[k, 1:M[i]] <-
              v_i1 * sin(2 * pi * k * (1:M[i]) / 100) + v_i2 * cos(2 * pi * k * (1:M[i]) / 100)
          } else {
            temp[k, (M[i, 1]:M[i, 2])] <-
              v_i1 * sin(2 * pi * k * (M[i, 1]:M[i, 2]) / 100) + v_i2 * cos(2 * pi * k * (M[i, 1]:M[i, 2]) / 100)
          }
        }

        B <- apply(temp,2,sum)


        B  <- B + u1
        B2 <- B + rnorm(1, sd = 0.02) + (t / 10)

        X_s[i,]  <- B
        X_se[i,] <- B + rnorm(T, 0, 1) # WE ADD NOISE
        Y_s[i,]  <- B2
        Y_se[i,] <- B2 + rnorm(T, 0, 1) # WE ADD NOISE

      }

      # SOME SAVE CHECKS FOR UNWANTED NAs
      for (i in 1:T)
        if (length(which(is.na(X_s[,i]))) == N)
          stop("some NA columns exist", call. = FALSE)

      Beta <- array(dim = c(N, T, 4))
      nu <- y <- rep(0, N)
      x1 <- runif(N)
      # f0 <- \(x) cos((x-0.5)^2)
      f0 <- \(x) 2*sin(pi*x)

      for (i in 1:N) {
        # TRUE FUNCTIONAL COEFFICIENTS

        Beta[i, 1:(M[i]), 1] <- ((10 * t[1:(M[i])] / M[i]) - 5) / 10
        Beta[i, 1:(M[i]), 2] <- ((1 - (2 * M[i] / T)) * (5 - 40 * ((t[1:(M[i])] / M[i]) - 0.5) ^ 2)) / 10

        # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)

        if (multivariate) {
          if (iter_out == 8) {
            nu[i] <- sum(X_se[i, ] * Beta[i, , 4], na.rm = 1) / (M[i]) + sum(Y_se[i, ] * Beta[i, , 4], na.rm = 1) / (M[i]) # NOISY
          } else if (iter_out <= 4) {
            nu[i] <-  sum(X_s[i, ] * Beta[i, , 1], na.rm = 1) / (M[i]) + sum(Y_s[i, ] * Beta[i, , 2], na.rm = 1) / (M[i]) #NOT NOISY
          } else if (iter_out > 4 && iter_out < 8) {
            nu[i] <- sum(X_se[i, ] * Beta[i, , iter_out %% 4], na.rm = 1) / M[i] + sum(Y_se[i, ] * Beta[i, , iter_out %% 4], na.rm = 1) / M[i] # NOISY
          }
        } else {
          if (iter_out == 8) {
            nu[i] <- sum(X_se[i, ] * Beta[i, , 4], na.rm = 1) / (M[i]) # NOISY
          } else if (iter_out <= 4) {
            nu[i] <-  sum(X_s[i, ] * Beta[i, , 1], na.rm = 1) / (M[i]) #NOT NOISY
          } else if (iter_out > 4 && iter_out < 8) {
            nu[i] <- sum(X_se[i, ] * Beta[i, , iter_out %% 4], na.rm = 1) / M[i] # NOISY
          }
        }

      }

      nu <- nu + f0(x1)
      var_e <- (1/Rsq - 1) * var(nu)
      # y <- x1 + nu + rnorm(N, sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL
      y <- nu + rnorm(N, sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL
    }
  }

  data = data.frame(y = y)
  data[["X_s"]]  <- X_s
  data[["X_se"]] <- X_se
  data[["Y_s"]]  <- Y_s
  data[["Y_se"]] <- Y_se
  data[["x1"]] <- x1
  data[["Beta"]] <- Beta

  # añadir rnorm. min(range(X)) sea el mean del rnorm y el std el std de las X o del range
  data
}


h1 <- function(t) {
  res <- c()
  for (value in t) {
    res <- c(res, max(25 - abs(value - 50), 0))
  }
  res
}
h2 <- function(t) h1(t - 20)
h3 <- function(t) h1(t + 20)

x <- function(t, cl = 1) {
  h <- if (cl == 1) h2 else h3
  u <- runif(1)
  e <- rnorm(length(t), mean = 0, sd = 0.1)

  u * h1(t) + (1 - u) * h(t) + e
}
