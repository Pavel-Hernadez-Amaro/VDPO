response_integrator <- function(..., fx, fb, sub = 50, index = 1) {
  if (!requireNamespace("matrixcalc"))
    stop("package 'matrixcalc' is required for this functionality", call. = FALSE)

  dots <- list(...)
  N <- dots[["N"]]
  M <- dots[["M"]]

  n_y <- 2 * sub
  x <- seq(min(M), max(M), by = n_y)

  dg <- fx(N, M, precision = n_y)
  bg <- fb(N, M, precision = n_y)

  h <- seq(min(M), max(M), by = ncol(dg$X_se))  # if y_b and y_a are functions of x_observations

  simp_w_y <- rep(1,ncol(dg$X_se))
  even_y <- seq(2,ncol(dg$X_se)-1,2)
  odd_y <- seq(3,ncol(dg$X_se)-1,2)

  simp_w_y[even_y] <- 2
  simp_w_y[odd_y] <- 4

  Sim_w_x_y <- (h / 3) * simp_w_y

  dg$X_se <- ifelse(is.na(dg$X_se),
                    0, dg$X_se)
  bg[, , index] <- ifelse(is.na(bg[, , index]),
                          0, bg[, , index])


  sum(dg$X_se %*% diag(Sim_w_x_y) %*% t(bg[, , index]))

  # as.double(t(matrixcalc::vec(dg$X_se)) %*% diag(Sim_w_x_y) %*% matrixcalc::vec(bg[, , index]))

}


data_generator <- function(N, M, precision = 1) {
  precision <- 1 / precision
  sj <- seq(min(M), max(M), by = precision)
  X_se <- matrix(NA,N,length(sj)) # NOISY
  X_s <- matrix(NA,N,length(sj)) # NOT NOISY

  for (i in 1:N) {
    s <- seq(M[i,1], M[i,2], by = precision)
    u <- rnorm(1)

    temp <- matrix(NA, 10, length(s))


    for (k in 1:10) {

      v_i1 <- rnorm(1,0,4/k^2)
      v_i2 <- rnorm(1,0,4/k^2)

      temp[k, ] <-
        v_i1*sin(2*pi*k* s / 100)+v_i2*cos(2*pi*k*s/100)
    }

    B <- colSums(temp)[which(!is.na(colSums(temp)))]
    B <- B+u

    X_s[i,s] <- B
    X_se[i,s] <- (B)+rnorm(length(B), 0, 1) # WE ADD NOISE

  }

  list(
    X_s  = X_s,
    X_se = X_se
  )
}

beta_generator <- function(N, M, precision = 1) {
  precision <- 1 / precision
  rangeM <- seq(min(M), max(M), by = precision)
  maxM <- max(M)
  Beta <- array(dim = c(N,length(rangeM),4))
  nu <- y <- rep(0,N)

  for (i in 1:N) {
    s <- seq(M[i,1], M[i,2], by = precision)

    # TRUE FUNCTIONAL COEFFICIENTS

    Beta[i, s, 1] <- ((10*rangeM[s]/precision)-5)/10
    Beta[i, s, 2] <- ((1-(2*precision/maxM))*(5-40*((rangeM[s]/precision)-0.5)^2))/10
    Beta[i, s, 3] <- (5-10*((precision-rangeM[s])/maxM))/10
    Beta[i, s, 4] <- (sin(2*pi*precision/maxM)*(5-10*((precision-rangeM[s])/maxM)))/10

  }

  Beta
}


#' data generator
#'
#' @param N Number of subjects.
#' @param J Number of maximum observations per subject.
#' @param R Number of simulations per the simulation study.
#' @param case Total number of scenarios simulated.
#' @param aligned If the data is aligned or not.
#'
#' @return
dg <- function(N = 10, J = 100, R = 1, case = 1, aligned = TRUE) {
  for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE
    print(c("case = ",iter_out))

    for (iter in 1:R) {
      print(c("iter = ",iter))
      set.seed(42+iter)

      M <-  round(runif(N, 10, J), digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS


      if (max(M) > J)
        M[which(M>J)] <- J


      if (min(M) <= 10) {
        M[which(M <= 10)] <- 10
      }

      T <- max(M)

      t <- 1:T

      M <- sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY

      ############ HERE WE GENERATE THE FUNCTIONAL DATA

      X_se <- matrix(NA, N, T) # NOISY

      X_s <- matrix(NA, N, T) # NOT NOISY

      for (i in 1:N) {
        u <- rnorm(1)

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

        B <- B+u

        X_s[i,] <- B
        X_se[i,] <- B+rnorm(T, 0, 1) # WE ADD NOISE

      }

      # SOME SAVE CHECKS FOR UNWANTED NAs
      for (i in 1:T)
        if (length(which(is.na(X_s[,i]))) == N)
          stop("some NA columns exist", call. = FALSE)



      Beta <- array(dim = c(N, T, 4))
      nu <- y <- rep(0, N)

      for (i in 1:N) {

        # TRUE FUNCTIONAL COEFFICIENTS

        Beta[i, 1:(M[i]), 1] <- ((10 * t[1:(M[i])] / M[i]) - 5) / 10
        Beta[i, 1:(M[i]), 2] <- ((1 - (2 * M[i] / T)) * (5 - 40 * ((t[1:(M[i])] / M[i]) - 0.5) ^ 2)) / 10
        Beta[i, 1:(M[i]), 3] <- (5 - 10 * ((M[i] - t[1:(M[i])]) / T)) / 10
        Beta[i, 1:(M[i]), 4] <- (sin(2 * pi * M[i] / T) * (5 - 10 * ((M[i] - t[1:(M[i])]) / T))) / 10
        #


        # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)

        if (iter_out == 8) {
          nu[i] <- sum(X_se[i, ] * Beta[i, , 4], na.rm = 1) / (M[i]) # NOISY
        } else if (iter_out <= 4) {
          nu[i] <-  sum(X_s[i, ] * Beta[i, , iter_out], na.rm = 1) / (M[i]) #NOT NOISY
        } else if(iter_out > 4 && iter_out < 8){
          nu[i] <- sum(X_se[i, ] * Beta[i, , iter_out %% 4], na.rm = 1) / M[i] # NOISY
        }

      }


      y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE GAUSSIAN MODEL
    }
  }

  data = data.frame(y = y)
  data[["X"]] <- X_se
  data[["Y"]] <- X_s

  data
}
