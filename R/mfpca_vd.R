#' Multivariate functional principal component analysis for variable domain data
#'
#' Performs a multivariate functional principal component analysis (MFPCA) for
#' functional data observed on subject-specific (variable) domains. Each
#' functional variable is first decomposed through a variable domain FPCA, and
#' the resulting univariate scores are combined into a domain-varying
#' multivariate decomposition.
#'
#' The observation times can be supplied in two mutually exclusive ways. With
#' \code{Times}, each variable carries a matrix (subjects in rows, observation
#' points in columns) holding the actual observation time of every measurement,
#' and the domain of subject \eqn{i} is taken as the maximum observed time. With
#' \code{M_grid}, the observations are assumed equidistant and \code{M_grid} is a
#' vector of length \eqn{N} giving the domain length of each subject, so subject
#' \eqn{i} is evaluated on \code{seq(0, M_grid[i], by = 1 / Hz)}. Exactly one of
#' \code{Times} or \code{M_grid} must be provided.
#'
#' @param Data A list (one entry per functional variable) of matrices with
#'   subjects in rows and observation points in columns, with \code{NA} on the
#'   unobserved part of each domain. All variables must have the same number of
#'   subjects. Currently two variables are supported.
#' @param Times A list of matrices, the same shape as \code{Data}, giving the
#'   actual observation time of each measurement. Mutually exclusive with
#'   \code{M_grid}.
#' @param M_grid A numeric vector of length \eqn{N} with the domain length of
#'   each subject, used when the observations are equidistant. Mutually exclusive
#'   with \code{Times}.
#' @param Hz Sampling rate used to build the equidistant grid when \code{M_grid}
#'   is supplied (default \code{1}).
#' @param m_npcs Number of multivariate principal components to retain. If
#'   \code{NULL}, the number explaining more than 90\% of the variance is used.
#' @param u_npcs Number of univariate principal components retained per variable
#'   (default \code{5}).
#' @param k_m Dimension of the basis used to model the score covariances along
#'   the domain (default \code{15}).
#' @param model_type Either \code{"gam"} (the default) or \code{"sop"}, the
#'   method used to smooth the score covariances along the domain.
#'
#' @return A list with the following items:
#'
#' - \code{scores_m}: matrix of multivariate scores.
#' - \code{efunctions_m}: multivariate eigenfunctions, as a list over domains,
#'   each a list over variables.
#' - \code{efunctions_u}: univariate eigenfunctions, as a list over variables.
#' - \code{scores_u}: matrix of univariate scores.
#' - \code{evalues_u}: univariate eigenvalues, as a list over variables.
#' - \code{evalues_m}: multivariate eigenvalues, as a list over domains.
#' - \code{var_u}: cumulative variance explained by the univariate components.
#' - \code{mean_model}: the fitted univariate mean models (one per variable).
#' - \code{M_grid}: the domain grid used.
#'
#' @importFrom mgcv gam
#' @importFrom stats predict
#' @export
mfpca_vd <- function(Data, Times = NULL, M_grid = NULL, Hz = 1,
                     m_npcs = NULL, u_npcs = 5, k_m = 15, model_type = "gam") {

  model_type <- match.arg(model_type, c("gam", "sop"))

  if (is.null(Times) && is.null(M_grid)) {
    stop("Either 'Times' or 'M_grid' must be provided.", call. = FALSE)
  }
  if (!is.null(Times) && !is.null(M_grid)) {
    stop("Provide only one of 'Times' or 'M_grid', not both.", call. = FALSE)
  }

  nobs_subjects <- sapply(Data, nrow)
  stopifnot("Same sample size for all variables required" = all(nobs_subjects == nobs_subjects[1]))

  n_var <- length(Data)
  N <- nrow(Data[[1]])

  fit_Variable <- e_len_u <- npc_u <- var_u <- univ_efunctions <-
    vector(mode = "list", length = n_var)
  scores <- NULL

  for (ind in seq_len(n_var)) {

    ## ---- long format (id, time, y) ----
    long <- .vd_long(Data[[ind]], if (is.null(Times)) NULL else Times[[ind]], M_grid, Hz)

    ## ---- domain per subject, reindexed by increasing domain ----
    maxT_by_id <- tapply(long$time, long$id, max)
    dom <- data.frame(id = as.integer(names(maxT_by_id)), maxT = as.numeric(maxT_by_id))
    dom <- dom[dom$maxT > 0, , drop = FALSE]
    dom <- dom[order(dom$maxT), , drop = FALSE]
    dom$newid <- seq_len(nrow(dom))

    long <- merge(long, dom, by = "id")
    long$id <- long$newid
    long <- long[order(long$id, long$time), c("id", "maxT", "time", "y")]

    ## ---- univariate mean and residuals ----
    fit_Variable[[ind]] <- mgcv::gam(y ~ s(time, maxT), data = long, method = "REML")
    long$mean  <- as.vector(predict(fit_Variable[[ind]]))
    long$ydiff <- long$y - long$mean

    ## ---- covariance setup and fit ----
    cov_dat  <- .datsetup_cov(long)
    cov_temp <- mgcv::gam(kprod ~ s(stime, ttime, maxT), data = cov_dat)

    ## ---- domain grid for the variable domain FPCA ----
    if (!is.null(M_grid)) {
      gridM    <- M_grid
      argvals  <- NULL
    } else {
      gridM    <- apply(Times[[ind]], 1, max, na.rm = TRUE)
      argvals  <- Times[[ind]]
    }

    pc_grid <- .get_pcs_vd(xtimes = gridM, covfx = cov_temp, argvals = argvals,
                           Hz = Hz, includezero = TRUE, npcs = u_npcs)

    univ_efunctions[[ind]] <- pc_grid$efunctions
    npc_aux <- sapply(pc_grid$evalues, length)
    stopifnot("Same number of principal components required for all the sample curves" =
                all(npc_aux == npc_aux[1]))
    npc_u[[ind]]   <- npc_aux
    e_len_u[[ind]] <- pc_grid$evalues
    var_u[[ind]]   <- pc_grid$cvar

    scores_list <- .get_scores_vd(Y = Data[[ind]], xtimes = gridM,
                                  univ_efunctions = univ_efunctions[[ind]],
                                  npc = npc_u[[ind]], argvals = argvals, Hz = Hz,
                                  mean_func = fit_Variable[[ind]], includezero = TRUE)

    score_matrix <- matrix(NA_real_, nrow = N, ncol = pc_grid$npc[[1]])
    for (i in seq_len(N)) score_matrix[i, ] <- scores_list[[i]]

    scores <- cbind(scores, score_matrix)
  }

  ## ---- per-subject outer products of stacked scores ----
  Cov_e <- lapply(seq_len(nrow(scores)), function(x) outer(scores[x, ], scores[x, ]))
  for (s in seq_along(Cov_e)) Cov_e[[s]][Cov_e[[s]] < 0] <- 0

  if (!is.null(M_grid)) {
    Domain <- M_grid
  } else {
    Domain <- unname(apply(cbind(Times[[1]], Times[[2]]), 1, max, na.rm = TRUE))
  }

  p <- nrow(Cov_e[[1]])

  ## ---- one score-covariance model per (i, j), j <= i ----
  comp_models <- list()
  pos <- 0
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      yij <- vapply(seq_len(N), function(x) Cov_e[[x]][i, j], numeric(1))
      df  <- data.frame(y = yij, M = Domain)
      m   <- if (model_type == "gam") {
        mgcv::gam(y ~ s(M, bs = "ps", k = k_m), data = df)
      } else {
        fml <- stats::as.formula(sprintf("y ~ f(M, nseg = %s)", k_m))
        environment(fml) <- asNamespace("SOP")
        SOP::sop(formula = fml, data = df)
      }
      pos <- pos + 1
      comp_models[[pos]] <- list(i = i, j = j, model = m)
    }
  }

  ## ---- domain-varying multivariate covariance ----
  Cov_T <- vector(mode = "list", length = length(Domain))
  for (d in seq_along(Domain)) {
    Ct <- matrix(0, nrow = p, ncol = p)
    for (cm in comp_models) {
      pr <- as.numeric(predict(cm$model, data.frame(M = Domain[d])))
      if (pr < 0) pr <- 0
      Ct[cm$i, cm$j] <- pr
    }
    Ct <- Ct + t(Ct) - diag(diag(Ct))
    Cov_T[[d]] <- Ct
  }

  ## ---- multivariate eigendecomposition per domain ----
  v_m <- c_m <- vector(mode = "list", length = length(Cov_T))
  for (i in seq_along(Cov_T)) {
    e <- eigen(Cov_T[[i]])
    sub_npcs <- if (is.null(m_npcs)) min(which(cumsum(e$values) / sum(e$values) > 0.9)) else m_npcs
    sub_npcs <- min(sub_npcs, Domain[i])
    c_m[[i]] <- e$vectors[, seq_len(sub_npcs), drop = FALSE]
    v_m[[i]] <- e$values[seq_len(sub_npcs)]
  }

  ## ---- multivariate eigenfunctions and scores ----
  npc_per_var <- vapply(npc_u, function(v) v[1], integer(1))
  offsets     <- cumsum(c(0L, npc_per_var))

  Phi_m   <- vector(mode = "list", length = length(Domain))
  score_m <- matrix(nrow = N, ncol = if (is.null(m_npcs)) ncol(c_m[[1]]) else m_npcs)

  for (i in seq_along(Domain)) {
    Phi_m[[i]] <- vector(mode = "list", length = n_var)
    for (j in seq_len(n_var)) {
      block <- (offsets[j] + 1):offsets[j + 1]
      Phi_m[[i]][[j]] <- univ_efunctions[[j]][[i]] %*% as.matrix(c_m[[i]][block, ])
    }
    score_m[i, ] <- t(as.matrix(scores[i, ])) %*% c_m[[i]]
  }

  list(
    scores_m     = score_m,
    efunctions_m = Phi_m,
    efunctions_u = univ_efunctions,
    scores_u     = scores,
    evalues_u    = e_len_u,
    evalues_m    = v_m,
    var_u        = var_u,
    mean_model   = fit_Variable,
    M_grid       = gridM
  )
}

#' Build the long format (id, time, y) for one variable
#' @noRd
.vd_long <- function(Dmat, Tmat, M_grid, Hz) {
  N <- nrow(Dmat)
  out <- vector("list", N)
  for (i in seq_len(N)) {
    obs <- which(!is.na(Dmat[i, ]))
    if (!is.null(Tmat)) {
      obs <- obs[!is.na(Tmat[i, obs])]
      tt  <- Tmat[i, obs]
    } else if (!is.null(M_grid)) {
      tt  <- seq(0, by = 1 / Hz, length.out = length(obs))
    } else {
      tt  <- obs
    }
    out[[i]] <- data.frame(id = i, time = as.numeric(tt), y = as.numeric(Dmat[i, obs]))
  }
  long <- do.call(rbind, out)
  long <- long[!is.na(long$y), , drop = FALSE]
  long[order(long$id, long$time), , drop = FALSE]
}

#' Pairwise residual products for the covariance fit
#' @noRd
.datsetup_cov <- function(d) {
  ids <- unique(d$id)
  out <- vector("list", length(ids))
  for (k in seq_along(ids)) {
    di <- d[d$id == ids[k], , drop = FALSE]
    m  <- nrow(di)
    out[[k]] <- data.frame(
      id    = ids[k],
      maxT  = di$maxT[1],
      stime = rep(di$time, m),
      ttime = rep(di$time, each = m),
      kprod = as.vector(kronecker(di$ydiff, di$ydiff))
    )
  }
  do.call(rbind, out)
}

#' Trapezoidal quadrature weights
#' @noRd
.quad_weights <- function(argvals) {
  D <- length(argvals)
  if (D == 2) {
    0.5 * c(argvals[2] - argvals[1], argvals[D] - argvals[D - 1])
  } else {
    0.5 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 2)],
            argvals[D] - argvals[D - 1])
  }
}

#' Prediction data frame for a covariance slice
#' @noRd
.create_new_data <- function(maxT, Hz = NULL, argvals = NULL, includezero = TRUE) {
  stopifnot(!is.null(Hz) || !is.null(argvals))
  if (is.null(argvals)) {
    time_vals <- if (includezero) seq(0, maxT, by = 1 / Hz) else seq(0, maxT, by = 1 / Hz)[-1]
  } else {
    time_vals <- argvals
  }
  npts <- length(time_vals)
  data.frame(stime = rep(time_vals, npts),
             ttime = rep(time_vals, each = npts),
             maxT  = rep(maxT, npts * npts))
}

#' Covariance slice at a given domain length
#' @noRd
.pcov_m <- function(covx, times, argvals = NULL, Hz = NULL, includezero = TRUE) {
  stopifnot(!is.null(Hz) || !is.null(argvals))
  out <- vector("list", length(times))
  counter <- 1
  for (i in times) {
    temp      <- .create_new_data(i, Hz = Hz, argvals = argvals, includezero = includezero)
    temp_pred <- predict(covx, temp)
    npts      <- if (is.null(argvals)) {
      length(if (includezero) seq(0, i, by = 1 / Hz) else seq(0, i, by = 1 / Hz)[-1])
    } else {
      length(argvals)
    }
    pred_mat <- matrix(temp_pred, ncol = npts)
    out[[counter]] <- (pred_mat + t(pred_mat)) / 2
    counter <- counter + 1
  }
  out
}

#' Variable domain eigenfunctions and eigenvalues at each domain length
#' @noRd
.get_pcs_vd <- function(xtimes, covfx, argvals = NULL, Hz = NULL, npcs = NULL,
                        pve = 0.99, includezero = TRUE) {
  nsub <- length(xtimes)
  efx_list <- evals_list <- cumvar_list <- vector("list", nsub)
  npc_v <- integer(nsub)
  for (i in seq_len(nsub)) {
    maxt <- xtimes[i]
    if (!is.null(argvals)) {
      aux <- argvals[i, ]
      argvals_i <- aux[seq_len(sum(!is.na(aux)))]
      if (!includezero) argvals_i <- argvals_i[argvals_i > 0]
    } else {
      argvals_i <- if (includezero) seq(0, maxt, by = 1 / Hz) else seq(0, maxt, by = 1 / Hz)[-1]
    }
    npc.0    <- .pcov_m(covx = covfx, times = maxt, argvals = argvals_i,
                        Hz = Hz, includezero = includezero)[[1]]
    qw       <- .quad_weights(argvals_i)
    Wsqrt    <- diag(sqrt(qw))
    Winvsqrt <- diag(1 / sqrt(qw))
    V        <- Wsqrt %*% npc.0 %*% Wsqrt
    evalues  <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues2 <- replace(evalues, which(evalues <= 0), 0)
    var_pct  <- evalues2 / sum(evalues2)
    cumvar   <- cumsum(var_pct)
    sub_npcs <- if (is.null(npcs)) min(which(cumsum(evalues2) / sum(evalues2) > pve)) else npcs
    sub_npcs <- min(sub_npcs, length(argvals_i))
    efunctions <- matrix(
      Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq_len(sub_npcs)],
      nrow = length(argvals_i), ncol = sub_npcs
    )
    npc_v[i]         <- sub_npcs
    efx_list[[i]]    <- efunctions
    evals_list[[i]]  <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[seq_len(sub_npcs)]
    cumvar_list[[i]] <- cumvar
  }
  list(evalues = evals_list, efunctions = efx_list, npc = npc_v, cvar = cumvar_list)
}

#' Variable domain scores by numerical integration
#' @noRd
.get_scores_vd <- function(Y, xtimes, univ_efunctions, npc, Hz = NULL, argvals = NULL,
                           mean_func = NULL, includezero = TRUE) {
  stopifnot(!is.null(Hz) || !is.null(argvals))
  stopifnot(!is.null(mean_func))
  nsub <- length(xtimes)
  scores_list <- vector("list", nsub)
  for (i in seq_len(nsub)) {
    maxt <- xtimes[i]
    if (!is.null(argvals)) {
      aux <- argvals[i, ]
      subj_argvals <- aux[seq_len(sum(!is.na(aux)))]
      if (!includezero) subj_argvals <- subj_argvals[subj_argvals > 0]
    } else {
      subj_argvals <- if (includezero) seq(0, maxt, by = 1 / Hz) else seq(0, maxt, by = 1 / Hz)[-1]
    }
    efunctions <- univ_efunctions[[i]]
    K  <- npc[i]
    qw <- .quad_weights(subj_argvals)
    mu_i       <- predict(mean_func, newdata = data.frame(time = subj_argvals, maxT = maxt))
    centered_Y <- Y[i, seq_along(subj_argvals)] - mu_i
    subj_scores <- numeric(K)
    for (k in seq_len(K)) subj_scores[k] <- sum(centered_Y * efunctions[, k] * qw)
    scores_list[[i]] <- subj_scores
  }
  scores_list
}
