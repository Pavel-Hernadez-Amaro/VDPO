#' @export
plot.vd_fit <- function(x, beta_index = 1, ...) {
  Beta <- NULL

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for this functionality", call. = FALSE)
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("package 'RColorBrewer' is required for this functionality", call. = FALSE)
  }

  if (beta_index < 1 || beta_index > length(x$M)) {
    stop(
      "'beta_index' should be between 1 and the number of variable domain
         functional variables used in the formula",
      call. = FALSE
    )
  }

  N <- attr(x, "N")

  max_M <- max(x$M[[beta_index]])
  T_dat <- IND <- NULL
  t_dat <- rep(1:max_M, N)

  Beta_estimated <- t(x$Beta[[beta_index]])
  dim(Beta_estimated) <- c(nrow(x$Beta[[beta_index]]) * ncol(x$Beta[[beta_index]]), 1)

  for (ind in 1:N) {
    T_dat <- c(T_dat, rep(x$M[[beta_index]][ind, 2], max_M))
    IND <- c(IND, rep(ind, x$M[[beta_index]][ind, 2]))
  }

  Heat_map_data <- data.frame(t = t_dat, M = T_dat, Beta = Beta_estimated)
  Heat_map_data <- Heat_map_data[t_dat <= T_dat, ]
  Heat_map_data[["IND"]] <- IND

  lims <- range(Heat_map_data$Beta)

  ggplot2::ggplot(Heat_map_data, ggplot2::aes(x = t, y = IND)) +
    ggplot2::geom_tile(ggplot2::aes(colour = Beta, fill = Beta)) +
    ggplot2::scale_fill_gradientn(
      name     = "",
      limits   = lims,
      colours  = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      na.value = "white"
    ) +
    ggplot2::scale_colour_gradientn(
      name    = "",
      limits  = lims,
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "T")
}


plot_beta_with_ci <- function(res, curve = 1) {
  if (curve > nrow(res$Beta[[1]])) {
    stop("Requested curve index is out of range.")
  }

  N <- length(res$fit$fitted.values)
  Beta_FFVD_inf <- Beta_FFVD_sup <- matrix(nrow = nrow(res$Beta[[1]]), ncol = ncol(res$Beta[[1]]))
  M <- res$M[[1]][,2]

  for (aux_ind in 1:N) {
    prod = as.matrix(kronecker(res$ffvd_evals[[1]]$L_Phi$B[1:M[aux_ind],], t(res$ffvd_evals[[1]]$B_T$B[aux_ind,])))
    var_curve = diag(prod %*% res$covar_theta %*% t(prod))
    std_curve = sqrt(var_curve)
    Beta_FFVD_inf[aux_ind,1:M[aux_ind]] = res$Beta[[1]][aux_ind,1:M[aux_ind]] - 1.96 * std_curve
    Beta_FFVD_sup[aux_ind,1:M[aux_ind]] = res$Beta[[1]][aux_ind,1:M[aux_ind]] + 1.96 * std_curve
  }

  plot_data <- data.frame(
    index = rep(1:length(res$Beta[[1]][curve,]), 3),
    value = c(res$Beta[[1]][curve,], Beta_FFVD_inf[curve,], Beta_FFVD_sup[curve,]),
    group = rep(c("Estimate", "Lower CI", "Upper CI"), each = length(res$Beta[[1]][curve,]))
  )

  # Filter out missing values and values outside the scale range
  plot_data <- na.omit(plot_data)
  plot_data <- subset(plot_data, value >= min(plot_data$value) & value <= max(plot_data$value))

  ggplot(plot_data, aes(x = index, y = value, color = group, linetype = group)) +
    geom_line() +
    scale_color_manual(values = c("black", "#0072B2", "#009E73")) +
    scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "res$Beta[[1]][100,]",
      title = paste0("Confidence Intervals for the ", curve, "th Row of res$Beta[[1]]")
    )
}

plot <- plot_beta_with_ci(res, 10)
plot
