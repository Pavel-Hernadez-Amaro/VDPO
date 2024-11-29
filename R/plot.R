#' @export
plot.vd_fit <- function(x, beta_index = 1, ...) {
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

  Beta <- NULL

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

#' Plot Functional Curves with Confidence Intervals
#'
#' Generates a plot of functional Beta estimates for specified curves, along with their 95% confidence intervals.
#' This function computes the 95% confidence intervals for each curve based on the covariance matrix and the fitted values
#' from the provided object. The resulting plot includes estimated curves, confidence interval ribbons, and a legend
#' distinguishing the curves.
#'
#' @param object An object of class `'vd_fit'` or similar, containing the fitted model results, Beta estimates, and evaluation details.
#' @param beta_index An integer specifying which Beta coefficient matrix to use. Default is 1.
#' @param curves A numeric vector specifying the indices of the curves (rows) to plot.
#'
#' @return A `ggplot2` object displaying the Beta estimates and confidence intervals for the specified curves.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Assuming `model_object` is an object of class 'vd_fit'
#'   plot_functional_curves_combined(model_object, beta = 1, curves = c(50, 70, 100))
#' }
#' }
#'
#' @export
plot_ci <- function(object, beta_index = 1, curves) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for this functionality", call. = FALSE)
  }

  if (beta_index < 1 || beta_index > length(x$M)) {
    stop(
      "'beta_index' should be between 1 and the number of variable domain
         functional variables used in the formula",
      call. = FALSE
    )
  }

  # Bindings for global variables
  x <- NULL
  Domain <- NULL
  lower <- NULL
  upper <- NULL
  estimated <- NULL

  N <- length(object$fit$fitted.values)
  Beta_FFVD_inf <- Beta_FFVD_sup <- matrix(nrow = nrow(object$Beta[[beta_index]]), ncol = ncol(object$Beta[[beta_index]]))
  M <- object$M[[1]][, 2]

  for (aux_ind in 1:N) {
    prod <- as.matrix(kronecker(object$ffvd_evals[[1]]$L_Phi$B[1:M[aux_ind], ], t(object$ffvd_evals[[1]]$B_T$B[aux_ind, ])))
    var_curve <- diag(prod %*% object$covar_theta %*% t(prod))
    std_curve <- sqrt(var_curve)
    Beta_FFVD_inf[aux_ind, 1:M[aux_ind]] <- object$Beta[[1]][aux_ind, 1:M[aux_ind]] - 1.96 * std_curve
    Beta_FFVD_sup[aux_ind, 1:M[aux_ind]] <- object$Beta[[1]][aux_ind, 1:M[aux_ind]] + 1.96 * std_curve
  }

  plot_data_list <- lapply(curves, function(curve) {
    x_seq <- seq_len(ncol(object$Beta[[beta_index]]))
    data.frame(
      t = x_seq,
      estimated = object$Beta[[beta_index]][curve,],
      lower = Beta_FFVD_inf[curve,],
      upper = Beta_FFVD_sup[curve,],
      Domain = paste("Domain:", curve)
    )
  })

  plot_data <- do.call(rbind, plot_data_list)

  case_colors <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd", "#8c564b")[1:length(curves)]
  names(case_colors) <- paste("Domain:", curves)

  # Create the plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = t, group = Domain)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = Domain),
                alpha = 0.2, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = estimated, color = Domain),
              size = 1, linetype = "solid", na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray50", na.rm = TRUE) +
    ggplot2::scale_color_manual(name = "Domains:",
                       values = case_colors) +
    ggplot2::scale_fill_manual(name = "Domains:",
                      values = case_colors) +
    ggplot2::labs(x = "t",
         y = "Value") +
    ggplot2::annotate("text", x = 50, y = 0.6,
             label = "- Estimated",
             size = 3, hjust = 0.5, color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top",
          panel.grid.minor = ggplot2::element_blank(),
          legend.box = "horizontal",
          plot.title = ggplot2::element_blank(),
          plot.subtitle = ggplot2::element_blank())
}
