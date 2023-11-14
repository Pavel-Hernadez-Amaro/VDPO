#' @export
plot.VDFO <- function(x, plot_index = 1,...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for this functionality", call. = FALSE)
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("package 'RColorBrewer' is required for this functionality", call. = FALSE)
  }

  if (plot_index < 1 || plot_index > length(x$M_ffvd)) {
    stop("'plot_index' should be between 1 and the number of variable domain
         functional variables used in the formula", call. = FALSE)
  }

  N <- attr(x, "N")

  max_M <- max(x$M_ffvd[[plot_index]])
  T_dat <- IND <- NULL
  t_dat <- rep(1:max_M, N)

  Beta_estimated <- t(x$Beta_ffvd[[plot_index]])
  dim(Beta_estimated) <- c(nrow(x$Beta_ffvd[[plot_index]]) * ncol(x$Beta_ffvd[[plot_index]]), 1)

  for (ind in 1:N) {
    T_dat <- c(T_dat, rep(x$M_ffvd[[plot_index]][ind], max_M))
    IND   <- c(IND, rep(ind, x$M_ffvd[[plot_index]][ind]))
  }

  Heat_map_data <- data.frame(t = t_dat, M = T_dat, Beta = Beta_estimated)
  Heat_map_data <- Heat_map_data[t_dat <= T_dat, ]
  Heat_map_data[["IND"]] <- IND

  lims <- range(Heat_map_data$Beta)

  ggplot2::ggplot(Heat_map_data, ggplot2::aes(x = t, y = IND)) +
    ggplot2::geom_tile(ggplot2::aes(colour=Beta, fill=Beta)) +
    ggplot2::scale_fill_gradientn(name="", limits=lims,
                                  colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    ggplot2::scale_colour_gradientn(name="", limits=lims,
                                    colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::theme_bw()+
    ggplot2::labs(y="T") +
    ggplot2::ggtitle("Beta FFVD") -> plot

  plot
}
