heatmap_betas <- function(data, beta_index = 1, ...) {
  Beta <- NULL

  data <- data$Beta[,,beta_index]

  N <- nrow(data)

  M <- t(apply(data, 1, function(x) range(which(!is.na(x)))))

  max_M <- max(M)
  T_dat <- IND <- NULL
  t_dat <- rep(1:max_M, N)

  Beta_estimated <- t(data)
  dim(Beta_estimated) <- c(nrow(data) * ncol(data), 1)

  for (ind in 1:N) {
    T_dat <- c(T_dat, rep(M[ind, 2], max_M))
    IND <- c(IND, rep(ind, M[ind, 2]))
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
    ggplot2::labs(y = "T") +
    ggplot2::ggtitle(paste0("FFVD HeatMap (Beta ", beta_index, ")"))
}

heatmap_betas(data, beta_index = 1)

sum((data$Beta[,,1] - res$Beta_ffvd[[1]])^2, na.rm = TRUE)/(100*101)
sum((data$y - res$fit$fitted.values)^2, na.rm = TRUE)/(100)

plotly::plot_ly(z = data$Beta[,,1], type = "surface")
plotly::plot_ly(z = res$Beta_ffvd[[1]], type = "surface")


