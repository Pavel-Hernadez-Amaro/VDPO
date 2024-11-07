data <- generate_2d_po_functional_data(noise_sd = 0.1)

data$noisy_surfaces_miss[[1]] -> data_miss
data$noisy_surfaces_miss[[2]] -> miss_points
data$noisy_surfaces_miss[[3]] -> missing_points

res_2 = ffpo_2d(data_miss, miss_points = miss_points, missing_points = missing_points)

case=100

plotly::plot_ly(z=data$surfaces[[case]], type="surface")
plotly::plot_ly(z=data$noisy_surfaces_miss$X_miss[[case]], type="surface")
plotly::plot_ly(z=res_2$X_hat[[case]], type="surface")

data_1d <- generate_1d_po_functional_data()
data_1d |> names()

res = ffpo(X = data_1d$curves, grid = data_1d$grid)

case=1

plot(data_1d$curves[case,],type="l")
plot(data_1d$noisy_curves_miss[case,],type="l")
plot(res$[case,],type="l")




