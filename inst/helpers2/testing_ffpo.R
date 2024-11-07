data <- generate_2d_po_functional_data()

data$noisy_surfaces_miss[[1]] -> data_miss
data$noisy_surfaces_miss[[2]] -> miss_points
data$noisy_surfaces_miss[[3]] -> missing_points

res_2 = ffpo_2d(data_miss, miss_points = miss_points, missing_points = missing_points)

data_1d <- generate_1d_po_functional_data()
data_1d |> names()

res = ffpo(X = data_1d$curves, grid = data_1d$grid)
