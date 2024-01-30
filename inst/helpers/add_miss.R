#' Set missing values to the surfaces
#'
#' @param X List of surfaces.
#' @param n_missing Number of holes in every surface.
#' @param min_distance_x Length of the holes in the x-axis.
#' @param min_distance_y Length of the holes in the y-axis.
#'
#' @return List of surfaces with missing values.
#'
#' @noRd
add_miss <- function(X, n_missing = 1, min_distance_x = 9, min_distance_y = 9) {
  N <- length(X)
  missing_points=miss_points=vector(mode = "list", length = N)

  x_b <- nrow(X[[1]])
  y_b <- ncol(X[[1]])

  for (i in 2:N) {
    if (nrow(X[[i]]) != x_b || ncol(X[[i]]) != y_b) {
      stop("The dimension of all the surfaces must be the same.", call. = FALSE)
    }
  }

  for (j in seq_along(missing_points)) {
    x_missing=sort(sample(1:(x_b-min_distance_x),n_missing))
    y_missing=sort(sample(1:(y_b-min_distance_y),n_missing))

    x_pos=x_missing:(x_missing+min_distance_x-1)
    y_pos=y_missing:(y_missing+min_distance_y-1)

    missing_points[[j]]=expand.grid(x_pos,y_pos)

    for (add_miss in 1:nrow(missing_points[[j]])) {
      X[[j]][missing_points[[j]][add_miss, 1], missing_points[[j]][add_miss, 2]] <-
        NA
    }

    miss_points[[j]]=vector(mode="list", length = ncol(X[[j]]))

    for (i in 1:ncol(X[[j]])) {

      miss_spots=NULL

      for (j_row in 1:nrow(X[[j]])) {

        if (is.na(X[[j]][j_row,i])) {

          miss_spots=c(miss_spots,j_row)

        }

      }

      if (!is.null(miss_spots)) {

        miss_points[[j]][[i]]=miss_spots

      }

    }
  }

  list(
    X_miss = X,
    miss_point = miss_points
  )

}


for (ind_miss_points in 1:N) {

  if (length(miss_points[[ind_miss_points]])!=ncol(X[[j]])) {
    stop(paste("Revisa que el length de miss point en ", ind_miss_points, "sea", ncol(X[[j]])), call. = FALSE)
  }

}

var_e <-  (1/Rsq - 1) * var(n) #(1-Rsq)*var(nu[ind,])
