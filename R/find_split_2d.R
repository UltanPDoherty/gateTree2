find_split_2d <- function(
    x, y, slope_type, mechanism, degrees = NULL, boundary_noise_comp = TRUE) {
  if (is.null(degrees)) {
    if (slope_type == "incline") {
      degrees <- seq(15, 75, by = 10)
    } else if (slope_type == "vertical") {
      degrees <- seq(60, 120, by = 10)
    } else if (slope_type == "decline") {
      degrees <- seq(105, 165, by = 10)
    } else {
      stop("slope_type must be one of incline, decline, or vertical.")
    }
  }

  theta <- degrees * (pi / 180)
  slope <- ifelse(theta == pi / 2, NA, tan(theta))
  alpha <- ifelse(is.na(slope), 1, abs(slope / sqrt(1 + slope^2)))
  beta <- ifelse(is.na(slope), 0, -sign(slope) / sqrt(1 + slope^2))

  split_mat <- matrix(nrow = length(theta), ncol = 7)
  colnames(split_mat) <-
    c("split", "score", "intercept", "slope", "alpha", "beta", "degree")
  for (i in seq_along(theta)) {
    z <- alpha[i] * x + beta[i] * y
    if (mechanism == "valley") {
      split_mat[i, 1:2] <- find_valley(z)
    } else if (mechanism == "boundary") {
      split_mat[i, 1:2] <- find_boundary(z)
    }
  }
  intercept <- split_mat[, 1] / beta

  split_mat[, "intercept"] <- intercept
  split_mat[, "slope"] <- slope
  split_mat[, "alpha"] <- alpha
  split_mat[, "beta"] <- beta
  split_mat[, "degree"] <- degrees

  max_ind <- which.max(split_mat[, 2])

  split_mat[max_ind, ]
}

#' Choose the slope and intercept for a two-dimensional valley split.
#'
#' @param x x axis data vector.
#' @param y y axis data vector.
#' @param slope_type `"incline"`, `"vertical"`, or `"decline"`.
#' @param degrees
#' The range of degrees to be searched over. Overrides `slope_type`.
#'
#' @returns A 1-row `data.frame` containing:
#' \describe{
#'   \item{`split`}{Valley split value.}
#'   \item{`score`}{Valley depth.}
#'   \item{`intercept`}{Intercept for two-dimensional line.}
#'   \item{`slope`}{Slope for two-dimensional line.}
#'   \item{`alpha`}{Mixing coefficient for x-axis variable.}
#'   \item{`beta`}{Mixing coefficient for y-axis variable.}
#'   \item{`degree`}{Degree of chosen two-dimensional line.}
#' }
#'
#' @export
find_valley_2d <- function(x, y, slope_type, degrees = NULL) {
  find_split_2d(x, y, slope_type, "valley", degrees)
}

#' Choose the slope and intercept for a two-dimensional boundary split.
#'
#' @inheritParams find_valley_2d
#' @inheritParams find_boundary
#'
#' @returns A 1-row `data.frame` containing:
#' \describe{
#'   \item{`split`}{Valley split value.}
#'   \item{`score`}{Valley depth.}
#'   \item{`intercept`}{Intercept for two-dimensional line.}
#'   \item{`slope`}{Slope for two-dimensional line.}
#'   \item{`alpha`}{Mixing coefficient for x-axis variable.}
#'   \item{`beta`}{Mixing coefficient for y-axis variable.}
#'   \item{`degree`}{Degree of chosen two-dimensional line.}
#' }
#'
#' @export
find_boundary_2d <- function(
    x, y, slope_type, degrees = NULL, boundary_noise_comp = TRUE) {
  find_split_2d(x, y, slope_type, "boundary", degrees, boundary_noise_comp)
}
