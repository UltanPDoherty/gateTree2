#' Combine a pair of variables to find a line at the given angle.
#'
#' @inheritParams find_valley_2d
#' @param min_max_scaling Vector giving the new minimum and new maximum.
#'
#' @returns A numeric vector.
#' @export
combine_variables <- function(
    x, y, degrees, min_max_scaling = NULL) {
  stopifnot(length(degrees) == 1)
  stopifnot(length(x) == length(y))

  theta <- degrees * (pi / 180)
  slope <- ifelse(theta == pi / 2, NA, tan(theta))
  alpha <- ifelse(is.na(slope), 1, abs(slope / sqrt(1 + slope^2)))
  beta <- ifelse(is.na(slope), 0, -sign(slope) / sqrt(1 + slope^2))

  z <- alpha * x + beta * y

  if (!is.null(min_max_scaling)) {
    z <- z - min(z) + min_max_scaling[1]
    z <- min_max_scaling[2] * (z / max(z))
  }

  z
}
