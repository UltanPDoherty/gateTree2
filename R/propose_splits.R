
#===============================================================================

#' Wrapper for `find_valley`.
#'
#' @param x Data.
#' @param subsetter_g Column g of the subsetting matrix.
#' @param splittable_vars_g Row g of the splittable_vars matrix.
#' @param min_depth `min_depth` for `find_valley` function.
#' @param min_height `min_height` for `find_valley` function.
#'
#' @return valleys
#' @export
propose_valleys <- function(
    x,
    subsetter_g,
    splittable_vars_g = rep(TRUE, ncol(x)),
    min_depth, min_height
) {
  valleys <- matrix(nrow = 2, ncol = ncol(x))

  # loop over all variables to propose splits
  for (p in which(splittable_vars_g)){
    scale01_gp <- scale01(x[subsetter_g, p])
    dens01_gp <- stats::density(scale01_gp$y)

    valleys[, p] <- find_valley(
      dens01_gp,
      depth = TRUE,
      min_depth = min_depth,
      min_height = min_height
    )

    valleys[1, p] <- unscale01(valleys[1, p], scale01_gp$min, scale01_gp$max)
  }

  return(valleys)
}

#===============================================================================

#' Wrapper for `find_valley`.
#'
#' @param x Data.
#' @param subsetter_g Column g of the subsetting matrix.
#' @param splittable_vars_g Row g of the splittable_vars matrix.
#'
#' @return valleys
#' @export
propose_boundaries <- function(
    x,
    subsetter_g,
    splittable_vars_g = rep(TRUE, ncol(x))
) {
  boundaries <- matrix(nrow = 2, ncol = ncol(x))

  # loop over all variables to propose splits
  for (p in which(splittable_vars_g)){
    boundaries[, p] <- find_boundary(x[subsetter_g, p])
  }

  return(boundaries)
}

