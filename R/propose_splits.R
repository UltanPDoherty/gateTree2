#' Title
#'
#' @inheritParams gatetree
#' @param subsetter_g Column g of the subsetting matrix.
#' @param splittable_vars_g Row g of the splittable_vars matrix.
#'
#' @return List: matrix, scenario
#' @export
propose_splits <- function(x, subsetter_g, splittable_vars_g,
                           min_size, min_depth, min_height = min_depth,
                           min_scaled_bic_diff = 0,
                           use_boundaries) {
  if (sum(subsetter_g) < min_size) {
    found_valley <- FALSE
    found_boundary <- FALSE
    proposals <- NA
  } else {
    proposals <- propose_valleys(
      x, subsetter_g, splittable_vars_g,
      min_depth, min_height
    )
    found_valley <- any(!is.na(proposals[1, ]))
    found_boundary <- FALSE

    if (!found_valley && use_boundaries) {
      proposals <- propose_boundaries(x, min_scaled_bic_diff,
                                      subsetter_g, splittable_vars_g)

      found_boundary <- any(!is.na(proposals[1, ]))
    }
  }

  scenario <- ifelse(found_valley,
                     "valley",
                     ifelse(found_boundary, "boundary", "nothing"))

  return(list(matrix = proposals,
              scenario = scenario))
}

#===============================================================================

#' Wrapper for `find_valley`.
#'
#' @inheritParams propose_splits
#'
#' @return valleys
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
      return_depth = TRUE,
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
#' @inheritParams propose_splits
#'
#' @return valleys
propose_boundaries <- function(
  x,
  min_scaled_bic_diff = 0,
  subsetter_g,
  splittable_vars_g = rep(TRUE, ncol(x))
) {
  boundaries <- matrix(nrow = 2, ncol = ncol(x))

  # loop over all variables to propose splits
  for (p in which(splittable_vars_g)){
    boundaries[, p] <- find_boundary(x[subsetter_g, p], min_scaled_bic_diff)
  }

  return(boundaries)
}
