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

