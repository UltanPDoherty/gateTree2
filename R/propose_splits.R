#' @title Propose a density valley or GMM boundary.
#'
#' @description
#' Find a density valley or GMM boundary given a subset of cells and a set of
#' splittable variables.
#'
#' @inheritParams gatetree
#' @param subsetter_g Column g of the subsetting matrix.
#' @param splittable_vars_g Row g of the splittable_vars matrix.
#'
#' @return List
#' * splits: the split locations for each variable.
#' * scores: the valley depth percentage or scaled BIC difference of that split.
#' * scenario: a character string indicating the type of the split (`"valley"`
#' or `"boundary"`), if no split was found (`"nothing"`), or if the split was
#' identified in the context of the `explore` option (`"explore"`).
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
    found_valley <- any(!is.na(proposals$splits))
    found_boundary <- FALSE

    if (!found_valley && use_boundaries) {
      proposals <- propose_boundaries(
        x, min_scaled_bic_diff,
        subsetter_g, splittable_vars_g
      )

      found_boundary <- any(!is.na(proposals$splits))
    }
  }

  scenario <- ifelse(found_valley,
    "valley",
    ifelse(found_boundary, "boundary", "nothing")
  )

  return(list(
    splits = proposals$splits,
    scores = proposals$scores,
    scenario = scenario
  ))
}
