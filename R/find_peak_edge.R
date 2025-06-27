#' Find the border between a Gaussian component a Uniform noise component.
#'
#' @inheritParams find_valley
#' @param after_peak Whether the algorithm should look for an edge to the right
#'                   of the peak.
#'
#' @returns A scalar: the edge location.
#'
#' @export
find_peak_edge <- function(x, after_peak = TRUE) {
  gmm1 <- Mclust(x, G = 1, modelNames = "E", verbose = FALSE)

  dens1 <- dens(x, modelName = "X", parameters = gmm1$parameters)

  noise1 <- dens1 < hypvol(x, TRUE)

  gmm1n <- Mclust(
    x,
    G = 1, modelNames = "E", initialization = list("noise" = noise1),
    verbose = FALSE
  )

  if (after_peak) {
    edge <- max(x[gmm1n$classification == 1])
  } else {
    edge <- min(x[gmm1n$classification == 1])
  }

  return(edge)
}
