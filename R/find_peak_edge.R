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
  gmm1 <- densityMclust(
    x, G = 1, modelNames = "E", verbose = FALSE, plot = FALSE
  )

  dens1 <- gmm1$dens

  noise1 <- dens1 < hypvol(x, TRUE)

  gmm1n <- Mclust(
    x,
    G = 1, modelNames = "E", initialization = list("noise" = noise1),
    verbose = FALSE
  )

  if (after_peak) {
    after <- x > x[tail(which(dens1 == max(dens1)), 1)]
    edge <- min(x[after][noise1[after]])
  } else {
    before <- x < x[head(which(dens1 == max(dens1)), 1)]
    edge <- max(x[before][noise1[before]])
  }

  return(edge)
}
