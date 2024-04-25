find_peaks <- function(dens, width_percent = 0.01, min_height = 0.01) {

  dens$y <- dens$y / max(dens$y) * 100

  # w is the number of density points in width_percent of the range
  w <- round(length(dens$x) * width_percent)

  is_peak <- peak_left <- peak_right <- c()
  for (i in seq(w + 1, length(dens$y) - w)) {
    peak_left[i - w] <- all(dens$y[i] > dens$y[(i - w - 1):(i - 1)])
    peak_right[i - w] <- all(dens$y[i] > dens$y[(i + 1):(i + w)])
    is_peak[i - w] <- peak_left[i - w] & peak_right[i - w]
  }

  is_peak <- c(rep(FALSE, w), is_peak, rep(FALSE, w)) & dens$y >= min_height

  return(is_peak)
}

#===============================================================================

#' @title Find the optimal density valley.
#'
#' @description
#' Find the optimal kernel density valley and check if its depth percentage is
#' greater than `min_depth`.
#'
#' @inheritParams gatetree
#' @param dens Output from `stats::density` function.
#'
#' @return Vector consisting of the valley and its depth percentage. The
#' valley will be NA if its depth percentage is less than `min_depth`.
find_valley <- function(dens, min_depth = 0.01, min_height = 0.01) {

  is_peak <- find_peaks(dens, min_height = min_height)

  if (sum(is_peak) == 1) {
    return(c(NA, NA))
  }

  maxpeak_ind <- which.max(dens$y)
  maxpeak <- data.frame(x = dens$x[maxpeak_ind],
                        y = dens$y[maxpeak_ind])

  dens$y <- dens$y / maxpeak$y * 100
  maxpeak$y <- 100.0

  otherpeaks_ind <- which(is_peak & dens$y != max(dens$y))
  otherpeaks <- data.frame(x = dens$x[otherpeaks_ind],
                           y = dens$y[otherpeaks_ind])

  depths <- valleys <- valley_ind <- c()
  interpeak <- matrix(nrow = length(dens$x), ncol = length(otherpeaks$x))
  for (i in seq_along(otherpeaks$x)) {
    peakpair_x <- sort(c(otherpeaks$x[i], maxpeak$x))
    peakpair_ind <- sort(c(otherpeaks_ind[i], maxpeak_ind))
    interpeak[, i] <- dens$x > peakpair_x[1] & dens$x < peakpair_x[2]
    valley_ind[i] <- peakpair_ind[1] + which.min(dens$y[interpeak[, i]])
    valleys[i] <- dens$y[valley_ind[i]]
    depths[i] <- otherpeaks$y[i] - valleys[i]
  }

  if (max(depths) < min_depth) {
    return(c(NA, NA))
  } else {
    return(c(dens$x[valley_ind[which.max(depths)]], max(depths)))
  }
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
      min_depth = min_depth,
      min_height = min_height
    )

    valleys[1, p] <- unscale01(valleys[1, p], scale01_gp$min, scale01_gp$max)
  }

  return(valleys)
}
