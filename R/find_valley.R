#' Find the deepest univariate density valley.
#'
#' @param x Data vector
#'
#' @returns Vector of length two: valley location and depth
#'
#' @export
find_valley <- function(x) {
  if (length(x) < 10) {
    return(c(NA, NA))
  }

  dens <- stats::density(x)
  locations <- dens$x
  counts <- dens$y * length(x)

  location_num <- length(locations)
  w <- round(location_num / 100)

  is_peak <- peak_left <- peak_right <- rep(FALSE, location_num)
  for (i in seq(1, w)) {
    peak_left[i] <- all(counts[i] >= counts[1:i])
    peak_right[i] <- all(counts[i] >= counts[(i + 1):(i + w)])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  for (i in seq(w + 1, location_num - w)) {
    peak_left[i] <- all(counts[i] >= counts[(i - w):(i - 1)])
    peak_right[i] <- all(counts[i] >= counts[(i + 1):(i + w)])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  for (i in seq(location_num - w + 1, location_num)) {
    peak_left[i] <- all(counts[i] >= counts[(i - w):(i - 1)])
    peak_right[i] <- all(counts[i] >= counts[i:location_num])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  is_peak

  if (sum(is_peak) == 1) {
    best_valley <- NA
    best_depth <- NA
  } else {
    maxpeak_ind <- which.max(counts)
    maxpeak <- data.frame(x = locations[maxpeak_ind], y = counts[maxpeak_ind])

    otherpeaks_ind <- which(is_peak & counts != max(counts))
    otherpeaks <- data.frame(
      x = locations[otherpeaks_ind],
      y = counts[otherpeaks_ind]
    )

    depths <- valleys <- valley_ind <- c()
    interpeak <- matrix(nrow = location_num, ncol = length(otherpeaks$x))
    for (i in seq_along(otherpeaks$x)) {
      peakpair_x <- sort(c(otherpeaks$x[i], maxpeak$x))
      peakpair_ind <- sort(c(otherpeaks_ind[i], maxpeak_ind))
      interpeak[, i] <- locations > peakpair_x[1] & locations < peakpair_x[2]
      valley_ind[i] <- peakpair_ind[1] + which.min(counts[interpeak[, i]])
      valleys[i] <- counts[valley_ind[i]]
      depths[i] <- otherpeaks$y[i] - valleys[i]
    }

    if (is.null(depths)) {
      best_depth <- NA
      best_valley <- NA
    } else {
      best_depth <- max(depths)
      best_valley <- locations[valley_ind[which.max(depths)]]
    }
  }

  c(best_valley, best_depth)
}
