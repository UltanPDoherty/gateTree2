#' Find the deepest univariate density valley.
#'
#' @inheritParams gatetree
#' @param x Data vector
#'
#' @returns Vector of length two: valley location and depth
#'
#' @export
find_valley <- function(x, min_kde_size = 100) {
  if (length(x) < min_kde_size) {
    return(c(NA, NA))
  }

  dens <- stats::density(x)
  locations <- dens$x
  dens01 <- dens$y / max(dens$y)

  location_num <- length(locations)
  w <- round(location_num / 100)

  is_peak <- peak_left <- peak_right <- rep(FALSE, location_num)
  for (i in seq(1, w)) {
    peak_left[i] <- all(dens01[i] >= dens01[1:i])
    peak_right[i] <- all(dens01[i] >= dens01[(i + 1):(i + w)])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  for (i in seq(w + 1, location_num - w)) {
    peak_left[i] <- all(dens01[i] >= dens01[(i - w):(i - 1)])
    peak_right[i] <- all(dens01[i] >= dens01[(i + 1):(i + w)])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  for (i in seq(location_num - w + 1, location_num)) {
    peak_left[i] <- all(dens01[i] >= dens01[(i - w):(i - 1)])
    peak_right[i] <- all(dens01[i] >= dens01[i:location_num])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  is_peak

  if (sum(is_peak) == 1) {
    best_valley <- NA
    best_depth <- NA
  } else {
    maxpeak_ind <- which.max(dens01)
    maxpeak <- data.frame(x = locations[maxpeak_ind], y = dens01[maxpeak_ind])

    otherpeaks_ind <- which(is_peak & dens01 != max(dens01))
    otherpeaks <- data.frame(
      x = locations[otherpeaks_ind],
      y = dens01[otherpeaks_ind]
    )

    depths <- valleys <- valley_ind <- c()
    interpeak <- matrix(nrow = location_num, ncol = length(otherpeaks$x))
    for (i in seq_along(otherpeaks$x)) {
      peakpair_x <- sort(c(otherpeaks$x[i], maxpeak$x))
      peakpair_ind <- sort(c(otherpeaks_ind[i], maxpeak_ind))
      interpeak[, i] <- locations > peakpair_x[1] & locations < peakpair_x[2]
      valley_ind[i] <- peakpair_ind[1] + which.min(dens01[interpeak[, i]])
      valleys[i] <- dens01[valley_ind[i]]
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
