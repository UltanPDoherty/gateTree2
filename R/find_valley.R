#' Title
#'
#' @param x Data vector
#'
#' @returns Vector of length two: valley location and depth
find_valley <- function(x) {
  nbin <- 100
  x_ash <- ash::bin1(x, nbin = nbin)
  counts <- x_ash$nc
  breaks <- seq(x_ash$ab[1], x_ash$ab[2], length.out = nbin + 1)
  midpoints <- (breaks[1:nbin] + breaks[-1]) / 2
  
  w <- 1
  is_peak <- peak_left <- peak_right <- rep(FALSE, nbin)
  for (i in seq(w + 1, nbin - w)) {
    peak_left[i] <- all(counts[i] > counts[(i - w):(i - 1)])
    peak_right[i] <- all(counts[i] > counts[(i + 1):(i + w)])
    is_peak[i] <- peak_left[i] & peak_right[i]
  }
  is_peak
  
  if (sum(is_peak) == 1) {
    best_valley <- NA
    best_depth <- NA
  } else {
    maxpeak_ind <- which.max(counts)
    maxpeak <- data.frame(x = midpoints[maxpeak_ind], y = counts[maxpeak_ind])
    
    otherpeaks_ind <- which(is_peak & counts != max(counts))
    otherpeaks <- data.frame(
      x = midpoints[otherpeaks_ind],
      y = counts[otherpeaks_ind]
    )
    
    depths <- valleys <- valley_ind <- c()
    interpeak <- matrix(nrow = nbin, ncol = length(otherpeaks$x))
    for (i in seq_along(otherpeaks$x)) {
      peakpair_x <- sort(c(otherpeaks$x[i], maxpeak$x))
      peakpair_ind <- sort(c(otherpeaks_ind[i], maxpeak_ind))
      interpeak[, i] <- midpoints > peakpair_x[1] & midpoints < peakpair_x[2]
      valley_ind[i] <- peakpair_ind[1] + which.min(counts[interpeak[, i]])
      valleys[i] <- counts[valley_ind[i]]
      depths[i] <- otherpeaks$y[i] - valleys[i]
    }
    
    if (is.null(depths)) {
      best_depth <- NA
      best_valley <- NA
    } else {
      best_depth <- max(depths)
      best_valley <- midpoints[valley_ind[which.max(depths)]]
    }
  }
  
  return(c(best_valley, best_depth))
  
}
