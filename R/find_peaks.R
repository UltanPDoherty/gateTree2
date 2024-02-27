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
