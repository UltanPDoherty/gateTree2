find_valley <- function(dens, is_peak = NULL, score = FALSE, min_score = 0.01, min_height = 0.01) {

  if(is.null(is_peak)) {
    is_peak <- find_peaks(dens, min_height = min_height)
  }

  if(sum(is_peak) == 1) {
    if (score){
      return(c(NA, NA))
    } else {
      return(NA)
    }
  }

  maxpeak_ind <- which.max(dens$y)
  maxpeak <- data.frame(x = dens$x[maxpeak_ind],
                        y = dens$y[maxpeak_ind])

  otherpeaks_ind <- which(is_peak & dens$y != max(dens$y))
  otherpeaks <- data.frame(x = dens$x[otherpeaks_ind],
                           y = dens$y[otherpeaks_ind])

  scores <- valleys <- valley_ind <- c()
  interpeak <- matrix(nrow = length(dens$x), ncol = length(otherpeaks$x))
  for (i in seq_along(otherpeaks$x)) {
    peakpair_x <- sort(c(otherpeaks$x[i], maxpeak$x))
    peakpair_ind <- sort(c(otherpeaks_ind[i], maxpeak_ind))
    interpeak[, i] <- dens$x > peakpair_x[1] & dens$x < peakpair_x[2]
    valley_ind[i] <- peakpair_ind[1] + which.min(dens$y[interpeak[, i]])
    valleys[i] <- dens$y[valley_ind[i]]
    scores[i] <- otherpeaks$y[i] - valleys[i]
  }

  if (max(scores) < min_score) {
    if(score) {
      return(c(NA, NA))
    } else {
      return(NA)
    }
  } else {
    if (score) {
      return(c(dens$x[valley_ind[which.max(scores)]],max(scores)))
    } else {
      return(dens$x[valley_ind[which.max(scores)]])
    }
  }
}
