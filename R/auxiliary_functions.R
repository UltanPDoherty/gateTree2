#===============================================================================

#' 0-1 scaling.
#'
#' @param x Data.
#' @param other_min Minimum to be used for 0-1 scaling.
#' @param other_max Maximum to be used for 0-1 scaling.
#'
#' @return Scaled version of `x`.
#' @export
scale01 <- function(x, other_min = NULL, other_max = NULL) {

  if (is.null(other_min)) {
    minimum <- min(x)
  } else {
    minimum <- other_min
  }

  if (is.null(other_max)) {
    maximum <- max(x)
  } else {
    maximum <- other_max
  }

  stopifnot("Maximum and minimum are equal" = maximum != minimum)
  y <- (x - minimum) / (maximum - minimum)

  return(list(y = y, min = minimum, max = maximum))
}

#===============================================================================

#' Undo 0-1 scaling.
#'
#' @param x Data.
#' @param unscaled_min Minimum used for 0-1 scaling.
#' @param unscaled_max Maximum used for 0-1 scaling.
#'
#' @return Unscaled version of `x`.
#' @export
unscale01 <- function(x, unscaled_min, unscaled_max) {
  return(x * (unscaled_max - unscaled_min) + unscaled_min)
}

#===============================================================================

#' Check if final subsets are identical.
#'
#' @param subsetter The subsetting matrix.
#'
#' @return `is_a_duplicate`.
#' @export
check_duplicates <- function(subsetter) {

  path_num <- ncol(subsetter)

  is_a_duplicate <- rep(FALSE, path_num)

  if (path_num > 1) {
    equal_subsets <- matrix(nrow = path_num, ncol = path_num)
    for (g in 1:(path_num - 1)) {
      for (h in (g + 1):path_num) {
        equal_subsets[g, h] <- all(subsetter[, g] == subsetter[, h])
        if (equal_subsets[g, h]) {
          print(paste0("Failed to distinguish between populations ",
                       g, " & ", h, "."))
          is_a_duplicate[h] <- TRUE
        }
      }
    }
  }

  return(is_a_duplicate)
}

#===============================================================================

#' Find which observations are inside the cutoffs for each variable.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#'
#' @return inside_cutoffs
#' @export
find_inside_cutoffs <- function(x, min_val_cutoff, max_val_cutoff) {
  # find which observations are outside either of the cutoffs for each variable
  if (is.null(min_val_cutoff)) {
    below_cutoff <- array(FALSE, dim = dim(x))
  } else {
    below_cutoff <- array(dim = dim(x))
    for (j in seq_len(ncol(x))) {
      below_cutoff[, j] <- x[, j] <= min_val_cutoff[j]
    }
  }
  if (is.null(max_val_cutoff)) {
    above_cutoff <- array(FALSE, dim = dim(x))
  } else {
    above_cutoff <- array(dim = dim(x))
    for (j in seq_len(ncol(x))) {
      above_cutoff[, j] <- x[, j] >= max_val_cutoff[j]
    }
  }
  # inside_cutoffs is a matrix of the same dimensions as the data
  # an observation-variable pair is TRUE if it is inside the two cutoffs
  inside_cutoffs <- !below_cutoff & !above_cutoff

  return(inside_cutoffs)
}
