# ==============================================================================

#' Check if final subsets are identical.
#'
#' @param subsetter The subsetting matrix.
#'
#' @return Vector indicating whether a cluster contains the same observations as
#' an existing cluster. If a pair of clusters are identical, the second one is
#' identified as a duplicate.
check_duplicates <- function(subsetter) {
  path_num <- ncol(subsetter)

  is_a_duplicate <- rep(FALSE, path_num)

  if (path_num > 1) {
    equal_subsets <- matrix(nrow = path_num, ncol = path_num)
    for (g in 1:(path_num - 1)) {
      for (h in (g + 1):path_num) {
        equal_subsets[g, h] <- all(subsetter[, g] == subsetter[, h])
        if (equal_subsets[g, h]) {
          print(paste0(
            "Failed to distinguish between populations ",
            g, " & ", h, "."
          ))
          is_a_duplicate[h] <- TRUE
        }
      }
    }
  }

  is_a_duplicate
}

# ==============================================================================

#' Find which observations are inside the cutoffs for each variable.
#'
#' @param x Matrix
#' @param min_val_cutoff Minimum value vector
#' @param max_val_cutoff Maximum value vector
#'
#' @return Matrix of the same dimensions as `x` containing logical values
#' indicating whether that observation (row) is inside the cutoffs for that
#' variable (column).
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

  inside_cutoffs
}
