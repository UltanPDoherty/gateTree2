#' Extract labels from `gatetree` output.
#'
#' @param gatetree_out Output from `gatetree` function.
#'
#' @returns List of label vectors for each sample.
#' @export
extract_labels <- function(gatetree_out) {
  bool_mats <- make_bool_mats(gatetree_out)

  label_list <- list()
  for (s in seq_along(bool_mats)) {
    pop_names <- colnames(bool_mats[[s]])
    pop_num <- ncol(bool_mats[[s]])
    obs_num <- nrow(bool_mats[[s]])
    label_list[[s]] <- rep("other", obs_num)
    for (p in seq_len(pop_num)) {
      label_list[[s]][bool_mats[[s]][, p]] <- pop_names[p]
    }
    label_list[[s]] <- as.factor(label_list[[s]])
  }
  names(label_list) <- names(bool_mats)

  label_list
}


make_bool_mats <- function(gatetree_out) {
  pop_num <- length(gatetree_out)
  samp_num <- length(gatetree_out[[1]]$subsetter)
  obs_num <- integer(samp_num)
  subset_num <- integer(samp_num)

  pop_names <- names(gatetree_out)

  bool_mats <- list()
  for (s in seq_len(samp_num)) {
    obs_num[s] <- nrow(gatetree_out[[1]]$subsetter[[s]])
    bool_mats[[s]] <- matrix(nrow = obs_num[s], ncol = pop_num)
    for (p in seq_len(pop_num)) {
      subset_num <- ncol(gatetree_out[[p]]$subsetter[[s]])
      bool_mats[[s]][, p] <- gatetree_out[[p]]$subsetter[[s]][, subset_num]
    }
    colnames(bool_mats[[s]]) <- pop_names

    bool_mats[[s]] <- merge_identical_columns(bool_mats[[s]])

    if (max(rowSums(bool_mats[[s]])) > 1) {
      stop("Populations are overlapping.")
    }
  }
  names(bool_mats) <- names(gatetree_out[[1]]$subsetter)

  bool_mats
}

merge_identical_columns <- function(mat) {
  col_num <- ncol(mat)
  col_names <- colnames(mat)

  if (col_num > 1) {
    for (i in seq(1, col_num - 1, by = 1)) {
      for (j in seq(2, col_num, by = 1)) {
        same_bool <- identical(mat[, i], mat[, j])
        if (same_bool) {
          colnames(mat)[i] <- paste0(col_names[i], "_", col_names[j])
          mat <- mat[, -j, drop = FALSE]

          merge_identical_columns(mat)
        }
      }
    }
  }

  mat
}
