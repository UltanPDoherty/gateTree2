#' Extract labels from `gatetree` output.
#'
#' @param gatetree_out Output from `gatetree` function.
#'
#' @returns List of label vectors for each sample.
#' @export
extract_labels <- function(gatetree_out) {
  if (!is.null(gatetree_out$output)) {
    gatetree_out <- gatetree_out$output
  }

  bool_mats <- make_bool_mats(gatetree_out)

  label_list <- list()
  for (b in seq_along(bool_mats)) {
    label_list[[b]] <- list()
    for (s in seq_along(bool_mats[[b]])) {
      pop_names <- colnames(bool_mats[[b]][[s]])
      pop_num <- ncol(bool_mats[[b]][[s]])
      obs_num <- nrow(bool_mats[[b]][[s]])
      label_list[[b]][[s]] <- rep("other", obs_num)
      for (p in seq_len(pop_num)) {
        label_list[[b]][[s]][bool_mats[[b]][[s]][, p]] <- pop_names[p]
      }
      label_list[[b]][[s]] <- as.factor(label_list[[b]][[s]])
    }
    names(label_list[[b]]) <- names(bool_mats[[b]])
  }
  names(label_list) <- names(bool_mats)

  label_list
}

make_bool_mats <- function(gatetree_out) {
  pop_num <- length(gatetree_out)
  batch_num <- length(gatetree_out[[1]]$subsetter)
  samp_num <- integer(batch_num)

  obs_num <- list()
  subset_num <- list()

  pop_names <- names(gatetree_out)

  bool_mats <- list()
  for (b in seq_len(batch_num)) {
    samp_num[b] <- length(gatetree_out[[1]]$subsetter[[b]])
    obs_num[[b]] <- integer(samp_num[b])
    subset_num[[b]] <- integer(samp_num[b])
    bool_mats[[b]] <- list()
    for (s in seq_len(samp_num[b])) {
      obs_num[[b]][s] <- nrow(gatetree_out[[1]]$subsetter[[b]][[s]])
      bool_mats[[b]][[s]] <- matrix(nrow = obs_num[[b]][s], ncol = pop_num)
      for (p in seq_len(pop_num)) {
        subset_num[[b]][s] <- ncol(gatetree_out[[p]]$subsetter[[b]][[s]])
        bool_mats[[b]][[s]][, p] <-
          gatetree_out[[p]]$subsetter[[b]][[s]][, subset_num[[b]][s]]
      }
      colnames(bool_mats[[b]][[s]]) <- pop_names

      bool_mats[[b]][[s]] <- merge_identical_columns(bool_mats[[b]][[s]])

      if (max(rowSums(bool_mats[[b]][[s]])) > 1) {
        stop("Populations are overlapping.")
      }
    }
    names(bool_mats[[b]]) <- names(gatetree_out[[1]]$subsetter[[b]])
  }
  names(bool_mats) <- names(gatetree_out[[1]]$subsetter)

  bool_mats
}

merge_identical_columns <- function(mat) {
  if (ncol(mat) > 1) {
    i <- 1
    while (i < ncol(mat)) {
      j <- i + 1
      while (j <= ncol(mat)) {
        same_bool <- identical(mat[, i], mat[, j])
        if (same_bool) {
          colnames(mat)[i] <- paste0(colnames(mat)[i], "_", colnames(mat)[j])
          mat <- mat[, -j, drop = FALSE]

          merge_identical_columns(mat)
        }
        j <- j + 1
      }
      i <- i + 1
    }
  }

  mat
}
