#' Compute f1 given two vectors of labels
#'
#' @param clust_labels Vector of cluster labels to evaluate.
#' @param true_labels Vector of true / reference labels to compare to.
#' @param exclude_from_true Vector of true labels to be excluded.
#'                          Ignored if NULL.
#' @param prec_rec Logical: Should vectors of precision and recall values,
#'                 with assignment based on f1, be returned?
#'
#' @return List of mean f1, mean precision, mean recall, and f1 matrix.
#' @export

compute_f1 <- function(clust_labels,
                       true_labels,
                       exclude_from_true = NULL,
                       do_not_assign_to = NULL,
                       prec_rec = FALSE) {

  clust_nas <- is.na(clust_labels)
  clust_labels[clust_nas] <- max(clust_labels[!clust_nas]) + 1

  if (!is.null(exclude_from_true)) {
    excluded     <- true_labels %in% exclude_from_true
    true_labels  <- true_labels[!excluded]
    clust_labels <- clust_labels[!excluded]

    if (is.factor(true_labels)) {
      true_labels <- droplevels(true_labels)
    }

    if (is.factor(clust_labels)) {
      clust_labels <- droplevels(clust_labels)
    }
  }

  true_num  <- length(unique(true_labels))
  clust_num <- length(unique(clust_labels))

  tab <- table(true_labels, clust_labels)

  re_mat <- sweep(tab, STATS = rowSums(tab), FUN = "/", MARGIN = 1)
  pr_mat <- sweep(tab, STATS = colSums(tab), FUN = "/", MARGIN = 2)

  f1_mat <- 2 * re_mat * pr_mat / (re_mat + pr_mat)
  f1_mat[f1_mat == "NaN"] <- 0

  if (true_num <= clust_num) {
    true_to_clust <- as.numeric(
      clue::solve_LSAP(f1_mat[1:true_num, ], maximum = TRUE)
    )
    other_clust   <- setdiff(1:clust_num, true_to_clust)
    true_to_clust <- c(true_to_clust, other_clust)
  } else {
    clust_to_true <- as.numeric(
      clue::solve_LSAP(t(f1_mat[1:true_num, ]), maximum = TRUE)
    )
    other_true    <- setdiff(1:true_num, clust_to_true)
    true_to_clust <- (1:true_num)[order(c(clust_to_true, other_true))]
    true_to_clust[true_to_clust > clust_num] <- NA
  }

  f1_mat2 <- f1_mat[, true_to_clust]
  f1_mat2[, is.na(true_to_clust)] <- 0

  f1_vec <- pr_vec <- re_vec <- rep(NA, true_num)
  for (i in 1:true_num) {
    if (!is.na(true_to_clust)[i]) {
      f1_vec[i] <- f1_mat[i, true_to_clust[i]]
      pr_vec[i] <- pr_mat[i, true_to_clust[i]]
      re_vec[i] <- re_mat[i, true_to_clust[i]]
    } else {
      f1_vec[i] <- 0
      pr_vec[i] <- 0
      re_vec[i] <- 0
    }
  }

  vec_sums <- c(sum(f1_vec), sum(pr_vec), sum(re_vec))
  matched_vals <- vec_sums / sum(!is.na(true_to_clust))
  names(matched_vals) <- c("f1", "pr", "re")

  out <- list(f1_mat = f1_mat2,
              f1_vec = f1_vec)

  if (prec_rec) {
    out$pr_vec <- pr_vec
    out$re_vec <- re_vec
  }

  return(out)
}
