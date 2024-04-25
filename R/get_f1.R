#' @title Compute the F1 measure of a clustering solution.
#'
#' @description
#' Compare a vector of cluster labels to another vector of class labels using
#' the F1 measure with cluster - class matching carried out using the Hungarian
#' algorithm. Particular classes or clusters may be excluded from this matching.
#' Observations of a certain class may also be removed completely.
#'
#' @param class_labels Vector of true / reference labels to compare to.
#' @param clust_labels Vector of cluster labels to evaluate.
#' @param remove_class Vector of class labels to be removed from the data.
#' @param no_match_class Vector of class labels to which cluster labels should
#'                       not be matched.
#' @param no_match_clust Vector of cluster labels which should not be matched
#'                       to class labels.
#' @param prec_rec Logical: whether precision and recall vectors are returned.
#'
#' @return List:
#' * f1_mat: matrix of pairwise F1 measure values for every class and cluster
#' with matched values on the diagonal.
#' * f1_vec: vector of matched F1 measure values for each class.
#' * pr_vec (if `prec_rec = TRUE`): vector of matched precision values for each
#' class.
#' * re_vec (if `prec_rec = TRUE`): vector of matched recall values for each
#' class.
#' @export
get_f1 <- function(
  class_labels,
  clust_labels,
  remove_class = NULL,
  no_match_class = NULL,
  no_match_clust = NULL,
  prec_rec = FALSE
) {
  clust_labels[is.na(clust_labels)] <- max(clust_labels, na.rm = TRUE) + 1
  class_labels[is.na(class_labels)] <- max(class_labels, na.rm = TRUE) + 1

  not_removed  <- !(class_labels %in% remove_class)
  class_labels <- class_labels[not_removed]
  clust_labels <- clust_labels[not_removed]

  class_num <- length(unique(class_labels))
  clust_num <- length(unique(clust_labels))

  cross_tab <- table(class_labels, clust_labels)
  class_tab <- table(class_labels)
  clust_tab <- table(clust_labels)

  re_mat <- pr_mat <- f1_mat <- matrix(nrow = class_num, ncol = clust_num,
                                       dimnames = dimnames(cross_tab))
  for (i in seq_len(class_num)) {
    for (j in seq_len(clust_num)) {
      re_mat[i, j] <- cross_tab[i, j] / class_tab[i]
      pr_mat[i, j] <- cross_tab[i, j] / clust_tab[j]

      f1_mat[i, j] <- ifelse(
        cross_tab[i, j] != 0,
        2 * re_mat[i, j] * pr_mat[i, j] / (re_mat[i, j] + pr_mat[i, j]),
        0
      )
    }
  }

  matchable_classes <- setdiff(unique(class_labels), no_match_class)
  matchable_class_rows <- which(rownames(f1_mat) %in% matchable_classes)
  match_class_num <- length(matchable_classes)
  class_to_clust <- rep(NA, class_num)

  matchable_clusts <- setdiff(unique(clust_labels), no_match_clust)
  matchable_clust_rows <- which(colnames(f1_mat) %in% matchable_clusts)
  match_clust_num <- length(matchable_clusts)

  if (match_class_num <= match_clust_num) {
    class_to_clust_match <- as.numeric(
      clue::solve_LSAP(f1_mat[matchable_class_rows, matchable_clust_rows],
                       maximum = TRUE)
    )
    # clust_to_class gives column numbers of the selected f1_mat submatrix
    class_to_clust[matchable_class_rows] <- class_to_clust_match
  } else {
    clust_to_class <- as.numeric(
      clue::solve_LSAP(t(f1_mat[matchable_class_rows, matchable_clust_rows]),
                       maximum = TRUE)
    )
    # clust_to_class gives row numbers of the selected f1_mat submatrix
    other_class    <- setdiff(matchable_class_rows, clust_to_class)
    clust_to_class <- c(clust_to_class, other_class)
    class_to_clust[matchable_class_rows] <- order(clust_to_class)
    class_to_clust[!is.na(class_to_clust)
                   & class_to_clust > match_clust_num] <- NA
    # match_clust_num is the number of columns of the selected f1_mat submatrix
  }

  # we need to convert from column numbers of the f1_mat submatrix to column
  # numbers of the full f1_mat matrix
  class_to_clust2 <- (1:clust_num)[matchable_clust_rows][class_to_clust]

  pr_mat_ord <- pr_mat[, class_to_clust2]
  pr_vec <- diag(pr_mat_ord)

  re_mat_ord <- re_mat[, class_to_clust2]
  re_vec <- diag(re_mat_ord)

  f1_mat_ord <- f1_mat[, class_to_clust2]
  f1_vec <- diag(f1_mat_ord)

  out <- list(f1_mat = f1_mat_ord,
              f1_vec = f1_vec)

  if (prec_rec) {
    out$pr_vec <- pr_vec
    out$re_vec <- re_vec
  }

  return(out)
}
