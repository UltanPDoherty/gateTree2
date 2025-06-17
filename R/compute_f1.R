compute_f1 <- function(class_labels, cluster_labels) {
  class_num <- length(unique(class_labels))
  clust_num <- length(unique(cluster_labels))

  class_count <- table(class_labels)
  clust_count <- table(cluster_labels)
  contingency <- table(class_labels, cluster_labels)

  precision <- recall <- f1 <- matrix(nrow = class_num, ncol = clust_num)
  for (i in seq_len(class_num)) {
    for (j in seq_len(clust_num)) {
      precision[i, j] <- contingency[i, j] / clust_count[j]
      recall[i, j] <- contingency[i, j] / class_count[i]
      f1[i, j] <- 2 * contingency[i, j] / (class_count[i] + clust_count[j])
    }
  }
  rownames(precision) <- rownames(recall) <- rownames(f1) <-
    rownames(contingency)
  colnames(precision) <- colnames(recall) <- colnames(f1) <-
    colnames(contingency)

  list(
    "class_num" = class_num,
    "clust_num" = clust_num,
    "class_count" = class_count,
    "clust_count" = clust_count,
    "contingency" = contingency,
    "precision" = precision,
    "recall" = recall,
    "f1" = f1
  )
}

#' Compute the maximum F1 for a single class with respect to a single cluster.
#'
#' @description
#' Compute a matrix of F1 scores for each class-cluster pair, then optionally
#' obtain maximum F1 values for each class, either independently or according
#' to a one-to-one mapping by the Hungarian algorithm. Finally, a weighted or
#' unweighted mean may be computed across the classes.
#'
#' @param class_labels Vector of true class labels.
#' @param cluster_labels Vector of estimated cluster labels.
#' @param weighted Logical: whether the mean should be weighted with respect to
#'                 the size of the classes.
#' @param hungarian Logical: whether the Hungarian algorithm should be used to
#'                  obtain a one-to-one mapping of classes to clusters.
#' @param return_option Character: whether to return the full pairwise F1
#'                      matrix ("matrix"), the vector of maximum F1 scores for
#'                      each class ("vector", default), or the mean ("value").
#' @param no_match_class Classes to be excluded when matching. Their F1 values
#'                       are set to -Inf.
#' @param no_match_cluster Clusters to be excluded when matching. Their F1
#'                         values are set to -Inf.
#'
#' @returns
#' Either a matrix, vector, or mean value, depending on `return_option`.
#'
#' @export
compute_max_f1 <- function(
    class_labels, cluster_labels,
    weighted = TRUE, hungarian = FALSE,
    return_option = c("vector", "value", "matrix"),
    no_match_class = NULL, no_match_cluster = NULL) {
  return_option <- rlang::arg_match(
    return_option, c("value", "vector", "matrix")
  )

  f1_output <- compute_f1(class_labels, cluster_labels)
  class_num <- f1_output$class_num
  clust_num <- f1_output$clust_num
  class_count <- f1_output$class_count
  f1 <- f1_output$f1

  if (!is.null(no_match_class)) {
    bool_no_match_class <- rownames(f1) %in% no_match_class
    if (any(bool_no_match_class)) {
      f1 <- f1[-which(bool_no_match_class), , drop = FALSE]
      class_count <- class_count[-which(bool_no_match_class)]
      class_num <- class_num - sum(bool_no_match_class)
    }
  }
  if (!is.null(no_match_cluster)) {
    bool_no_match_clust <- rownames(f1) %in% no_match_cluster
    if (any(bool_no_match_clust)) {
      f1 <- f1[, -which(bool_no_match_clust), drop = FALSE]
      clust_num <- clust_num - sum(bool_no_match_clust)
    }
  }

  if (clust_num < class_num) {
    f1 <- cbind(
      f1, matrix(0, nrow = class_num, ncol = class_num - clust_num)
    )
  }

  if (hungarian) {
    hungarian_order <- clue::solve_LSAP(f1, maximum = TRUE)

    f1 <- f1[, hungarian_order]

    f1_vec <- diag(f1)
  } else {
    f1_vec <- apply(f1, 1, max)
  }

  if (weighted) {
    f1_val <- stats::weighted.mean(f1_vec, class_count, na.rm = TRUE)
  } else {
    f1_val <- mean(f1_vec, na.rm = TRUE)
  }

  if (return_option == "matrix") {
    return(f1)
  } else if (return_option == "vector") {
    return(f1_vec)
  } else if (return_option == "value") {
    return(f1_val)
  } else {
    stop("Invalid return_option")
  }
}

#' Compute the F1 for each class of the optimal merged set of cluster.
#'
#' @description
#' For each class, start with the highest F1 cluster, then repeatedly merge it
#' with whichever cluster results in the best F1 score. Retrospectively identify
#' when the F1 of the merged cluster with respect to that class is maximised.
#'
#' @inheritParams compute_max_f1
#' @param merge_limit Maximum number of cluster merges to perform per class.
#'
#' @returns `compute_merged_f1` returns a list with the following elements:
#' \describe{
#'  \item{`f1`}{Vector: max merged F1 values for each class.}
#'  \item{`merges`}{Vector: number of merges performed for each class.}
#' }
#'
#' @export
compute_merged_f1 <- function(
    class_labels, cluster_labels, merge_limit = NULL) {
  f1_output <- compute_f1(class_labels, cluster_labels)
  class_num <- f1_output$class_num
  clust_num <- f1_output$clust_num
  class_count <- f1_output$class_count
  clust_count <- f1_output$clust_count
  contingency <- f1_output$contingency
  f1 <- f1_output$f1

  if (is.null(merge_limit)) {
    merge_limit <- clust_num - 1
  } else {
    merge_limit <- min(clust_num - 1, merge_limit)
  }

  merging_f1 <- matrix(nrow = class_num, ncol = clust_num)
  for (i in seq_len(class_num)) {
    f1_order <- order(f1[i, ], decreasing = TRUE)
    cont_i <- contingency[i, f1_order]
    clust_count_i <- clust_count[f1_order]

    merging_f1[i, 1] <- max(f1[i, ])

    merge_count <- 1
    while (merge_count <= merge_limit) {
      temp_merging_f1 <- double(clust_num - merge_count)

      for (j in seq_len(clust_num - merge_count)) {
        x <- cont_i[1] + cont_i[1 + j]
        n <- clust_count_i[1] + clust_count_i[1 + j]
        temp_merging_f1[j] <- 2 * x / (class_count[i] + n)
      }
      merge_choice <- which.max(temp_merging_f1) + 1

      cont_i[1] <- cont_i[1] + cont_i[merge_choice]
      cont_i <- cont_i[-merge_choice]
      clust_count_i[1] <- clust_count_i[1] + clust_count_i[merge_choice]
      clust_count_i <- clust_count_i[-merge_choice]

      merging_f1[i, 1 + merge_count] <-
        2 * cont_i[1] / (class_count[i] + clust_count_i[1])

      merge_count <- merge_count + 1
    }
  }

  max_merged_f1 <- apply(merging_f1, 1, max, na.rm = TRUE)
  merge_counts <- apply(merging_f1, 1, which.max) - 1

  list("f1" = max_merged_f1, "merges" = merge_counts)
}

#' Compute the precision for a single class with respect to a single cluster.
#'
#' @description
#' Compute a matrix of precision scores for each class-cluster pair, then
#' optionally obtain recall-weighted mean precision values for each class,
#' Finally, a weighted or unweighted mean may be computed across the classes.
#'
#' @inheritParams compute_max_f1
#'
#' @returns
#' Either a matrix, vector, or mean value, depending on `return_option`.
#'
#' @export
compute_precision <- function(
    class_labels, cluster_labels,
    return_option = c("vector", "value", "matrix"),
    weighted = TRUE) {
  return_option <- rlang::arg_match(
    return_option, c("value", "vector", "matrix")
  )

  f1_output <- compute_f1(class_labels, cluster_labels)
  class_count <- f1_output$class_count
  precision <- f1_output$precision
  recall <- f1_output$recall

  precision_vec <- rowSums(precision * recall)

  if (weighted) {
    precision_val <- stats::weighted.mean(precision_vec, class_count)
  } else {
    precision_val <- mean(precision_vec)
  }

  if (return_option == "matrix") {
    return(precision)
  } else if (return_option == "vector") {
    return(precision_vec)
  } else if (return_option == "value") {
    return(precision_val)
  } else {
    stop("Invalid return_option")
  }
}

#' Compute the entropy for a single class across all clusters.
#'
#' @description
#' Compute a matrix of entropy scores for each class-cluster pair, then
#' optionally obtain the summed total entropy values for each class,
#' Finally, a weighted or unweighted mean may be computed across the classes.
#'
#' @inheritParams compute_max_f1
#'
#' @returns
#' Either a matrix, vector, or mean value, depending on `return_option`.
#'
#' @export
compute_recall_entropy <- function(
    class_labels, cluster_labels,
    return_option = c("vector", "value", "matrix"),
    weighted = TRUE) {
  return_option <- rlang::arg_match(
    return_option, c("vector", "value", "matrix")
  )

  f1_output <- compute_f1(class_labels, cluster_labels)
  class_count <- f1_output$class_count
  recall <- f1_output$recall

  entropy_mat <- ifelse(recall == 0, 0, -recall * log(recall))
  entropy_vec <- rowSums(entropy_mat)

  if (weighted) {
    entropy_val <- stats::weighted.mean(entropy_vec, class_count)
  } else {
    entropy_val <- mean(entropy_vec)
  }

  if (return_option == "matrix") {
    return(entropy_mat)
  } else if (return_option == "vector") {
    return(entropy_vec)
  } else if (return_option == "value") {
    return(entropy_val)
  } else {
    stop("Invalid return_option")
  }
}
