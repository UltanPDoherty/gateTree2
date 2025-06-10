compute_f1 <- function(class_labels, cluster_labels, i, j) {
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

  if (clust_num < class_num) {
    for (j in seq(clust_num + 1, class_num)) {
      f1[, j] <- rep(-Inf, class_num)
    }
  }

  if (!is.null(no_match_class)) {
    f1[which(rownames(f1) %in% no_match_class), ] <- rep(-Inf, clust_num)
  }
  if (!is.null(no_match_class)) {
    f1[, which(colnames(f1) %in% no_match_cluster)] <- rep(-Inf, class_num)
  }

  if (hungarian) {
    hungarian_order <- clue::solve_LSAP(f1, maximum = TRUE)

    f1 <- f1[, hungarian_order]

    f1_vec <- diag(f1)
  } else {
    f1_vec <- apply(f1, 1, max)
  }

  if (weighted) {
    f1_val <- stats::weighted.mean(f1_vec, class_count)
  } else {
    f1_val <- mean(f1_vec)
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

compute_merged_f1 <- function(
    class_labels, cluster_labels, merge_limit = NULL,
    no_match_class = NULL, no_match_cluster = NULL) {
  f1_output <- compute_f1(class_labels, cluster_labels)
  class_num <- f1_output$class_num
  clust_num <- f1_output$clust_num
  class_count <- f1_output$class_count
  clust_count <- f1_output$clust_count
  contingency <- f1_output$contingency
  precision <- f1_output$precision

  if (is.null(merge_limit)) {
    merge_limit <- clust_num
  } else {
    merge_limit <- max(clust_num, merge_limit)
  }

  merging_f1 <- matrix(nrow = class_num, ncol = clust_num)
  for (i in seq_len(class_num)) {
    precision_order <- order(precision[i, ], decreasing = TRUE)
    cont_temp_i <- contingency[i, precision_order]
    clust_count_temp <- clust_count[precision_order]

    cont_temp_i[1] <- cont_temp_i[1]
    merging_f1[i, 1] <-
      2 * cont_temp_i[1] / (class_count[i] + clust_count_temp[1])

    for (j in seq(2, merge_limit)) {
      cont_temp_i[j] <- cont_temp_i[j - 1] + cont_temp_i[j]
      clust_count_temp[j] <- clust_count_temp[j - 1] + clust_count_temp[j]
      merging_f1[i, j] <-
        2 * cont_temp_i[j] / (class_count[i] + clust_count_temp[j])
    }
  }
  max_merged_f1 <- apply(merging_f1, 1, max, na.rm = TRUE)
  merge_counts <- apply(merging_f1, 1, which.max) - 1

  list("f1" = max_merged_f1, "merges" = merge_counts)
}

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
