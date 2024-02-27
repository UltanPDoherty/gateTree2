#' Partition data set based on unsupervised marginal splits.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param plusminus_table Table indicated whether each group (row) is positive
#'                        (+1), negative (-1), or neutral / unknown (0) for each
#'                        variable (column).
#' @param min_height Minimum height, as a percentage of the height of the global
#'                   density maximum, for a peak to be recognised by find_peaks.
#' @param min_depth Minimum depth, as a percentage of the height of the global
#'                   density maximum, for a split to be returned by find_valley.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#' @param use_boundaries Logical value.
#' @param show_plot Logical value.
#'
#' @return List: splits, split_order, subsetter, edge_df, labels, signs.
#' @import ggplot2
#' @import ggraph
#' @importFrom ggpubr ggarrange
#' @importFrom utils tail
#' @export
targeted_tree <- function(
    x,
    plusminus_table,
    min_height = min_depth,
    min_depth = 1,
    min_val_cutoff = NULL,
    max_val_cutoff = NULL,
    use_boundaries = TRUE,
    show_plot = FALSE
) {
  var_num <- ncol(x)
  obs_num <- nrow(x)
  path_num <- 1

  splits      <- scores <- matrix(NA, nrow = path_num, ncol = var_num)
  split_order <- signs  <- matrix(NA, nrow = path_num, ncol = var_num)

  already_split <- matrix(FALSE, nrow = path_num, ncol = var_num)

  plot_list <- list(list())

  split_num <- c(0)
  start_layer <- c(1)
  node_num <- 1
  start_node <- c(1)
  current_node <- c(1)
  parent_node <- c(1)
  edge_name <- c("All")
  path_nodes <- list(c(1))
  node_name <- c("All")
  is_leaf <- c(TRUE)

  pop_num <- nrow(plusminus_table)
  same_branch <- matrix(TRUE, nrow = pop_num, ncol = pop_num)

  common_variables <- matrix(nrow = path_num, ncol = var_num)
  colnames(common_variables) <- colnames(x)
  common_variables[1, ] <- apply(plusminus_table == 0, 2, function(x) !any(x))

  inside_cutoffs <- find_inside_cutoffs(x, min_val_cutoff, max_val_cutoff)
  inside_common <- apply(
    inside_cutoffs[, common_variables[1, ], drop = FALSE], 1, all
  )
  subsetter <- matrix(inside_common, nrow = obs_num, ncol = path_num)

  pop_to_path <- rep(1, pop_num)
  path_num <- 1

  g <- 1
  k <- 1
  while (g <= path_num) {
    if (sum(subsetter[, g]) < 100) {
      found_valley <- FALSE
      found_boundary <- FALSE
    } else {
      proposals <- propose_valleys(x, g, var_num, subsetter, already_split,
                                   min_depth, min_height, common_variables)
      found_valley <- any(!is.na(proposals[1, ]))
      found_boundary <- FALSE

      if (!found_valley && use_boundaries) {
        proposals <- propose_boundaries(x, g, var_num, subsetter, already_split,
                                        common_variables)
        found_boundary <- any(!is.na(proposals[1, ]))
        if (found_boundary) {
          scenario <- "boundary"
        }
      } else if (found_valley) {
        found_boundary <- NA
        scenario <- "valley"
      }
    }

    if (found_valley || found_boundary) {

      p_choice <- which.max(proposals[2, ])
      splits[g, p_choice] <- proposals[1, p_choice]
      scores[g, p_choice] <- proposals[2, p_choice]
      already_split[g, p_choice] <- TRUE

      for (j in which(pop_to_path == g)) {
        same_sign <- plusminus_table[, p_choice] == plusminus_table[j, p_choice]
        same_branch[, j] <- same_branch[, j] & same_sign
      }

      same_sign <- plusminus_table[, p_choice] == plusminus_table[k, p_choice]
      no_new_branch <- all(same_sign[pop_to_path == g])
      pop_to_path[pop_to_path == g & !same_sign] <- path_num + 1

      x_gp <- x[subsetter[, g], p_choice]
      less_gp <- x[, p_choice] < splits[g, p_choice]

      is_negative <- plusminus_table[g, p_choice] == -1
      if (is_negative) {
        refine_current <- less_gp
      } else {
        refine_current <- !less_gp
      }

      split_num[g] <- split_num[g] + 1
      split_order[g, p_choice] <- split_num[g]

      if (!no_new_branch) {
        split_num[path_num + 1] <- split_num[g]
        plot_list[[path_num + 1]] <- plot_list[[g]]
        start_node[path_num + 1] <- node_num
        splits <- rbind(splits, splits[g, ])
        scores <- rbind(scores, scores[g, ])
        already_split <- rbind(already_split, already_split[g, ])
        signs <- rbind(signs, signs[g, ])
        signs[path_num + 1, p_choice] <- - plusminus_table[k, p_choice]
        split_order <- rbind(split_order, split_order[g, ])
        subsetter <- cbind(subsetter, subsetter[, g])
        subsetter[, path_num + 1] <- subsetter[, path_num + 1] & !refine_current
      }

      signs[g, p_choice] <- plusminus_table[k, p_choice]
      colnames(signs) <- colnames(plusminus_table)

      subsetter[, g] <- subsetter[, g] & refine_current

      plot_list[[g]][[split_num[g]]] <- plot_targeted_split(
        x_gp, g, p_choice, scores[g, p_choice],
        signs, scenario, splits[g, p_choice]
      )

      if (!no_new_branch) {
        plot_list[[path_num + 1]][[split_num[path_num + 1]]] <- plot_targeted_split(
          x_gp, path_num + 1, p_choice, scores[g, p_choice],
          signs, scenario, splits[g, p_choice]
        )
      }

      if (!no_new_branch) {
        parent_node[node_num + 2] <- utils::tail(path_nodes[[g]], 1)
        path_nodes[[path_num + 1]] <- append(path_nodes[[g]], node_num + 2)
        edge_name[node_num + 2] <- paste0(colnames(x)[p_choice],
                                          ifelse(is_negative, "+", "-"))
        node_name[node_num + 2] <- paste(edge_name[path_nodes[[path_num + 1]]],
                                         collapse = "/")
        is_leaf[node_num + 2] <- TRUE
        path_num <- path_num + 1
      }

      parent_node[node_num + 1] <- utils::tail(path_nodes[[g]], 1)
      path_nodes[[g]] <- append(path_nodes[[g]], node_num + 1)
      edge_name[node_num + 1] <- paste0(colnames(x)[p_choice],
                                        ifelse(is_negative, "-", "+"))
      node_name[node_num + 1] <- paste(edge_name[path_nodes[[g]]],
                                       collapse = "/")

      is_leaf[parent_node[node_num + 1]] <- FALSE
      is_leaf[node_num + 1] <- TRUE

      node_num <- node_num + 1 + !no_new_branch

      if (any(pop_to_path == g)) {
        k <- match(g, pop_to_path)
        current_node[g] <- start_node[g]
        common_variables[g, ] <- apply(
          plusminus_table[pop_to_path == g, , drop = FALSE],
          2,
          function(x) all(x != 0)
        )
        inside_common <- apply(
          inside_cutoffs[, common_variables[1, ], drop = FALSE], 1, all
        )
        subsetter[, g] <- subsetter[, g] & inside_common
      }

    } else {
      for (p in which(!(plusminus_table[g, ] == 0) & !already_split[g, ])) {
        plot_list[[g]][[split_num[g] + 1]] <- plot_targeted_split(
          x[subsetter[, g], p], g, p, depth = NA,
          signs, scenario = "nothing", split_gp = NA
        )
      }

      g <- g + 1

      if (any(pop_to_path == g)) {
        k <- match(g, pop_to_path)
        current_node[g] <- start_node[g]
        common_variables <- rbind(
          common_variables,
          apply(plusminus_table[pop_to_path == g, , drop = FALSE], 2,
                function(x) all(x != 0))
        )
        inside_common <- apply(
          inside_cutoffs[, common_variables[1, ], drop = FALSE], 1, all
        )
        subsetter[, g] <- subsetter[, g] & inside_common
      }
    }
  }

  leaf_name <- rep("", node_num)
  leaf_name[is_leaf] <- node_name[is_leaf]

  if (show_plot) {
    plot_paths(plot_list)
  }

  tree_plot <- make_tree_plot(parent_node, node_num, edge_name, node_name,
                              is_leaf, leaf_name, path_nodes)

  plot(tree_plot$graph)

  labels <- c()
  unassigned <- rowSums(subsetter) == 0
  # labels[inside_all & !unassigned] <- apply(
  #   subsetter[inside_all & !unassigned, , drop = FALSE], 1, which
  # )
  # labels[!inside_all] <- 0
  labels[!unassigned] <- apply(
    subsetter[!unassigned, , drop = FALSE], 1, which
  )
  labels[unassigned] <- 0

  signs[is.na(signs)] <- 0

  return(list(splits = splits,
              split_order = split_order,
              subsetter = subsetter,
              edge_df = tree_plot$df,
              labels = labels,
              signs = signs))
}



