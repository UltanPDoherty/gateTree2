#' Partition data set based on unsupervised marginal splits.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param min_height Minimum height, as a percentage of the height of the global
#'                   density maximum, for a peak to be recognised by find_peaks.
#' @param min_depth Minimum depth, as a percentage of the height of the global
#'                   density maximum, for a split to be returned by find_valley.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#' @param show_plot Logical value.
#'
#' @return List: splits, split_order, subsetter, edge_df, labels, signs.
#' @import ggplot2
#' @import ggraph
#' @importFrom ggpubr ggarrange
#' @importFrom utils tail
#' @export
exploratory_tree <- function(
  x,
  min_height = 0.5,
  min_depth = 0.5,
  min_val_cutoff = NULL,
  max_val_cutoff = NULL,
  show_plot = FALSE
) {
  var_num <- ncol(x)
  obs_num <- nrow(x)
  path_num <- 1

  splits      <- scores <- matrix(NA, nrow = path_num, ncol = var_num)
  split_order <- signs  <- matrix(NA, nrow = path_num, ncol = var_num)

  inside_cutoffs <- find_inside_cutoffs(x, min_val_cutoff, max_val_cutoff)
  inside_all <- apply(inside_cutoffs, 1, all)

  subsetter <- matrix(inside_all, nrow = obs_num, ncol = path_num)

  progress <- matrix(FALSE, nrow = path_num, ncol = var_num)

  plot_list <- list(list())

  split_num <- c(0)
  node_num <- 1
  start_node <- c(1)
  current_node <- c(1)
  parent_node <- c(1)
  edge_name <- c("All")
  path_nodes <- list(c(1))
  node_name <- c("All")
  is_leaf <- c(TRUE)

  g <- 1
  while (g <= path_num) {
    if (sum(subsetter[, g]) < 100) {
      found_valley <- FALSE
    } else {
      valleys <- propose_valleys(x, g, var_num, subsetter, progress,
                                 min_depth, min_height)
      found_valley <- any(!is.na(valleys[1, ]))
    }

    if (found_valley) {

      p_choice <- which.max(valleys[2, ])
      splits[g, p_choice] <- valleys[1, p_choice]
      scores[g, p_choice] <- valleys[2, p_choice]
      progress[g, p_choice] <- TRUE

      x_gp <- x[subsetter[, g], p_choice]
      less_gp <- x[, p_choice] < splits[g, p_choice]

      split_num[g] <- split_num[path_num + 1] <- split_num[g] + 1

      plot_list[[path_num + 1]] <- plot_list[[g]]
      start_node[path_num + 1] <- node_num

      splits <- rbind(splits, splits[g, ])
      scores <- rbind(scores, scores[g, ])
      progress <- rbind(progress, progress[g, ])

      signs <- rbind(signs, signs[g, ])
      signs[g, p_choice] <- -1
      signs[path_num + 1, p_choice] <- +1

      split_order <- rbind(split_order, split_order[g, ])

      subsetter <- cbind(subsetter, subsetter[, g])
      subsetter[, g] <- subsetter[, g] & less_gp
      subsetter[, path_num + 1] <- subsetter[, path_num + 1] & !less_gp

      plot_list[[g]][[split_num[g]]] <- exploratory_plot(
        x_gp, splits, g, p_choice, subsetter, scores,
        is_negative = TRUE, var_name = colnames(x)[p_choice]
      )
      plot_list[[path_num + 1]][[split_num[path_num + 1]]] <- exploratory_plot(
        x_gp, splits, path_num + 1, p_choice, subsetter, scores,
        is_negative = FALSE, var_name = colnames(x)[p_choice]
      )

      parent_node[node_num + 1] <- utils::tail(path_nodes[[g]], 1)
      parent_node[node_num + 2] <- utils::tail(path_nodes[[g]], 1)

      path_nodes[[path_num + 1]] <- path_nodes[[g]]
      path_nodes[[g]] <- append(path_nodes[[g]],
                                node_num + 1)
      path_nodes[[path_num + 1]] <- append(path_nodes[[path_num + 1]],
                                           node_num + 2)

      edge_name[node_num + 1] <- paste0(colnames(x)[p_choice], "-")
      edge_name[node_num + 2] <- paste0(colnames(x)[p_choice], "+")

      node_name[node_num + 1] <- paste(edge_name[path_nodes[[g]]],
                                       collapse = "/")
      node_name[node_num + 2] <- paste(edge_name[path_nodes[[path_num + 1]]],
                                       collapse = "/")

      is_leaf[parent_node[node_num + 1]] <- FALSE
      is_leaf[node_num + 1] <- TRUE
      is_leaf[node_num + 2] <- TRUE

      path_num <- path_num + 1

      node_num <- node_num + 2
    } else {
      g <- g + 1
      current_node[g] <- start_node[g]
    }
  }


  if (show_plot) {
    plot_paths(plot_list)
  }

  edge_df <- make_edge_df(parent_node, node_num, edge_name, node_name,
                          is_leaf, path_nodes)
  tree_plot <- make_tree_plot(edge_df)

  plot(tree_plot$graph)

  labels <- c()
  labels[inside_all] <- apply(subsetter[inside_all, , drop = FALSE], 1, which)
  labels[!inside_all] <- 0

  signs[is.na(signs)] <- 0

  return(list(splits = splits, split_order = split_order, subsetter = subsetter,
              edge_df = edge_df, labels = labels, signs = signs))
}
