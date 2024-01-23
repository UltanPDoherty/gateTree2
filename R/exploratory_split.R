#' Partition data set based on unsupervised marginal splits.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param min_height Minimum height for a peak to be recognised by find_peaks.
#' @param min_depth Minimum depth for a split to be returned by find_valley.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#' @param show_plot Logical value.
#'
#' @return splits, typemarker, subsetter
#' @import ggplot2
#' @import ggraph
#' @importFrom ggpubr ggarrange
#' @importFrom graphics par
#' @importFrom utils tail
#' @export
exploratory_split <- function(
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
  start_layer <- c(1)
  node_number <- 1
  start_node <- c(1)
  current_node <- c(1)
  parent_node <- c(1)
  edge_name <- c("All")
  path_nodes <- list(c(1))
  node_name <- c("All")
  is_leaf <- c(TRUE)

  g <- 1
  par(mfrow = c(2, 2))
  while (g <= path_num) {
    if (sum(subsetter[, g]) < 100) {
      found_valley <- FALSE
    } else {
      valleys <- propose_valleys(x, g, var_num, subsetter, progress,
                                 min_depth, min_height)
      found_valley <- any(!is.na(valleys[1, ]))
    }

    if (found_valley) {
      split_num[g] <- split_num[g] + 1

      splits <- rbind(splits, splits[g, ])
      scores <- rbind(scores, scores[g, ])
      signs <- rbind(signs, signs[g, ])
      split_order <- rbind(split_order, split_order[g, ])
      progress <- rbind(progress, progress[g, ])
      subsetter <- cbind(subsetter, subsetter[, g])
      plot_list[[path_num + 1]] <- plot_list[[g]]
      split_num[path_num + 1] <- split_num[g]
      start_layer[path_num + 1] <- split_num[g]
      start_node[path_num + 1] <- node_number

      p_choice <- which.max(valleys[2, ])

      splits[g, p_choice] <- valleys[1, p_choice]
      scores[g, p_choice] <- valleys[2, p_choice]
      signs[g, p_choice] <- -1
      split_order[g, p_choice] <- split_num[g]
      progress[g, p_choice] <- TRUE
      x_gp <- x[subsetter[, g], p_choice]
      less_gp <- x[, p_choice] < splits[g, p_choice]
      subsetter[, g] <- subsetter[, g] & less_gp

      splits[path_num + 1, p_choice] <- valleys[1, p_choice]
      scores[path_num + 1, p_choice] <- valleys[2, p_choice]
      signs[path_num + 1, p_choice] <- +1
      split_order[path_num + 1, p_choice] <- split_num[path_num + 1]
      progress[path_num + 1, p_choice] <- TRUE
      subsetter[, path_num + 1] <- subsetter[, path_num + 1] & !less_gp

      plot_list[[g]][[split_num[g]]] <- exploratory_plot(
        x_gp, splits, g, p_choice, subsetter, scores,
        is_negative = TRUE, var_name = colnames(x)[p_choice]
      )
      plot_list[[path_num + 1]][[split_num[path_num + 1]]] <- exploratory_plot(
        x_gp, splits, path_num + 1, p_choice, subsetter, scores,
        is_negative = FALSE, var_name = colnames(x)[p_choice]
      )

      parent_node[node_number + 1] <- tail(path_nodes[[g]], 1)
      parent_node[node_number + 2] <- tail(path_nodes[[g]], 1)


      path_nodes[[path_num + 1]] <- path_nodes[[g]]
      path_nodes[[g]] <- append(path_nodes[[g]],
                                node_number + 1)
      path_nodes[[path_num + 1]] <- append(path_nodes[[path_num + 1]],
                                           node_number + 2)

      edge_name[node_number + 1] <- paste0(colnames(x)[p_choice], "-")
      edge_name[node_number + 2] <- paste0(colnames(x)[p_choice], "+")

      node_name[node_number + 1] <- paste(edge_name[path_nodes[[g]]],
                                          collapse = "/")
      node_name[node_number + 2] <- paste(edge_name[path_nodes[[path_num + 1]]],
                                          collapse = "/")

      is_leaf[parent_node[node_number + 1]] <- FALSE
      is_leaf[node_number + 1] <- TRUE
      is_leaf[node_number + 2] <- TRUE

      path_num <- path_num + 1

      node_number <- node_number + 2
    } else {
      g <- g + 1
      current_node[g] <- start_node[g]
    }
  }

  leaf_name <- rep("", node_number)
  leaf_name[is_leaf] <- node_name[is_leaf]

  edge_df <- data.frame(
    parent = parent_node,
    node = 1:node_number,
    edge_name = edge_name,
    node_name = node_name,
    is_leaf = is_leaf,
    leaf_name = leaf_name
  )

  path_order <- unique(Reduce(c, path_nodes))

  ggraph_order <- rep(NA, node_number)
  ggraph_order[is_leaf] <- (sum(!is_leaf) + 1):node_number
  ggraph_order[!is_leaf] <- order(path_order[!is_leaf[path_order]])

  node_name <- gsub("All/", "", node_name)

  my_ggraph <- ggraph::ggraph(edge_df, layout = "tree") +
    ggraph::geom_edge_elbow(aes(label = edge_name),
                            angle_calc = "across",
                            label_dodge = grid::unit(-0.03, "npc")) +
    ggraph::geom_node_label(aes(label = leaf_name[order(ggraph_order)])) +
    ggraph::theme_graph()

  arranged <- list()
  for (g in 1:path_num) {

    arranged[[g]] <- ggpubr::ggarrange(plotlist = plot_list[[g]],
                                       ncol = 2, nrow = 2)
    if (show_plot && is.ggplot(arranged[[g]])) {
      plot(arranged[[g]])
    } else if (show_plot) {
      for (j in seq_along(arranged[[g]])) {
        plot(arranged[[g]][[j]])
      }
    }
  }

  plot(my_ggraph)

  labels <- c()
  labels[inside_all] <- apply(subsetter[inside_all, ], 1, which)
  labels[!inside_all] <- 0

  signs[is.na(signs)] <- 0

  return(list(splits = splits, split_order = split_order, subsetter = subsetter,
              edge_df = edge_df, labels = labels, signs = signs))
}

#===============================================================================

#' Plotting function for exploratory_split.
#'
#' @param x_gp Univariate data.
#' @param splits split matrix.
#' @param g The current pathway.
#' @param p_choice The chosen variable.
#' @param subsetter The subsetting matrix.
#' @param scores The scores matrix.
#' @param is_negative Whether to display the negative or positive side.
#' @param var_name The name of the variable.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats density
#' @export
exploratory_plot <- function(
  x_gp,
  splits,
  g,
  p_choice,
  subsetter,
  scores,
  is_negative,
  var_name
) {
  scale01_gp <- scale01(x_gp)
  dens_gp <- stats::density(scale01_gp$y)

  dens_gp$y <- dens_gp$y / max(dens_gp$y) * 100

  trans_split_gp <- scale01(splits[g, p_choice],
                            scale01_gp$min, scale01_gp$max)$y

  xleft <- ifelse(is_negative, 0, trans_split_gp)
  xright <- ifelse(is_negative, trans_split_gp, 1)

  rect_col <- ifelse(is_negative, "lightgreen", "lightcoral")

  size_before <- length(x_gp)
  size_after <- sum(subsetter[, g])

  dens_gp_x <- dens_gp$x
  dens_gp_y <- dens_gp$y
  dens_gp_df <- data.frame(dens_gp_x, dens_gp_y)

  ggplot(dens_gp_df, aes(x = dens_gp_x, y = dens_gp_y)) +
    geom_rect(aes(xmin = xleft, xmax = xright, ymin = 0, ymax = max(dens_gp$y)),
              fill = rect_col) +
    geom_line() +
    geom_vline(xintercept = trans_split_gp) +
    labs(x = paste0("N before = ", size_before, ", ", "N after = ", size_after),
         y = "Density %",
         title = paste0("g = ", g, ", p = ", p_choice, ", ",
                        "depth = ", round(scores[g, p_choice], 1), "%"),
         subtitle = var_name) +
    theme_bw()
}

#===============================================================================

#' Plotting function for exploratory_split.
#'
#' @param signs The signs matrix.
#' @param typemarker The cell type - marker matrix.
#' @param subsetter The subsetting matrix.
#'
#' @return A ggplot object.
#' @export
merge_subsets <- function(signs, typemarker, subsetter) {
  match_array <- array(dim = c(nrow(signs), nrow(typemarker), ncol(signs)))
  match_matrix <- matrix(nrow = nrow(signs), ncol = nrow(typemarker))
  for (i in seq_len(nrow(signs))) {
    for (j in seq_len(nrow(typemarker))) {
      for (p in seq_len(ncol(signs))) {
        match_array[i, j, p] <- ifelse(typemarker[j, p] != 0,
                                       signs[i, p] == typemarker[j, p],
                                       TRUE)
      }
      match_matrix[i, j] <- all(match_array[i, j, ])
    }
  }

  merged_subsetter <- matrix(nrow = nrow(subsetter), ncol = ncol(typemarker))
  match_count <- c()
  for (j in seq_len(nrow(typemarker))) {
    match_count[j] <- sum(match_matrix[, j])
    if (match_count[j] == 0) {
      merged_subsetter[, j] <- FALSE
    } else if (match_count[j] == 1) {
      merged_subsetter[, j] <- subsetter[, match_matrix[, j]]
    } else {
      merged_subsetter[, j] <- apply(subsetter[, match_matrix[, j]], 1, any)
    }
  }

  return(list(match_matrix = match_matrix,
              merged_subsetter = merged_subsetter))
}
