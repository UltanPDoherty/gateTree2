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
#' @importFrom ggpubr ggarrange
#' @importFrom utils tail
#' @export
targeted_tree <- function(
  x,
  plusminus_table = expand.grid(rep(list(c(-1, 1)), nrow(x))),
  min_height = min_depth,
  min_depth = 1,
  min_val_cutoff = NULL,
  max_val_cutoff = NULL,
  use_boundaries = TRUE,
  show_plot = FALSE
) {
  min_size <- 100

  var_num <- ncol(x)
  obs_num <- nrow(x)
  path_num <- 1

  splits      <- scores <- matrix(NA, nrow = path_num, ncol = var_num)
  split_order <- signs  <- matrix(NA, nrow = path_num, ncol = var_num)
  splittable_vars <- matrix(NA, nrow = path_num, ncol = var_num)
  colnames(signs) <- colnames(plusminus_table)

  already_split <- matrix(FALSE, nrow = path_num, ncol = var_num)

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

  splittable_vars[1, ] <- !already_split[1, ] & common_variables[1, ]

  g <- 1
  k <- 1
  while (g <= path_num) {

    proposals <- propose_splits(
      x, subsetter[, g], splittable_vars[g, ],
      min_size, min_depth, min_height,
      use_boundaries
    )

    scenario <- proposals$scenario

    if (scenario %in% c("valley", "boundary")) {

      p_choice <- which.max(proposals$matrix[2, ])
      splits[g, p_choice] <- proposals$matrix[1, p_choice]
      scores[g, p_choice] <- proposals$matrix[2, p_choice]
      already_split[g, p_choice] <- TRUE

      for (j in which(pop_to_path == g)) {
        same_sign <- plusminus_table[, p_choice] == plusminus_table[j, p_choice]
        same_branch[, j] <- same_branch[, j] & same_sign
      }

      same_sign <- plusminus_table[, p_choice] == plusminus_table[k, p_choice]
      no_new_branch <- all(same_sign[pop_to_path == g])
      pop_to_path[pop_to_path == g & !same_sign] <- path_num + 1

      x_gp <- x[subsetter[, g], p_choice]
      more_gp <- x[, p_choice] > splits[g, p_choice]

      is_negative <- plusminus_table[g, p_choice] == -1
      refine_current <- if (is_negative) !more_gp else more_gp

      split_num[g] <- split_num[g] + 1
      split_order[g, p_choice] <- split_num[g]

      signs[g, p_choice] <- plusminus_table[k, p_choice]
      rownames(signs)[g] <-
        rownames(plusminus_table)[which(pop_to_path == g)[1]]

      plot_list[[g]][[split_num[g]]] <- plot_targeted_split(
        x_gp, g, p_choice, scores[g, p_choice],
        signs, scenario, splits[g, p_choice]
      )

      if (!no_new_branch) {
        split_num[path_num + 1] <- split_num[g]

        plot_list[[path_num + 1]] <- plot_list[[g]]

        start_node[path_num + 1] <- node_num

        splits <- rbind(splits, splits[g, ])

        scores <- rbind(scores, scores[g, ])

        already_split <- rbind(already_split, already_split[g, ])

        signs <- rbind(signs, signs[g, ])
        signs[path_num + 1, p_choice] <- - plusminus_table[k, p_choice]
        rownames(signs)[path_num + 1] <-
          rownames(plusminus_table)[which(pop_to_path == path_num + 1)[1]]

        split_order <- rbind(split_order, split_order[g, ])

        subsetter <- cbind(subsetter, subsetter[, g])
        subsetter[, path_num + 1] <- subsetter[, path_num + 1] & !refine_current

        splittable_vars <- rbind(splittable_vars, splittable_vars[g, ])

        plot_list[[path_num + 1]][[split_num[path_num + 1]]] <-
          plot_targeted_split(
            x_gp, path_num + 1, p_choice, scores[g, p_choice],
            signs, scenario, splits[g, p_choice]
          )

        parent_node[node_num + 2] <- utils::tail(path_nodes[[g]], 1)

        path_nodes[[path_num + 1]] <- append(path_nodes[[g]], node_num + 2)

        edge_name[node_num + 2] <- paste0(colnames(x)[p_choice],
                                          ifelse(is_negative, "+", "-"))

        node_name[node_num + 2] <- paste(edge_name[path_nodes[[path_num + 1]]],
                                         collapse = "/")

        is_leaf[node_num + 2] <- TRUE

        path_num <- path_num + 1
      }

      subsetter[, g] <- subsetter[, g] & refine_current

      parent_node[node_num + 1] <- utils::tail(path_nodes[[g]], 1)
      path_nodes[[g]] <- append(path_nodes[[g]], node_num + 1)
      edge_name[node_num + 1] <- paste0(colnames(x)[p_choice],
                                        ifelse(is_negative, "-", "+"))
      node_name[node_num + 1] <- paste(edge_name[path_nodes[[g]]],
                                       collapse = "/")

      is_leaf[parent_node[node_num + 1]] <- FALSE
      is_leaf[node_num + 1] <- TRUE

      node_num <- node_num + 1 + !no_new_branch

    } else {

      if (sum(subsetter[, g]) > 1) {
        missed_splits <- 0
        for (p in which(splittable_vars[g, ])) {
          missed_splits <- missed_splits + 1
          plot_list[[g]][[split_num[g] + missed_splits]] <- plot_targeted_split(
            x[subsetter[, g], p], g, p, depth = NA,
            signs, scenario, split_gp = NA
          )
        }

        scenario <- "explore"
        explore_valleys <- propose_valleys(
          x, subsetter[, g], !splittable_vars[g, ],
          5 * min_depth, 5 * min_height
        )
        explore_splits <- 0
        explore_check <- !is.na(explore_valleys[1, ])
        explore_check <- explore_check & sum(subsetter[, g]) >= min_size
        for (p in which(explore_check)) {
          explore_splits <- explore_splits + 1
          plot_list[[g]][[split_num[g] + missed_splits + explore_splits]] <-
            plot_targeted_split(
              x[subsetter[, g], p], g, p, depth = explore_valleys[2, p],
              signs, scenario, split_gp = explore_valleys[1, p]
            )
        }
      }

      g <- g + 1
      common_variables <- rbind(common_variables, NA)
    }

    if (any(pop_to_path == g)) {
      # if a path terminates without fully isolating a single population, then
      # the final number of paths and populations will be different because
      # there will be a path with multiple populations on it.
      # If population 1 and 2 are on path 1 when it terminates, then path 2 will
      # correspond to population 3, and so on.
      # Would it be easier to create duplicate paths in this case, so that g
      # and k are identical?
      #
      # Hang on, the second population does not necessarily lie on the second
      # path. E.g. if the first and second populations are distinguished by the
      # second split on path 1, path 2 starts at the first split.
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

      splittable_vars[g, ] <- !already_split[g, ] & common_variables[g, ]
    }
  }

  if (show_plot) {
    plot_paths(plot_list)
  }

  edge_df <- make_edge_df(parent_node, node_num, edge_name, node_name,
                          is_leaf, path_nodes)
  tree_plot <- make_tree_plot(edge_df)

  plot(tree_plot)

  labels <- c()
  unassigned <- rowSums(subsetter) == 0
  labels[!unassigned] <- apply(
    subsetter[!unassigned, , drop = FALSE], 1, which
  )
  labels[unassigned] <- 0

  signs[is.na(signs)] <- 0

  return(list(splits = splits,
              split_order = split_order,
              subsetter = subsetter,
              edge_df = edge_df,
              labels = labels,
              signs = signs,
              tree_plot = tree_plot))
}
