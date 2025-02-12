#' @title User-informed clustering decision tree.
#'
#' @description
#' Construct a semi-supervised decision tree to identify user-described
#' clusters.
#'
#' @param x Dataset in `matrix` or `data.frame` form.
#' @param plusminus_table
#' Table indicating whether each population (row) is positive (+1),
#' negative (-1), or neutral / unknown (0) for each variable (column).
#' @param order_table Table influencing the order of splits.
#' @param min_height
#' Minimum height, as a percentage of the height of the global density maximum,
#' for a peak to be recognised by [find_peaks].
#' @param min_depth
#' Minimum depth, as a percentage of the height of the global density maximum,
#' for a split to be returned by [find_valley].
#' @param min_scaled_bic_diff
#' Minimum value of difference between one-component and two-component BIC
#' divided by 2 log(obs_num).
#' @param min_size Minimum number of observations for a subset to be split.
#' @param min_val_cutoff
#' Minimum values per variable for observations to be included.
#' @param max_val_cutoff
#' Maximum values per variable for observations to be included.
#' @param temp_min_cutoff
#' Minimum values per variable for observations to be included in finding a
#' split.
#' @param temp_max_cutoff
#' Maximum values per variable for observations to be included in finding a
#' split.
#' @param use_boundaries Logical value.
#' @param show_plot
#' Logical vector of length 2. Indicating whether the split plots (`[1]`) and
#' the tree plot (`[2]`) should be plotted.
#' @param explore Logical value.
#'
#' @return A `list` object:
#' * `splits`: matrix of split locations for all variables and described
#' populations.
#' * `split_order`: matrix of split order for all variables and described
#' populations.
#' * `edge_df`: data.frame describing the tree's edges and nodes.
#' * `labels`: vector of integer cluster labels. `0` represents "Unassigned".
#' * `signs`: matrix of split signs for all variables and described populations.
#' * `tree_plot`: the clustering tree as a `ggplot` object.
#'
#' @export
gatetree <- function(
    x,
    plusminus_table = expand.grid(rep(list(c(-1, 1)), ncol(x))),
    order_table = array(0, dim = dim(plusminus_table)),
    min_height = min_depth,
    min_depth = 1,
    min_scaled_bic_diff = 0,
    min_size = 50,
    min_val_cutoff = NULL,
    max_val_cutoff = NULL,
    temp_min_cutoff = NULL,
    temp_max_cutoff = NULL,
    use_boundaries = TRUE,
    show_plot = c(FALSE, FALSE),
    explore = TRUE) {
  # Set the minimum requirement thresholds for the 'explore' stage.
  # These are higher than the user-defined thresholds for making splits.
  explore_min_height <- min(c(100, 5 * min_height))
  explore_min_depth <- min(c(100, 5 * min_depth))
  explore_min_size <- min_size

  var_num <- ncol(x)
  obs_num <- nrow(x)
  pop_num <- nrow(plusminus_table)
  path_num <- 1

  split_num <- c(0)
  node_num <- 1
  parent_node <- c(1)
  edge_name <- c("All")
  path_nodes <- list(c(1))
  node_name <- c("All")
  is_leaf <- c(TRUE)

  # ----------------------------------------------------------------------------

  # These 8 matrices all have the same dimensions.
  # splits will store the location / value of each split.
  # scores will store the depth % or BIC associated with each split.
  # split_order will store the order in which the splits occur.
  # splittable_vars keeps track of which variables can be split.
  # order_vars is used to implement `order_table`'s split restrictions.
  # common_variables: vars that have info for all of the path's populations.
  # already_split keeps track of which variables have been partitioned.
  splits <- scores <- split_order <- signs <-
    splittable_vars <- order_vars <- common_variables <- already_split <-
    matrix(NA, nrow = path_num, ncol = var_num)

  colnames(signs) <- colnames(plusminus_table)

  colnames(common_variables) <- colnames(x)
  common_variables[1, ] <- apply(plusminus_table == 0, 2, function(x) !any(x))

  already_split[] <- FALSE

  order_vars[1, ] <- split_num[1] >= order_table[1, ]

  splittable_vars[1, ] <-
    !already_split[1, ] & common_variables[1, ] & order_vars[1, ]

  # ----------------------------------------------------------------------------

  inside_cutoffs <- find_inside_cutoffs(x, min_val_cutoff, max_val_cutoff)
  inside_common <-
    apply(inside_cutoffs[, common_variables[1, ], drop = FALSE], 1, all)
  subsetter <- matrix(inside_common, nrow = obs_num, ncol = path_num)

  temp_in_cutoffs <- find_inside_cutoffs(x, temp_min_cutoff, temp_max_cutoff)
  temp_in_common <-
    apply(temp_in_cutoffs[, common_variables[1, ], drop = FALSE], 1, all)
  temp_subsetter <-
    matrix(inside_common & temp_in_common, nrow = obs_num, ncol = path_num)

  # ----------------------------------------------------------------------------

  pop_to_path <- rep(1, pop_num)

  g <- 1
  k <- 1
  plot_list <- list(list())
  while (g <= path_num) {
    proposals <- propose_splits(
      x, temp_subsetter[, g], splittable_vars[g, ],
      min_size, min_depth, min_height,
      min_scaled_bic_diff, use_boundaries
    )

    if (proposals$scenario %in% c("valley", "boundary")) {
      var_choice <- which.max(proposals$scores)

      splits[g, var_choice] <- proposals$splits[var_choice]
      scores[g, var_choice] <- proposals$scores[var_choice]
      already_split[g, var_choice] <- TRUE

      same_sign <-
        plusminus_table[, var_choice] == plusminus_table[k, var_choice]
      new_branch_created <- !all(same_sign[pop_to_path == g])
      pop_to_path[pop_to_path == g & !same_sign] <- path_num + 1

      x_gp <- x[subsetter[, g], var_choice]
      more_gp <- x[, var_choice] > splits[g, var_choice]

      is_negative <- plusminus_table[k, var_choice] == -1
      refine_current <- if (is_negative) !more_gp else more_gp

      split_num[g] <- split_num[g] + 1
      split_order[g, var_choice] <- split_num[g]

      signs[g, var_choice] <- plusminus_table[k, var_choice]
      rownames(signs)[g] <-
        rownames(plusminus_table)[which(pop_to_path == g)[1]]

      plot_list[[g]][[split_num[g]]] <- plot_gatetree_split(
        x_gp, g, var_choice, scores[g, var_choice],
        signs, proposals$scenario, splits[g, var_choice]
      )

      if (new_branch_created) {
        split_num[path_num + 1] <- split_num[g]

        plot_list[[path_num + 1]] <- plot_list[[g]]

        splits <- rbind(splits, splits[g, ])

        scores <- rbind(scores, scores[g, ])

        already_split <- rbind(already_split, already_split[g, ])

        signs <- rbind(signs, signs[g, ])
        signs[path_num + 1, var_choice] <- -plusminus_table[k, var_choice]
        rownames(signs)[path_num + 1] <-
          rownames(plusminus_table)[which(pop_to_path == path_num + 1)[1]]

        split_order <- rbind(split_order, split_order[g, ])

        subsetter <- cbind(subsetter, subsetter[, g])
        subsetter[, path_num + 1] <- subsetter[, path_num + 1] & !refine_current

        splittable_vars <- rbind(splittable_vars, splittable_vars[g, ])
        order_vars <- rbind(order_vars, order_vars[g, ])

        plot_list[[path_num + 1]][[split_num[path_num + 1]]] <-
          plot_gatetree_split(
            x_gp, path_num + 1, var_choice, scores[g, var_choice],
            signs, proposals$scenario, splits[g, var_choice]
          )

        parent_node[node_num + 2] <- utils::tail(path_nodes[[g]], 1)

        path_nodes[[path_num + 1]] <- append(path_nodes[[g]], node_num + 2)

        edge_name[node_num + 2] <-
          paste0(colnames(x)[var_choice], ifelse(is_negative, "+", "-"))

        node_name[node_num + 2] <-
          paste(edge_name[path_nodes[[path_num + 1]]], collapse = "\n")

        is_leaf[node_num + 2] <- TRUE

        path_num <- path_num + 1
      }

      subsetter[, g] <- subsetter[, g] & refine_current

      parent_node[node_num + 1] <- utils::tail(path_nodes[[g]], 1)
      path_nodes[[g]] <- append(path_nodes[[g]], node_num + 1)
      edge_name[node_num + 1] <-
        paste0(colnames(x)[var_choice], ifelse(is_negative, "-", "+"))
      node_name[node_num + 1] <-
        paste(edge_name[path_nodes[[g]]], collapse = "\n")

      is_leaf[parent_node[node_num + 1]] <- FALSE
      is_leaf[node_num + 1] <- TRUE

      node_num <- node_num + 1 + new_branch_created
    } else {
      if (sum(subsetter[, g]) > min_size) {
        missed_splits <- 0
        for (p in which(splittable_vars[g, ])) {
          missed_splits <- missed_splits + 1
          plot_list[[g]][[split_num[g] + missed_splits]] <-
            plot_gatetree_split(
              x[subsetter[, g], p], g, p,
              score = NA,
              signs, proposals$scenario, split_gp = NA
            )
        }

        plot_list <- explore_plots(
          explore,
          x, g, subsetter, splittable_vars,
          explore_min_depth, explore_min_height, explore_min_size,
          plot_list, split_num, missed_splits, signs
        )
      }

      g <- g + 1
      common_variables <- rbind(common_variables, NA)
    }

    if (any(pop_to_path == g)) {
      k <- match(g, pop_to_path)

      common_variables[g, ] <- apply(
        plusminus_table[pop_to_path == g, , drop = FALSE],
        2,
        function(x) all(x != 0)
      )

      inside_common <-
        apply(inside_cutoffs[, common_variables[1, ], drop = FALSE], 1, all)
      subsetter[, g] <- subsetter[, g] & inside_common

      temp_subsetter <- subsetter
      temp_in_common <-
        apply(temp_in_cutoffs[, common_variables[1, ], drop = FALSE], 1, all)
      temp_subsetter[, g] <- subsetter[, g] & temp_in_common

      splittable_vars[g, ] <- !already_split[g, ] & common_variables[g, ]

      order_vars[g, ] <- split_num[g] >= order_table[k, ]
      splittable_vars[g, ] <- splittable_vars[g, ] & order_vars[g, ]
    }
  }

  plot_paths(plot_list, show_plot = show_plot[1])

  edge_df <-
    make_edge_df(
      parent_node, node_num, edge_name, node_name, is_leaf, path_nodes
    )
  tree_plot <- make_tree_plot(edge_df, show_plot[2])

  labels <- c()
  unassigned <- rowSums(subsetter) == 0
  labels[!unassigned] <-
    apply(subsetter[!unassigned, , drop = FALSE], 1, which)
  labels[unassigned] <- 0

  signs[is.na(signs)] <- 0

  return(list(
    splits = splits,
    split_order = split_order,
    subsetter = subsetter,
    edge_df = edge_df,
    labels = labels,
    signs = signs,
    tree_plot = tree_plot
  ))
}
