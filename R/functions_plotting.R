#===============================================================================

#' @title Density plot with a `gateTree` split.
#'
#' @description
#' Plot a univariate kernel density estimate of the data with the `gateTree`
#' split illustrated, and information about the split in the title and subtitle.
#'
#'
#' @inheritParams gatetree
#' @param x_gp The data to be displayed, should be only pathway g & variable p.
#' @param g The pathway number.
#' @param p The variable number.
#' @param score The valley depth percentage or the scaled BIC difference.
#' @param scenario `"valley"`, `"boundary"`, `"nothing"`, or `"explore"`.
#' @param split_gp The split value.
#'
#' @return `ggplot` object.
plot_gatetree_split <- function(x_gp, g, p, score, plusminus_table,
                                scenario, split_gp) {

  # colours from ggokabeito package
  rect_col <- switch(
    scenario,
    "valley" = "#F0E442",
    "boundary" = "#56B4E9",
    "nothing" = NA,
    "explore" = "#CC79A7"
  )
  score <- switch(
    scenario,
    "valley"   = score,
    "boundary" = score,
    "nothing"  = NA,
    "explore"  = score
  )
  score_title <- switch(
    scenario,
    "valley"   = paste0("depth = ", round(score, 1), "%"),
    "boundary" = paste0("scaled_BIC_diff = ", round(score, 1)),
    "nothing"  = NA,
    "explore"  = paste0("depth = ", round(score, 1), "%")
  )
  line_type <- switch(
    scenario,
    "valley"   = "solid",
    "boundary" = "dashed",
    "nothing"  = "blank",
    "explore"  = "dotted"
  )
  is_negative <- switch(
    scenario,
    "valley"   = (plusminus_table[g, p] == -1),
    "boundary" = (plusminus_table[g, p] == -1),
    "nothing"  = NA,
    "explore"  = FALSE
  )

  scale01_gp <- scale01(x_gp)
  dens_gp <- stats::density(scale01_gp$y)

  trans_split_gp <- scale01(split_gp,
                            scale01_gp$min, scale01_gp$max)$y

  xleft <- ifelse(is_negative, 0, trans_split_gp)
  xright <- ifelse(is_negative, trans_split_gp, 1)

  if (scenario == "nothing") {
    xleft <- 0
    xright <- 1
  }

  size_before <- length(x_gp)
  if (is.na(is_negative)) {
    size_after <- NA
  } else if (is_negative) {
    size_after <- sum(x_gp < split_gp)
  } else {
    size_after <- sum(x_gp > split_gp)
  }

  dens_gp_x <- dens_gp$x
  dens_gp_y <- dens_gp$y / max(dens_gp$y) * 100
  dens_gp_df <- data.frame(dens_gp_x, dens_gp_y)

  gg <- ggplot(dens_gp_df, aes(x = dens_gp_x, y = dens_gp_y)) +
    geom_rect(aes(xmin = xleft, xmax = xright, ymin = 0, ymax = max(dens_gp_y)),
              fill = rect_col, na.rm = TRUE) +
    geom_line() +
    geom_vline(xintercept = trans_split_gp, linetype = line_type,
               na.rm = TRUE) +
    labs(
      title = paste0("g = ", g, ", p = ", p, ", ", score_title),
      subtitle = paste0("Path: ", rownames(plusminus_table)[g], ", ",
                        "Var: ", colnames(plusminus_table)[p], ", ",
                        "\"", scenario, "\""),
      x = paste0("N before = ", size_before, ", N after = ", size_after),
      y = "Density %"
    ) +
    theme_bw()

  return(gg)
}

#===============================================================================

#' @title Compile edge and node information.
#'
#' @description
#' Compile edge and node information into a `data.frame` for [make_tree_plot].
#'
#'
#' @param parent_node Integer vector: The number of each node's parent node.
#' @param node_number Integer: The total number of nodes.
#' @param edge_name Character vector: The name of each edge.
#' @param node_name Character vector: The name of each node.
#' @param is_leaf Logical vector: whether each node is a leaf / terminal node.
#' @param path_nodes List of nodes on each path.
#'
#' @return `data.frame` with columns:
#' * parent
#' * node
#' * edge_name
#' * node_name
#' * is_leaf
#' * ggraph_order
make_edge_df <- function(parent_node, node_number, edge_name, node_name,
                         is_leaf, path_nodes) {

  leaf_name <- rep(NA, node_number)
  leaf_name[is_leaf] <- node_name[is_leaf]
  leaf_name <- gsub("All\n", "", leaf_name)

  path_order <- unique(Reduce(c, path_nodes))

  ggraph_order <- rep(NA, node_number)
  ggraph_order[is_leaf] <- (sum(!is_leaf) + 1):node_number
  ggraph_order[!is_leaf] <- order(path_order[!is_leaf[path_order]])

  edge_df <- data.frame(
    parent = parent_node,
    node = 1:node_number,
    edge_name = edge_name,
    node_name = gsub("All/", "", node_name),
    is_leaf = is_leaf,
    leaf_name = leaf_name,
    ggraph_order = ggraph_order
  )

  return(edge_df)
}

#===============================================================================

#' @title Construct tree diagram for `gateTree`.
#'
#' @description
#' Construct a tree diagram from the information contained in `edge_df`.
#'
#' @param edge_df Output from [make_edge_df].
#' @param show_plot Logical: should the tree diagram be plotted?
#'
#' @import ggraph
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#'
#' @return `ggraph` object.
make_tree_plot <- function(edge_df, show_plot = FALSE) {

  tree_graph <- igraph::graph_from_data_frame(d = edge_df[, 1:3],
                                              v = edge_df[, c(2, 4:7)])

  tree_tbl_graph <- tidygraph::as_tbl_graph(tree_graph)

  leaf_name <- NULL # to silence "no visible binding" check
  tree_plot <- tree_tbl_graph |>
    ggraph::ggraph(layout = "tree") +
    ggraph::geom_edge_elbow2(
      show.legend = FALSE,
      angle_calc = "across",
      label_push = grid::unit(0.075, "npc"),
      label_dodge = grid::unit(0.075, "npc")
    ) +
    ggraph::geom_node_label(
      aes(label = leaf_name,
          colour = leaf_name),
      show.legend = FALSE
    ) +
    ggraph::theme_graph() +
    scale_colour_manual(values = rep("black", nrow(edge_df)))

  if (show_plot) {
    plot(tree_plot)
  }

  return(tree_plot)
}

#===============================================================================

plot_paths <- function(plot_list, show_plot) {
  if (show_plot) {
    arranged <- list()
    for (g in seq_along(plot_list)) {

      arranged[[g]] <- ggpubr::ggarrange(plotlist = plot_list[[g]],
                                         ncol = 2, nrow = 2)
      if (is.ggplot(arranged[[g]])) {
        plot(arranged[[g]])
      } else {
        for (j in seq_along(arranged[[g]])) {
          plot(arranged[[g]][[j]])
        }
      }
    }
  }
}

#===============================================================================

explore_plots <- function(
  explore,
  x,
  g,
  subsetter,
  splittable_vars,
  explore_min_depth,
  explore_min_height,
  explore_min_size,
  plot_list,
  split_num,
  missed_splits,
  signs
) {

  if (explore) {
    explore_valleys <- propose_valleys(
      x, subsetter[, g], !splittable_vars[g, ],
      explore_min_depth, explore_min_height
    )

    explore_splits <- 0

    explore_check <- !is.na(explore_valleys[1, ])
    explore_check <- explore_check & sum(subsetter[, g]) >= explore_min_size

    for (p in which(explore_check)) {
      explore_splits <- explore_splits + 1
      plot_list[[g]][[split_num[g] + missed_splits + explore_splits]] <-
        plot_gatetree_split(
          x[subsetter[, g], p], g, p, score = explore_valleys[2, p],
          signs, scenario = "explore", split_gp = explore_valleys[1, p]
        )
    }
  }

  return(plot_list)
}
