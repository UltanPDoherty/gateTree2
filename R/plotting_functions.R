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

#' Plotting function for `targeted_split`.
#'
#' @param x_gp The data to be displayed, should be only pathway g & variable p.
#' @param g The pathway number.
#' @param p The variable number.
#' @param depth The depth of the split.
#' @param plusminus_table The cell-type variable table.
#' @param scenario "valley", "boundary", "nothing", or "undiscovered".
#' @param split_gp The split value.
#'
#' @return NULL
#' @export
plot_targeted_split <- function(x_gp, g, p, depth, plusminus_table,
                                scenario, split_gp) {

  # colours from ggokabeito package
  rect_col <- switch(scenario,
                     "valley" = "#F0E442",
                     "boundary" = "#56B4E9",
                     "nothing" = NA,
                     "undiscovered" = "#CC79A7")
  depth <- switch(scenario,
                  "valley" = depth,
                  "boundary" = NA,
                  "nothing" = NA,
                  "undiscovered" = depth)
  line_type <- switch(scenario,
                      "valley" = "solid",
                      "boundary" = "dashed",
                      "nothing" = "blank",
                      "undiscovered" = "dotted")
  is_negative <- switch(scenario,
                        "valley" = (plusminus_table[g, p] == -1),
                        "boundary" = (plusminus_table[g, p] == -1),
                        "nothing" = NA,
                        "undiscovered" = FALSE)

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
      title = paste0("g = ", g, ", p = ", p, ", ",
                     "depth = ", round(depth, 1), "%"),
      subtitle = paste0("Pathway: ", rownames(plusminus_table)[g], ", ",
                        "Variable: ", colnames(plusminus_table)[p]),
      x = paste0("N before = ", size_before, ", N after = ", size_after),
      y = "Density %"
    ) +
    theme_bw()

  return(gg)
}

#===============================================================================

make_edge_df <- function(parent_node, node_number, edge_name, node_name,
                         is_leaf, path_nodes) {

  leaf_name <- rep("", node_number)
  leaf_name[is_leaf] <- node_name[is_leaf]
  leaf_name <- gsub("All/", "", leaf_name)

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

#' @import ggraph
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
make_tree_plot <- function(edge_df) {

  tree_graph <- igraph::graph_from_data_frame(d = edge_df[, 1:3],
                                              v = edge_df[, c(2, 4:7)])

  tree_tbl_graph <- tidygraph::as_tbl_graph(tree_graph)

  edge_name <- leaf_name <- NULL # to silence "no visible binding" check
  tree_plot <- tree_tbl_graph |>
    ggraph::ggraph(layout = "tree") +
    ggraph::geom_edge_elbow2(
      aes(label = edge_name),
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

  return(tree_plot)
}

#===============================================================================

plot_paths <- function(plot_list) {
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
