#===============================================================================

#' Wrapper for `find_valley`.
#'
#' @param x Data.
#' @param subsetter_g Column g of the subsetting matrix.
#' @param splittable_vars_g Row g of the splittable_vars matrix.
#' @param min_depth `min_depth` for `find_valley` function.
#' @param min_height `min_height` for `find_valley` function.
#'
#' @return valleys
#' @export
propose_valleys <- function(
  x,
  subsetter_g,
  splittable_vars_g = rep(TRUE, ncol(x)),
  min_depth, min_height
) {
  valleys <- matrix(nrow = 2, ncol = ncol(x))

  # loop over all variables to propose splits
  for (p in which(splittable_vars_g)){
    scale01_gp <- scale01(x[subsetter_g, p])
    dens01_gp <- stats::density(scale01_gp$y)

    valleys[, p] <- find_valley(
      dens01_gp,
      depth = TRUE,
      min_depth = min_depth,
      min_height = min_height
    )

    valleys[1, p] <- unscale01(valleys[1, p], scale01_gp$min, scale01_gp$max)
  }

  return(valleys)
}

#===============================================================================

#' Wrapper for `find_boundary`.
#'
#' @param x Data.
#' @param g The pathway number.
#' @param var_num Total number of variables.
#' @param subsetter The subsetting matrix.
#' @param already_split The already_split matrix.
#' @param common_variables Matrix
#'
#' @return boundaries
#' @export
propose_boundaries <- function(x, g, var_num, subsetter, already_split,
                               common_variables = NULL) {
  boundaries <- matrix(nrow = 2, ncol = var_num)

  if (is.null(common_variables)) {
    var_set <- 1:var_num
  } else {
    var_set <- which(common_variables[g, ])
  }

  # loop over all variables to propose splits
  for (p in var_set){
    # find variables that are not NA or TRUE in the already_split matrix
    if (!is.na(already_split[g, p]) && !already_split[g, p]) {
      boundaries[, p] <- find_boundary(x[subsetter[, g], p])
    }
  }

  return(boundaries)
}

#===============================================================================

#' 0-1 scaling.
#'
#' @param x Data.
#' @param other_min Minimum to be used for 0-1 scaling.
#' @param other_max Maximum to be used for 0-1 scaling.
#'
#' @return Scaled version of `x`.
#' @export
scale01 <- function(x, other_min = NULL, other_max = NULL) {

  if (is.null(other_min)) {
    minimum <- min(x)
  } else {
    minimum <- other_min
  }

  if (is.null(other_max)) {
    maximum <- max(x)
  } else {
    maximum <- other_max
  }

  y <- (x - minimum) / (maximum - minimum)

  return(list(y = y, min = minimum, max = maximum))
}

#===============================================================================

#' Undo 0-1 scaling.
#'
#' @param x Data.
#' @param unscaled_min Minimum used for 0-1 scaling.
#' @param unscaled_max Maximum used for 0-1 scaling.
#'
#' @return Unscaled version of `x`.
#' @export
unscale01 <- function(x, unscaled_min, unscaled_max) {
  return(x * (unscaled_max - unscaled_min) + unscaled_min)
}

#===============================================================================

#' Check if final subsets are identical.
#'
#' @param subsetter The subsetting matrix.
#'
#' @return `is_a_duplicate`.
#' @export
check_duplicates <- function(subsetter) {

  path_num <- ncol(subsetter)

  is_a_duplicate <- rep(FALSE, path_num)

  if (path_num > 1) {
    equal_subsets <- matrix(nrow = path_num, ncol = path_num)
    for (g in 1:(path_num - 1)) {
      for (h in (g + 1):path_num) {
        equal_subsets[g, h] <- all(subsetter[, g] == subsetter[, h])
        if (equal_subsets[g, h]) {
          print(paste0("Failed to distinguish between populations ",
                       g, " & ", h, "."))
          is_a_duplicate[h] <- TRUE
        }
      }
    }
  }

  return(is_a_duplicate)
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

#' Supervised merging.
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


#===============================================================================

#' Find which observations are inside the cutoffs for each variable.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#'
#' @return inside_cutoffs
#' @export
find_inside_cutoffs <- function(x, min_val_cutoff, max_val_cutoff) {
  # find which observations are outside either of the cutoffs for each variable
  if (is.null(min_val_cutoff)) {
    below_cutoff <- array(FALSE, dim = dim(x))
  } else {
    below_cutoff <- array(dim = dim(x))
    for (j in seq_len(ncol(x))) {
      below_cutoff[, j] <- x[, j] < min_val_cutoff[j]
    }
  }
  if (is.null(max_val_cutoff)) {
    above_cutoff <- array(FALSE, dim = dim(x))
  } else {
    above_cutoff <- array(dim = dim(x))
    for (j in seq_len(ncol(x))) {
      above_cutoff[, j] <- x[, j] > max_val_cutoff[j]
    }
  }
  # inside_cutoffs is a matrix of the same dimensions as the data
  # an observation-variable pair is TRUE if it is inside the two cutoffs
  inside_cutoffs <- !below_cutoff & !above_cutoff

  return(inside_cutoffs)
}

#===============================================================================

#' Plotting function for `targeted_split`.
#'
#' @param x_gp The data to be displayed, should be only pathway g and variable p.
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

  rect_col <- switch(scenario,
                     "valley" = "lightgreen",
                     "boundary" = "lightblue",
                     "nothing" = NA,
                     "undiscovered" = "lightcoral")
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
#' @importFrom ggokabeito scale_colour_okabe_ito
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
    ggokabeito::scale_colour_okabe_ito(order = c(9, 1:8))

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
