#' @title gateTree2 tree plot.
#'
#' @inheritParams plot_split
#' @param x_pad x-axis padding around plot via `ggplot2::coord_cartesian`.
#' @param y_pad y-axis padding around plot via `ggplot2::coord_cartesian`.
#' @param pm_pad Padding from vertical line to the plus / minus symbols.
#'
#' @return `ggplot` object.
#'
#' @export
plot_tree <- function(gatetree_out, x_pad = 0.2, y_pad = 0.2, pm_pad = 0.05) {
  if (!is.null(gatetree_out$output)) {
    gatetree_out <- gatetree_out$output
  }

  pop_num <- length(gatetree_out)
  pop_names <- names(gatetree_out)

  pm_list <- lapply(gatetree_out, \(x) x$pm_previous)
  order_list <- lapply(gatetree_out, \(x) x$order)
  node_type_list <- lapply(gatetree_out, \(x) x$method)

  max_splits <- vapply(order_list, max, double(1L), na.rm = TRUE)

  tree_list <- list()
  tree_list$Population <- c()
  tree_list$Split <- c()
  tree_list$node_label <- c()
  tree_list$Sign <- c()
  tree_list$node_type <- c()
  tree_count <- 1
  for (p in seq_len(pop_num)) {
    for (i in seq_len(max_splits[p])) {
      var <- which(order_list[[p]] == i)

      tree_list$Population[tree_count] <- p
      tree_list$Split[tree_count] <- i
      tree_list$node_label[tree_count] <- names(order_list[[p]])[var]
      tree_list$Sign[tree_count] <- pm_list[[p]][var]
      tree_list$node_type[tree_count] <- node_type_list[[p]][var]

      tree_count <- tree_count + 1
    }
  }
  tree_df <- as.data.frame(tree_list)
  node_num <- nrow(tree_df)

  pop_x <- pop_y <- double(pop_num)
  old_pop_x <- old_pop_y <- double(pop_num)
  for (j in seq_len(node_num)) {
    pop_j <- tree_df$Population[j]
    sign_j <- tree_df$Sign[j]
    split_j <- tree_df$Split[j]

    x_shift <- 1 / 2^split_j

    tree_df$node_x[j] <- pop_x[pop_j]
    tree_df$node_y[j] <- pop_y[pop_j]
    tree_df$last_node_x[j] <- old_pop_x[pop_j]
    tree_df$last_node_y[j] <- old_pop_y[pop_j]

    old_pop_x[pop_j] <- pop_x[pop_j]
    old_pop_y[pop_j] <- pop_y[pop_j]

    pop_x[pop_j] <- pop_x[pop_j] + sign_j * x_shift
    pop_y[pop_j] <- pop_y[pop_j] + 1
  }
  tree_df$before_x <- tree_df$last_node_x
  tree_df$before_y <- tree_df$node_y - (tree_df$Split != 1) * 0.5
  tree_df$after_y <- tree_df$node_y + 0.5

  for (p in seq_len(pop_num)) {
    tree_df <- rbind(tree_df, rep(NA, ncol(tree_df)))

    tree_df$Population[node_num + p] <- p
    tree_df$Split[node_num + p] <- NA
    tree_df$node_label[node_num + p] <- pop_names[p]
    tree_df$Sign[node_num + p] <- NA
    tree_df$node_type[node_num + p] <- "leaf"

    tree_df$node_x[node_num + p] <- pop_x[p]
    tree_df$node_y[node_num + p] <- pop_y[p]
    tree_df$last_node_x[node_num + p] <- old_pop_x[p]
    tree_df$last_node_y[node_num + p] <- old_pop_y[p]

    tree_df$before_x[node_num + p] <- old_pop_x[p]
    tree_df$before_y[node_num + p] <- pop_y[p] - 0.5
    tree_df$after_y[node_num + p] <- pop_y[p]
  }

  to_be_removed <- logical(node_num + pop_num)
  for (p in seq(node_num + 1, node_num + pop_num - 1)) {
    for (q in seq(p + 1, node_num + pop_num)) {
      same_leaf <- all(tree_df[p, -c(1:5)] == tree_df[q, -c(1:5)])
      if (same_leaf) {
        tree_df$node_label[p] <- paste0(
          tree_df$node_label[p], "_", tree_df$node_label[q]
        )
        to_be_removed[q] <- TRUE
      }
    }
  }
  if (any(to_be_removed)) {
    tree_df <- tree_df[-which(to_be_removed), ]
  }

  tree_df$plus_x <- tree_df$node_x + pm_pad
  tree_df$plus_y <- (tree_df$node_y + tree_df$after_y) / 2
  tree_df$plus_text <- ifelse(tree_df$node_type != "leaf", "+", "")

  tree_df$minus_x <- tree_df$node_x - pm_pad
  tree_df$minus_y <- (tree_df$node_y + tree_df$after_y) / 2
  tree_df$minus_text <- ifelse(tree_df$node_type != "leaf", "-", "")

  flip_tree_df <- tree_df
  flip_tree_df$node_y <- -flip_tree_df$node_y
  flip_tree_df$last_node_y <- -flip_tree_df$last_node_y
  flip_tree_df$before_y <- -flip_tree_df$before_y
  flip_tree_df$after_y <- -flip_tree_df$after_y
  flip_tree_df$plus_y <- -flip_tree_df$plus_y
  flip_tree_df$minus_y <- -flip_tree_df$minus_y

  colour_scheme <- c(
    "valley" = "#F0E442", "boundary" = "#56B4E9", "leaf" = "#009E73"
  )
  x_limits <- 
    c(min(flip_tree_df$node_x) - x_pad, max(flip_tree_df$node_x) + x_pad)
  y_limits <- 
    c(min(flip_tree_df$node_y) - y_pad, max(flip_tree_df$node_y) + y_pad)
  node_x <- node_y <- node_label <- after_y <- before_y <- before_x <-
    node_type <- plus_x <- plus_y <- plus_text <- minus_x <- minus_y <-
    minus_text <- NULL
  flip_tree_df |>
    ggplot2::ggplot(ggplot2::aes(
      x = node_x, y = node_y, label = node_label
    )) +
    ggplot2::geom_linerange(
      ggplot2::aes(x = node_x, ymin = node_y, ymax = after_y)
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(x = node_x, ymin = before_y, ymax = node_y)
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = before_x, xmax = node_x, y = before_y)
    ) +
    ggplot2::geom_label(
      ggplot2::aes(fill = node_type),
      colour = "white", fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = plus_x, y = plus_y, label = plus_text),
      size = 8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = minus_x, y = minus_y, label = minus_text),
      size = 8
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) +
    ggplot2::scale_fill_manual(values = colour_scheme) +
    ggplot2::theme_void()
}
