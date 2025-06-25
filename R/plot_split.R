#' @title Density plot with a `gateTree` split.
#'
#' @description
#' Plot a univariate density of the data with the `gateTree` split
#' illustrated, and information about the split in the title and subtitle.
#'
#' @inheritParams gatetree
#' @param gatetree_out Output from gatetree.
#' @param pop The population number or name.
#' @param batch The batch number or name.
#' @param samp The sample number or name.
#' @param var The variable number or name.
#'
#' @return `ggplot` object.
#'
#' @export
plot_split <- function(
    matrices, gatetree_out, pop = NULL, batch = NULL, samp = NULL, var = NULL) {
  if (!is.null(pop) && is.character(pop)) {
    if (pop %in% names(gatetree_out$output)) {
      pop <- which(pop == names(gatetree_out$output))
    }
  }
  if (!is.null(batch) && is.character(batch)) {
    if (batch %in% names(matrices)) {
      batch <- which(batch == names(matrices))
    }
  }
  if (!is.null(samp) && is.character(samp)) {
    if (samp %in% names(matrices[[batch]])) {
      samp <- which(samp == names(matrices[[batch]]))
    }
  }
  if (!is.null(var) && is.character(var)) {
    if (var %in% colnames(matrices[[1]][[1]])) {
      var <- which(var == colnames(matrices[[1]][[1]]))
    }
  }

  plots <- list()
  if (!is.null(pop) && !is.null(batch) && !is.null(samp) && !is.null(var)) {
    plots[[1]] <- plot_single_split(
      matrices, gatetree_out, pop, batch, samp, var, c(1, 1)
    )
  } else if (
    !is.null(pop) && !is.null(batch) && !is.null(samp) && is.null(var)
  ) {
    max_split <- max(gatetree_out$output[[pop]]$order, na.rm = TRUE)
    for (i in seq_len(max_split)) {
      var <- which(gatetree_out$output[[pop]]$order == i)
      plots[[i]] <- plot_single_split(
        matrices, gatetree_out, pop, batch, samp, var, c(i, max_split)
      )
    }
  } else if (
    !is.null(pop) && !is.null(batch) && is.null(samp) && !is.null(var)
  ) {
    samp_num <- length(matrices[[batch]])
    for (j in seq_len(samp_num)) {
      plots[[j]] <- plot_single_split(
        matrices, gatetree_out, pop, batch, j, var, c(j, samp_num)
      )
    }
  } else if (
    !is.null(pop) && is.null(batch) && is.null(samp) && !is.null(var)
  ) {
    batch_num <- length(matrices)
    samp_num <- vapply(matrices, length, integer(1L))
    samp_count <- 1
    for (b in seq_len(batch_num)) {
      for (j in seq_len(samp_num[b])) {
        plots[[samp_count]] <- plot_single_split(
          matrices, gatetree_out, pop, b, j, var, c(samp_count, sum(samp_num))
        )
        samp_count <- samp_count + 1
      }
    }
  } else if (
    is.null(pop) && !is.null(batch) && !is.null(samp) && !is.null(var)
  ) {
    pop_num <- length(gatetree_out$output)
    for (k in seq_len(pop_num)) {
      plots[[k]] <- plot_single_split(
        matrices, gatetree_out, k, batch, samp, var, c(k, pop_num)
      )
    }
  } else {
    stop(paste0(
      "At least three of pop, batch, samp, and var are required",
      " (or only pop and var)."
    ))
  }

  plots
}

#' @title Density plot with a `gateTree` split.
#'
#' @description
#' Plot a univariate density of the data with the `gateTree` split
#' illustrated, and information about the split in the title and subtitle.
#'
#' @inheritParams gatetree
#' @inheritParams plot_split
#'
#' @return `ggplot` object.
#'
#' @export
plot_explore <- function(
    matrices, gatetree_out, pop = NULL, batch = NULL, samp = NULL, var = NULL) {
  plots <- list()
  if (!is.null(pop) && !is.null(batch) && !is.null(samp) && !is.null(var)) {
    plots[[1]] <- plot_single_split(
      matrices, gatetree_out, pop, batch, samp, var, c(1, 1), TRUE
    )
  } else if (
    !is.null(pop) && !is.null(batch) && !is.null(samp) && is.null(var)
  ) {
    var_num <- ncol(matrices[[1]][[1]])
    for (i in seq_len(var_num)) {
      plots[[i]] <- plot_single_split(
        matrices, gatetree_out, pop, batch, samp, i, c(i, var_num), TRUE
      )
    }
  } else if (
    !is.null(pop) && !is.null(batch) && is.null(samp) && !is.null(var)
  ) {
    samp_num <- length(matrices)
    for (j in seq_len(samp_num)) {
      plots[[j]] <- plot_single_split(
        matrices, gatetree_out, pop, batch, j, var, c(j, samp_num), TRUE
      )
    }
  } else if (
    !is.null(pop) && is.null(batch) && is.null(samp) && !is.null(var)
  ) {
    batch_num <- length(matrices)
    samp_num <- vapply(matrices, length, integer(1L))
    samp_count <- 1
    for (b in seq_len(batch_num)) {
      for (j in seq_len(samp_num[b])) {
        plots[[samp_count]] <- plot_single_split(
          matrices, gatetree_out, pop, b, j, var, c(samp_count, sum(samp_num)),
          TRUE
        )
        samp_count <- samp_count + 1
      }
    }
  } else if (
    is.null(pop) && !is.null(batch) && !is.null(samp) && !is.null(var)
  ) {
    pop_num <- length(gatetree_out$output)
    for (k in seq_len(pop_num)) {
      plots[[k]] <- plot_single_split(
        matrices, gatetree_out, k, batch, samp, var, c(k, pop_num), TRUE
      )
    }
  } else {
    stop(paste0(
      "At least three of pop, batch, samp, and var are required",
      " (or only pop and var)."
    ))
  }

  plots
}

plot_single_split <- function(
    matrices, gatetree_out, pop, batch, samp, var, plot_num, explore = FALSE) {
  gatetree_call <- gatetree_out$call
  gatetree_out <- gatetree_out$output

  if (explore) {
    scenario <- "explore"
    split_num <- max(gatetree_out[[pop]]$order, na.rm = TRUE) + 1
  } else if (is.na(gatetree_out[[pop]]$method[[batch]][samp, var])) {
    scenario <- "nothing"
    split_num <- max(gatetree_out[[pop]]$order, na.rm = TRUE) + 1
  } else {
    scenario <- gatetree_out[[pop]]$method[[batch]][samp, var]
    split_num <- gatetree_out[[pop]]$order[var]
  }

  x <- matrices[[batch]][[samp]][
    gatetree_out[[pop]]$subsetter[[batch]][[samp]][, split_num],
  ]

  if (!is.null(names(matrices))) {
    batch_name <- names(matrices)[batch]
  } else {
    batch_name <- batch
  }
  if (!is.null(names(matrices[[batch]]))) {
    samp_name <- names(matrices[[batch]])[samp]
  } else {
    samp_name <- samp
  }
  if (!is.null(colnames(x))) {
    var_name <- colnames(x)[var]
  } else {
    var_name <- var
  }
  if (!is.null(names(gatetree_out))) {
    pop_name <- names(gatetree_out)[pop]
  } else {
    pop_name <- pop
  }

  x <- x[, var]
  min_cut <- gatetree_call$min_cutoffs[var]
  max_cut <- gatetree_call$max_cutoffs[var]
  x <- x[x > min_cut]
  x <- x[x < max_cut]

  if (explore) {
    valley <- find_valley(x)
    split_val <- valley[1]
    explore_valley_depth <- valley[2]
    if (is.na(valley[1])) {
      scenario <- "nothing"
    }
  } else {
    split_val <- gatetree_out[[pop]]$splits[[batch]][samp, var]
    explore_valley_depth <- NA
  }

  if (scenario == "nothing") {
    is_negative <- NA
  } else if (scenario == "explore") {
    is_negative <- FALSE
  } else {
    is_negative <- gatetree_out[[pop]]$pm_previous[var] == -1
  }

  min_x <- min(x)
  max_x <- max(x)
  if (scenario == "nothing") {
    xleft <- min_x
    xright <- max_x
  } else {
    xleft <- ifelse(is_negative, min_x, split_val)
    xright <- ifelse(is_negative, split_val, max_x)
  }

  size_before <- length(x)
  if (is.na(is_negative)) {
    size_after <- NA
  } else if (is_negative) {
    size_after <- sum(x < split_val)
  } else {
    size_after <- sum(x > split_val)
  }

  dens <- stats::density(x)
  dens_x <- dens$x
  dens01_y <- dens$y / max(dens$y)
  dens_df <- data.frame(dens_x, dens01_y)

  if (scenario == "nothing") {
    rect_col <- NA
    score <- NA
    score_title <- NA
  } else if (scenario == "explore") {
    rect_col <- "#CC79A7"
    score <- explore_valley_depth
    score_title <- paste0("valley depth = ", round(score, 3))
  } else if (scenario == "valley") {
    rect_col <- "#F0E442"
    score <- gatetree_out[[pop]]$depths[[batch]][samp, var]
    score_title <- paste0("valley depth = ", round(score, 3))
  } else if (scenario == "boundary") {
    rect_col <- "#56B4E9"
    score <- gatetree_out[[pop]]$diffs[[batch]][samp, var]
    score_title <- paste0("boundary diff = ", round(score, 3))
  } else {
    rect_col <- "#E69F00"
    score <- NA
    score_title <- scenario
  }
  # colours from ggokabeito package

  gg <- ggplot2::ggplot(
    dens_df, ggplot2::aes(x = dens_x, y = dens01_y)
  ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xleft, xmax = xright, ymin = 0, ymax = max(dens01_y)
      ),
      fill = rect_col, na.rm = TRUE
    ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = paste0(
        "(", plot_num[1], "/", plot_num[2], "). ",
        "Batch: ", batch_name, ", Sample: ", samp_name, ", ", score_title
      ),
      subtitle = paste0(
        "Population: ", pop_name, ", ",
        "Variable: ", var_name
      ),
      x = paste0("N before = ", size_before, ", N after = ", size_after),
      y = "Density (0-1 Scaled)"
    ) +
    ggplot2::theme_bw()

  if (!is.na(split_val)) {
    gg <- gg + ggplot2::geom_vline(xintercept = split_val, na.rm = TRUE)
  }
  if (is.finite(min_cut)) {
    gg <- gg + ggplot2::geom_vline(
      ggplot2::aes(xintercept = min_cut),
      linetype = "dotted", linewidth = 2, colour = "red"
    )
  }
  if (is.finite(max_cut)) {
    gg <- gg + ggplot2::geom_vline(
      ggplot2::aes(xintercept = max_cut),
      linetype = "dotted", linewidth = 2, colour = "red"
    )
  }

  gg
}
