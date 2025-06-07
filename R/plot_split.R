# ==============================================================================

#' @title Density plot with a `gateTree` split.
#'
#' @description
#' Plot a univariate kernel density estimate of the data with the `gateTree`
#' split illustrated, and information about the split in the title and subtitle.
#'
#' @inheritParams gatetree
#' @param gatetree_out Output from gatetree.
#' @param pop The population number.
#' @param samp The sample number.
#' @param var The variable number.
#'
#' @return `ggplot` object.
plot_split <- function(
    samples, gatetree_out, pop = NULL, samp = NULL, var = NULL) {
  plots <- list()
  if (!is.null(pop) && !is.null(samp) && !is.null(var)) {
    plots[[1]] <- plot_single_split(
      samples, gatetree_out, pop, samp, var, c(1, 1)
    )
  } else if (!is.null(pop) && !is.null(samp) && is.null(var)) {
    max_split <- max(gatetree_out[[pop]]$order, na.rm = TRUE)
    for (i in seq_len(max_split)) {
      var <- which(gatetree_out[[pop]]$order == i)
      plots[[i]] <- plot_single_split(
        samples, gatetree_out, pop, samp, var, c(i, max_split)
      )
    }
  } else if (!is.null(pop) && is.null(samp) && !is.null(var)) {
    samp_num <- length(samples)
    for (j in seq_len(samp_num)) {
      plots[[j]] <- plot_single_split(
        samples, gatetree_out, pop, j, var, c(j, samp_num)
      )
    }
  } else if (is.null(pop) && !is.null(samp) && !is.null(var)) {
    pop_num <- length(gatetree_out)
    for (k in seq_len(pop_num)) {
      plots[[k]] <- plot_single_split(
        samples, gatetree_out, k, samp, var, c(k, pop_num)
      )
    }
  } else {
    stop("At least two of pop, samp, and var are required.")
  }

  plots
}

plot_single_split <- function(
    samples, gatetree_out, pop, samp, var, plot_num) {
  scenario <- ifelse(
    is.na(gatetree_out[[pop]]$method[var]),
    "nothing",
    gatetree_out[[pop]]$method[var]
  )
  score <- switch(scenario,
    "valley" = gatetree_out[[pop]]$depths[[samp]][var],
    "boundary" = gatetree_out[[pop]]$diffs[[samp]][var],
    "nothing" = NA
  )

  # colours from ggokabeito package
  rect_col <- switch(scenario,
    "valley" = "#F0E442",
    "boundary" = "#56B4E9",
    "nothing" = NA,
    "explore" = "#CC79A7"
  )
  score <- switch(scenario,
    "valley"   = score,
    "boundary" = score,
    "nothing"  = NA,
    "explore"  = score
  )
  score_title <- switch(scenario,
    "valley"   = paste0("depth = ", round(score, 1), " events"),
    "boundary" = paste0("scaled BIC diff. = ", round(score, 3)),
    "nothing"  = NA,
    "explore"  = paste0("depth = ", round(score, 1))
  )
  is_negative <- switch(scenario,
    "valley"   = gatetree_out[[pop]]$pm_previous[var] == -1,
    "boundary" = gatetree_out[[pop]]$pm_previous[var] == -1,
    "nothing"  = NA,
    "explore"  = FALSE
  )

  split_num <- gatetree_out[[pop]]$order[var]
  if (is.na(split_num)) {
    split_num <- max(gatetree_out[[pop]]$order, na.rm = TRUE)
  }
  x <- samples[[samp]][gatetree_out[[pop]]$subsetter[[samp]][, split_num], ]

  if (!is.null(names(samples))) {
    samp_name <- names(samples)[samp]
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

  split_val <- gatetree_out[[pop]]$splits[[samp]][var]

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

  nbin <- 100
  x_ash <- ash::bin1(x, nbin = nbin)
  counts <- x_ash$nc
  breaks <- seq(x_ash$ab[1], x_ash$ab[2], length.out = nbin + 1)
  midpoints <- (breaks[1:nbin] + breaks[-1]) / 2
  width <- (x_ash$ab[2] - x_ash$ab[1]) / nbin

  hist_x <- midpoints
  hist_y <- counts
  hist_df <- data.frame(hist_x, hist_y)

  gg <- ggplot2::ggplot(
    hist_df, ggplot2::aes(x = hist_x, y = hist_y)
  ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xleft, xmax = xright, ymin = 0, ymax = max(counts)
      ),
      fill = rect_col, na.rm = TRUE
    ) +
    ggplot2::geom_col(width = width) +
    ggplot2::labs(
      title = paste0(
        "(", plot_num[1], "/", plot_num[2], "). ",
        "Sample: ", samp_name, ", ", score_title
      ),
      subtitle = paste0(
        "Population: ", pop_name, ", ",
        "Variable: ", var_name
      ),
      x = paste0("N before = ", size_before, ", N after = ", size_after),
      y = "# of Events"
    ) +
    ggplot2::theme_bw()

  if (!is.na(split_val)) {
    gg <- gg + ggplot2::geom_vline(xintercept = split_val, na.rm = TRUE)
  }

  gg
}
