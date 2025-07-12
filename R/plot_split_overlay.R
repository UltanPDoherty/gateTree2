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
plot_split_overlay <- function(
    matrices, gatetree_out, var, pop, batch = NULL, samp = NULL) {
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

  gatetree_call <- gatetree_out$call
  gatetree_out <- gatetree_out$output

  if (!is.null(names(gatetree_out[[pop]]$pm_previous))) {
    var_name <- names(gatetree_out[[pop]]$pm_previous)[var]
  } else {
    var_name <- var
  }
  if (!is.null(names(gatetree_out))) {
    pop_name <- names(gatetree_out)[pop]
  } else {
    pop_name <- pop
  }

  if (is.na(gatetree_out[[pop]]$method[var])) {
    split_num <- max(gatetree_out[[pop]]$order, na.rm = TRUE)
  } else {
    split_num <- gatetree_out[[pop]]$order[var]
  }

  dens_dfs <- list()
  if (!is.null(var) && !is.null(pop)) {
    batch_num <- length(matrices)
    samp_num <- vapply(matrices, length, integer(1L))
    samp_count <- 1
    for (b in seq_len(batch_num)) {
      for (j in seq_len(samp_num[b])) {
        dens_dfs[[samp_count]] <- dens_single_split(
          matrices, gatetree_out, gatetree_call,
          pop, b, j, var, split_num
        )
        samp_count <- samp_count + 1
      }
    }
  } else {
    stop("At least pop and var are required.")
  }

  batch_names <- names(matrices)
  samp_names <- lapply(matrices, names)
  dens_df <- Reduce(rbind, dens_dfs)
  dens_df$batch_name <- batch_names[dens_df$batch]
  dens_df$samp_name <- sapply(
    seq_len(nrow(dens_df)),
    \(x) samp_names[[c(dens_df$batch[x], dens_df$samp[x])]]
  )

  dens_x <- dens01_y <- samp_name <- NULL
  gg <- dens_df |>
    ggplot2::ggplot(ggplot2::aes(
      x = dens_x, y = dens01_y, group = samp_name, colour = samp_name
    )) +
    ggplot2::geom_line(show.legend = FALSE) +
    ggplot2::labs(
      title = paste0(
        "Population: ", pop_name, ", ",
        "Variable: ", var_name
      ),
      x = var_name,
      y = "Density (0-1 Scaled)"
    ) +
    ggplot2::theme_bw()

  if (!all(is.na(dens_df$split_val))) {
    split_val <- NULL
    gg <- gg + ggplot2::geom_vline(
      ggplot2::aes(xintercept = split_val, colour = samp_name),
      na.rm = TRUE, show.legend = FALSE
    )
  }

  min_cut <- gatetree_call$min_split_cutoffs[var]
  max_cut <- gatetree_call$max_split_cutoffs[var]
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

dens_single_split <- function(
    matrices, gatetree_out, gatetree_call, pop, batch, samp, var, split_num) {
  x <- matrices[[batch]][[samp]][
    gatetree_out[[pop]]$subsetter[[batch]][[samp]][, split_num],
  ]

  x <- x[, var]
  min_cut <- gatetree_call$min_split_cutoffs[var]
  max_cut <- gatetree_call$max_split_cutoffs[var]
  x <- x[x > min_cut]
  x <- x[x < max_cut]

  split_val <- gatetree_out[[pop]]$splits[[batch]][samp, var]

  dens <- stats::density(x)
  dens_x <- dens$x
  dens01_y <- dens$y / max(dens$y)
  dens_df <- data.frame(
    batch, samp, dens_x, dens01_y,
    "split_val" = as.numeric(split_val)
  )

  dens_df
}
