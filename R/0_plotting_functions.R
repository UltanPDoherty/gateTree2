# ==============================================================================

#' @title Density plot with a `gateTree` split.
#'
#' @description
#' Plot a univariate kernel density estimate of the data with the `gateTree`
#' split illustrated, and information about the split in the title and subtitle.
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
    "valley"   = paste0("depth % = ", round(score, 1)),
    "boundary" = paste0("sc. BIC diff. = ", round(score, 3)),
    "nothing"  = NA,
    "explore"  = paste0("depth % = ", round(score, 1))
  )
  line_type <- switch(scenario,
    "valley"   = "solid",
    "boundary" = "dashed",
    "nothing"  = "blank",
    "explore"  = "dotted"
  )
  is_negative <- switch(scenario,
    "valley"   = (plusminus_table[g, p] == -1),
    "boundary" = (plusminus_table[g, p] == -1),
    "nothing"  = NA,
    "explore"  = FALSE
  )

  dens_gp <- stats::density(x_gp)

  min_x_gp <- min(x_gp)
  max_x_gp <- max(x_gp)
  if (scenario == "nothing") {
    xleft <- min_x_gp
    xright <- max_x_gp
  } else {
    xleft <- ifelse(is_negative, min_x_gp, split_gp)
    xright <- ifelse(is_negative, split_gp, max_x_gp)
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

  gg <- ggplot2::ggplot(
    dens_gp_df, ggplot2::aes(x = dens_gp_x, y = dens_gp_y)
  ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xleft, xmax = xright, ymin = 0, ymax = max(dens_gp_y)
      ),
      fill = rect_col, na.rm = TRUE
    ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = paste0("g = ", g, ", p = ", p, ", ", score_title),
      subtitle = paste0(
        "Path: ", rownames(plusminus_table)[g], ", ",
        "Var: ", colnames(plusminus_table)[p], ", ",
        "\"", scenario, "\""
      ),
      x = paste0("N before = ", size_before, ", N after = ", size_after),
      y = "Density %"
    ) +
    ggplot2::theme_bw()

  if (!is.na(split_gp)) {
    gg <- gg + ggplot2::geom_vline(
      xintercept = split_gp, linetype = line_type,
      na.rm = TRUE
    )
  }

  gg
}
