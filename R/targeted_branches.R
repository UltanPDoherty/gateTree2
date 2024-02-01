#' Partition data set based on marginal splits described by +/- table.
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
#' @param explore Logical value.
#'
#' @return List: splits, plusminus_table, subsetter, plot_list
#' @importFrom ggpubr ggarrange
#' @export
targeted_branches <- function(
  x,
  plusminus_table,
  min_height = 0.5,
  min_depth = 0.5,
  min_val_cutoff = NULL,
  max_val_cutoff = NULL,
  use_boundaries = TRUE,
  explore = TRUE
) {

  path_num <- nrow(plusminus_table)
  var_num <- ncol(plusminus_table)
  obs_num <- nrow(x)

  splits <- scores <- array(dim = dim(plusminus_table))
  plot_list <- vector("list", path_num)

  subsetter <- matrix(TRUE, nrow = obs_num, ncol = path_num)
  colnames(subsetter) <- rownames(plusminus_table)

  no_valleys <- already_split <- array(FALSE, dim = dim(plusminus_table))
  no_valleys[plusminus_table == 0] <- NA
  already_split[plusminus_table == 0] <- NA

  inside_cutoffs <- find_inside_cutoffs(x, min_val_cutoff, max_val_cutoff)

  # loop over the gating pathways
  for (g in 1:path_num){
    # exclude observations that are outside the cutoffs for any variable used in
    # this gating pathway
    subsetter[, g] <- apply(inside_cutoffs[, plusminus_table[g, ] != 0], 1, all)

    # continue inner loop until this pathway's row of the already_split & no_valleys
    # matrices do not contain any FALSE values.
    # that is, move onto the next pathway only when every required variable for
    # this pathway is TRUE
    while (any(!already_split[g, ] & !no_valleys[g, ], na.rm = TRUE)) {

      proposals <- propose_valleys(x, g, var_num, subsetter, already_split,
                                   min_depth, min_height)
      found_valley <- any(!is.na(proposals[1, ]))

      if (found_valley) {
        found_boundary <- NA
        scenario <- "valley"
      } else if (use_boundaries) {
        proposals <- propose_boundaries(x, g, var_num, subsetter, already_split)
        found_boundary <- any(!is.na(proposals[1, ]))
        scenario <- "boundary"
      } else {
        found_boundary <- FALSE
        scenario <- "nothing"
      }

      if (found_valley || found_boundary) {
        # find the proposed valley with the highest score
        p_choice <- which.max(proposals[2, ])
        # put the valley and its score into the splits and scores matrices
        splits[g, p_choice] <- proposals[1, p_choice]
        scores[g, p_choice] <- proposals[2, p_choice]

        already_split[g, p_choice] <- TRUE
        x_gp <- x[subsetter[, g], p_choice]
        # find which observations have values less than the chosen split
        less_gp <- x[, p_choice] < splits[g, p_choice]
        is_neg_gp <- plusminus_table[g, p_choice] == -1
        refine_subset <- (is_neg_gp & less_gp) | (!is_neg_gp & !less_gp)
        subsetter[, g] <- subsetter[, g] & refine_subset

        g_length <- length(plot_list[[g]])
        plot_list[[g]][[g_length + 1]] <- plot_targeted_split(
          x_gp, g, p_choice, scores[g, p_choice],
          plusminus_table, scenario, splits[g, p_choice]
        )
      } else {
        no_valleys[g, !is.na(already_split[g, ]) & !already_split[g, ]] <- TRUE
        for (p in which(!is.na(already_split[g, ]) & !already_split[g, ])) {
          x_gp <- x[subsetter[, g], p]
          g_length <- length(plot_list[[g]])
          plot_list[[g]][[g_length + 1]] <- plot_targeted_split(
            x_gp, g, p, scores[g, p],
            plusminus_table, scenario, splits[g, p]
          )
        }
      }
    }

    if (explore) {
      false_already_split <- array(FALSE, dim = dim(plusminus_table))

      proposals <- propose_valleys(x, g, var_num, subsetter, false_already_split,
                                   2 * min_depth, 2 * min_height)
      for (p in 1:var_num) {
        if (!is.na(proposals[1, p])) {
          x_gp <- x[subsetter[, g], p]
          if (length(x_gp) > 100) {
            scenario <- "undiscovered"

            g_length <- length(plot_list[[g]])
            plot_list[[g]][[g_length + 1]] <- plot_targeted_split(
              x_gp, g, p, proposals[2, p],
              plusminus_table, scenario, proposals[1, p]
            )
          }
        }
      }
    }
  }

  plot_paths(plot_list)

  is_a_duplicate <- check_duplicates(subsetter)
  subsetter <- subsetter[, !is_a_duplicate, drop = FALSE]
  already_split <- already_split[!is_a_duplicate, , drop = FALSE]
  splits <- splits[!is_a_duplicate, , drop = FALSE]
  scores <- scores[!is_a_duplicate, , drop = FALSE]
  no_valleys <- no_valleys[!is_a_duplicate, , drop = FALSE]
  plusminus_table <- plusminus_table[!is_a_duplicate, , drop = FALSE]
  path_num <- path_num - sum(is_a_duplicate)

  plusminus_table <- plusminus_table * !is.na(splits)

  return(list(splits = splits,
              plusminus_table = plusminus_table,
              subsetter = subsetter,
              plot_list = plot_list))
}

