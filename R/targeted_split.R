#' Partition data set based on marginal splits described by +/- table.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param typemarker Cell-type marker table.
#' @param min_height Minimum height for a peak to be recognised by find_peaks.
#' @param min_depth Minimum depth for a split to be returned by find_valley.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#' @param use_boundaries Logical value.
#' @param explore Logical value.
#'
#' @return splits, typemarker, subsetter
#' @export
targeted_split <- function(
  x,
  typemarker,
  min_height = 0.5,
  min_depth = 0.5,
  min_val_cutoff = NULL,
  max_val_cutoff = NULL,
  use_boundaries = TRUE,
  explore = TRUE
) {

  path_num <- nrow(typemarker)
  var_num <- ncol(typemarker)
  obs_num <- nrow(x)

  splits <- scores <- array(dim = dim(typemarker))

  subsetter <- matrix(TRUE, nrow = obs_num, ncol = path_num)
  paused <- array(FALSE, dim = dim(typemarker))
  progress <- typemarker == 0

  inside_cutoffs <- find_inside_cutoffs(x, min_val_cutoff, max_val_cutoff)

  # loop over the gating pathways
  for (g in 1:path_num){
    # exclude observations that are outside the cutoffs for any marker used in
    # this gating pathway
    subsetter[, g] <- apply(inside_cutoffs[, typemarker[g, ] != 0], 1, all)

    graphics::par(mfrow = c(2, 2))
    # continue inner loop until this pathway's row of the progress & paused
    # matrices do not contain any FALSE values.
    # that is, move onto the next pathway only when every required variable for
    # this pathway is TRUE
    while (any(!progress[g, ] & !paused[g, ], na.rm = TRUE)) {

      proposals <- propose_valleys(x, g, var_num, subsetter, progress,
                                   min_depth, min_height)
      found_valley <- any(!is.na(proposals[1, ]))

      if (found_valley) {
        found_boundary <- NA
        scenario <- "valley"
      } else if (use_boundaries) {
        proposals <- propose_boundaries(x, g, var_num, subsetter, progress)
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

        progress[g, p_choice] <- TRUE
        x_gp <- x[subsetter[, g], p_choice]
        # find which events have values less than the chosen split
        less_gp <- x[, p_choice] < splits[g, p_choice]
        is_neg_gp <- typemarker[g, p_choice] == -1
        refine_subset <- (is_neg_gp & less_gp) | (!is_neg_gp & !less_gp)
        subsetter[, g] <- subsetter[, g] & refine_subset

        split_size <- c(length(x_gp), sum(subsetter[, g]))

        plot_targeted_split(x_gp, g, p_choice,
                            scores[g, p_choice], typemarker,
                            scenario,
                            splits[g, p_choice], split_size)
      } else {
        paused[g, !is.na(progress[g, ]) & !progress[g, ]] <- TRUE
        for (p in which(!is.na(progress[g, ]) & !progress[g, ])) {
          x_gp <- x[subsetter[, g], p]
          split_size <- c(length(x_gp), NA)
          plot_targeted_split(x_gp, g, p,
                              scores[g, p], typemarker,
                              scenario,
                              splits[g, p], split_size)
        }
      }
    }

    if (explore) {
      false_progress <- array(FALSE, dim = dim(typemarker))

      proposals <- propose_valleys(x, g, var_num, subsetter, false_progress,
                                   2 * min_depth, 2 * min_height)
      for (p in 1:var_num) {
        if (!is.na(proposals[1, p])) {
          x_gp <- x[subsetter[, g], p]
          split_size <- c(length(x_gp), sum(x_gp > proposals[1, p]))
          if (split_size[1] > 100) {
            scenario <- "undiscovered"
            plot_targeted_split(x_gp, g, p,
                                proposals[2, p], typemarker,
                                scenario,
                                proposals[1, p], split_size)
          }
        }
      }
    }
  }

  is_a_duplicate <- check_duplicates(subsetter)
  subsetter <- subsetter[, !is_a_duplicate]
  progress <- progress[!is_a_duplicate, ]
  splits <- splits[!is_a_duplicate, ]
  scores <- scores[!is_a_duplicate, ]
  paused <- paused[!is_a_duplicate, ]
  typemarker <- typemarker[!is_a_duplicate, ]
  path_num <- path_num - sum(is_a_duplicate)

  typemarker <- typemarker * !is.na(splits)

  return(list(splits = splits,
              typemarker = typemarker,
              subsetter = subsetter))
}

#===============================================================================

#' Find which observations are inside the cutoffs for each marker.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#'
#' @return inside_cutoffs
#' @export
find_inside_cutoffs <- function(x, min_val_cutoff, max_val_cutoff) {
  # find which observations are outside either of the cutoffs for each marker
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
  # an observation-marker pair is TRUE if it is inside the two cutoffs
  inside_cutoffs <- !below_cutoff & !above_cutoff

  return(inside_cutoffs)
}

#===============================================================================

#' Plotting function for `targeted_split`.
#'
#' @param x_gp The data to be displayed, should be only pathway g and marker p.
#' @param g The pathway number.
#' @param p The marker number.
#' @param depth The depth of the split.
#' @param typemarker The cell-type marker table.
#' @param scenario "valley", "boundary", "nothing", or "undiscovered".
#' @param split_gp The split value.
#' @param split_size The size of the subset before and after the split.
#'
#' @return NULL
#' @export
plot_targeted_split <- function(x_gp, g, p, depth, typemarker,
                                scenario, split_gp, split_size) {

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
  linetype <- switch(scenario,
                     "valley" = "solid",
                     "boundary" = "dashed",
                     "nothing" = "blank",
                     "undiscovered" = "dotted")
  is_negative <- switch(scenario,
                        "valley" = (typemarker[g, p] == -1),
                        "boundary" = (typemarker[g, p] == -1),
                        "nothing" = NA,
                        "undiscovered" = FALSE)

  scale01_gp <- scale01(x_gp)
  dens_gp <- stats::density(scale01_gp$y)

  trans_split_gp <- scale01(split_gp,
                            scale01_gp$min, scale01_gp$max)$y

  xleft <- ifelse(is_negative, 0, trans_split_gp)
  xright <- ifelse(is_negative, trans_split_gp, 1)

  plot(dens_gp,
       main = paste0("g = ", g, ", p = ", p,
                     ", depth = ", round(depth, 3)),
       sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p]),
       xlab = paste0("N before = ", split_size[1], ", ",
                     "N after = ", split_size[2]),
       panel.first = graphics::rect(xleft, 0, xright, max(dens_gp$y),
                                    col = rect_col, border =  NA),)
  graphics::abline(v = trans_split_gp, lty = linetype)

  return(NULL)
}

#===============================================================================

#' Wrapper for `find_valley`.
#'
#' @param x Data.
#' @param g The pathway number.
#' @param var_num Total number of markers.
#' @param subsetter The subsetting matrix.
#' @param progress The progress matrix.
#' @param min_score min_score for find_valley function.
#' @param min_height min_height for find_valley function.
#'
#' @return valleys
#' @export
propose_valleys <- function(x, g, var_num, subsetter, progress,
                            min_score, min_height) {
  valleys <- matrix(nrow = 2, ncol = var_num)

  # loop over all variables to propose splits
  for (p in 1:var_num){
    # find variables that are not NA or TRUE in the progress matrix
    if (!is.na(progress[g, p]) && !progress[g, p]) {
      # 0-1 scale this variable for the pathway's current subset
      scale01_gp <- scale01(x[subsetter[, g], p])
      dens01_gp <- stats::density(scale01_gp$y)

      valleys[, p] <- find_valley(
        dens01_gp,
        score = TRUE,
        min_score = min_score,
        min_height = min_height
      )

      valleys[1, p] <- unscale01(valleys[1, p], scale01_gp$min, scale01_gp$max)
    }
  }

  return(valleys)
}

#===============================================================================

#' Wrapper for `find_boundary`.
#'
#' @param x Data.
#' @param g The pathway number.
#' @param var_num Total number of markers.
#' @param subsetter The subsetting matrix.
#' @param progress The progress matrix.
#'
#' @return boundaries
#' @export
propose_boundaries <- function(x, g, var_num, subsetter, progress) {
  boundaries <- matrix(nrow = 2, ncol = var_num)

  # loop over all variables to propose splits
  for (p in 1:var_num){
    # find variables that are not NA or TRUE in the progress matrix
    if (!is.na(progress[g, p]) && !progress[g, p]) {
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
