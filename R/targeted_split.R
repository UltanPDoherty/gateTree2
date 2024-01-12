#' Partition data set based on marginal splits described by +/- table.
#'
#' @param x Dataset in matrix or data.frame form.
#' @param typemarker Cell-type marker table.
#' @param min_height Minimum height for a peak to be recognised by find_peaks.
#' @param min_score Minimum score for a split to be returned by find_valley.
#' @param min_val_cutoff Minimum value for an observation to be included.
#' @param max_val_cutoff Maximum value for an observation to be included.
#' @param plot Logical value.
#'
#' @return splits, typemarker, subsetter
#' @export
targeted_split <- function(
  x,
  typemarker,
  min_height = 0.1,
  min_score = 0.1,
  min_val_cutoff = NULL,
  max_val_cutoff = NULL,
  plot = TRUE
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
                                   min_score, min_height)
      found_valley <- any(!is.na(proposals[1, ]))

      if (found_valley) {
        found_boundary <- NA
        rect_col <- "green"
        score_name <- "depth"
      } else {
        proposals <- propose_boundaries(x, g, var_num, subsetter, progress)
        found_boundary <- any(!is.na(proposals[1, ]))
        rect_col <- "lightblue"
        score_name <- "bic"
      }

      # find the proposed valley with the highest score
      p_choice <- which.max(proposals[2, ])
      # put the valley and its score into the splits and scores matrices
      splits[g, p_choice] <- proposals[1, p_choice]
      scores[g, p_choice] <- proposals[2, p_choice]

      # if all of the current round's proposals are NA, use split_gmm,
      # otherwise, choose the proposed split with the highest score
      if (found_valley || found_boundary) {
        progress[g, p_choice] <- TRUE

        scale01_gp <- scale01(x[subsetter[, g], p_choice])

        # find which events have values less than the chosen split
        less_gp <- x[, p_choice] < splits[g, p_choice]
        is_neg_gp <- typemarker[g, p_choice] == -1
        subsetter[, g] <- subsetter[, g] & ((is_neg_gp & less_gp) | (!is_neg_gp & !less_gp))

        trans_split_gp <- scale01(splits[g, p_choice],
                                  scale01_gp$min, scale01_gp$max)$y
        xleft <- ifelse(is_neg_gp, 0, trans_split_gp)
        xright <- ifelse(is_neg_gp, trans_split_gp, 1)
        plot_targeted_split(stats::density(scale01_gp$y), g, p_choice,
                            scores[g, p_choice], typemarker,
                            xleft, xright, rect_col, trans_split_gp, score_name)
      } else {
        paused[g, !is.na(progress[g, ]) & !progress[g, ]] <- TRUE
        trans_split_gp <- xleft <- xright <- rect_col <- NA
        score_name <- "score"
        for (p in which(!is.na(progress[g, ]) & !progress[g, ])) {
          scale01_gp <- scale01(x[subsetter[, g], p])

          plot_targeted_split(stats::density(scale01_gp$y), g, p,
                              scores[g, p], typemarker,
                              xleft, xright, rect_col, trans_split_gp, score_name)
        }
      }
    }
  }

  if (path_num > 1) {
    equal_subsets <- matrix(nrow = path_num, ncol = path_num)
    is_a_duplicate <- rep(FALSE, path_num)
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

  for (g in 1:path_num) {
    for (p in 1:var_num) {
      if (is.na(splits[g, p])) {
        typemarker[g, p] <- 0
      }
    }
  }
  is_a_duplicate <- check_duplicates(subsetter, path_num)
  subsetter <- subsetter[, !is_a_duplicate]
  progress <- progress[!is_a_duplicate, ]
  splits <- splits[!is_a_duplicate, ]
  scores <- scores[!is_a_duplicate, ]
  paused <- paused[!is_a_duplicate, ]
  typemarker <- typemarker[!is_a_duplicate, ]
  path_num <- path_num - sum(is_a_duplicate)

  return(list(splits = splits,
              typemarker = typemarker,
              subsetter = subsetter))
}

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
      above_cutoff[, j] <- x[, j] < max_val_cutoff[j]
    }
  }
  # inside_cutoffs is a matrix of the same dimensions as the data
  # an observation-marker pair is TRUE if it is inside the two cutoffs
  inside_cutoffs <- !below_cutoff & !above_cutoff

  return(inside_cutoffs)
}

plot_targeted_split <- function(dens_gp, g, p, score, typemarker,
                                xleft, xright, rect_col, trans_split_gp,
                                score_name) {
  plot(dens_gp,
       main = paste0("g = ", g, ", p = ", p,
                     ", ", score_name, " = ", round(score, 3)),
       sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p]),
       panel.first = graphics::rect(xleft, 0, xright, max(dens_gp$y),
                                    col = rect_col, border =  NA))
  graphics::abline(v = trans_split_gp)
}

propose_valleys <- function(x, g, var_num, subsetter, progress, min_score, min_height) {
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

unscale01 <- function(x, unscaled_min, unscaled_max) {
  return(x * (unscaled_max - unscaled_min) + unscaled_min)
}

check_duplicates <- function(subsetter, path_num) {

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
