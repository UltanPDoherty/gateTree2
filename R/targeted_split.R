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
){

  G <- nrow(typemarker)
  P <- ncol(typemarker)
  N <- nrow(x)

  splits <- scores <- array(dim = dim(typemarker))

  subsetter <- matrix(TRUE, nrow = N, ncol = G)
  paused <- array(FALSE, dim = dim(typemarker))
  progress <- typemarker == 0

  inside_cutoffs <- find_inside_cutoffs(x, min_val_cutoff, max_val_cutoff)

  # loop over the gating pathways
  for (g in 1:G){
    # exclude observations that are outside the cutoffs for any marker used in
    # this gating pathway
    subsetter[, g] <- apply(inside_cutoffs[, typemarker[g, ] != 0], 1, all)

    graphics::par(mfrow = c(2, 2))
    # continue inner loop until this pathway's row of the progress & paused
    # matrices do not contain any FALSE values.
    # that is, move onto the next pathway only when every required variable for
    # this pathway is TRUE
    while (any(!progress[g, ] & !paused[g, ], na.rm = TRUE)) {

      proposals <- propose_splits(x, g, P, subsetter, progress, min_score, min_height)

      # if all of the current round's proposals are NA, use split_gmm,
      # otherwise, choose the proposed split with the highest score
      if (all(is.na(proposals[1, ]))) {
        less_gp <- x[, p_choice] < splits[g, p_choice]
        is_neg_gp <- typemarker[g, p_choice] == -1
        xleft <- ifelse(is_neg_gp, 0, trans_split_gp)
        xright <- ifelse(is_neg_gp, trans_split_gp, 1)
        subsetter[, g] <- subsetter[, g] & ((is_neg_gp & less_gp) | (!is_neg_gp & !less_gp))
        paused[g, !is.na(progress[g, ]) & !progress[g, ]] <- TRUE
        for (p in which(!is.na(progress[g, ]) & !progress[g, ])) {
          x_gp <- x[subsetter[, g], p]
          gmm_out <- split_gmm(x_gp)

          dens01_gp <- dens01(x_gp)
          if (gmm_out$bic_two < gmm_out$bic_one) {
            splits[g, p] <- gmm_out$split
            scores[g, p] <- -Inf
            less_gp <- x[, p] < splits[g, p]
            trans_split_gp <- (splits[g, p] - dens01_gp$min) / (dens01_gp$max - dens01_gp$min)

            is_neg_gp <- typemarker[g, p] == -1
            xleft <- ifelse(is_neg_gp, 0, trans_split_gp)
            xright <- ifelse(is_neg_gp, trans_split_gp, 1)
            subsetter[, g] <- subsetter[, g] & ((is_neg_gp & less_gp) | (!is_neg_gp & !less_gp))

            rect_col <- "lightblue"
          } else {
            trans_split_gp <- rect_col <- xleft <- xright <- NA
          }
          plot_targeted_split(dens01_gp$dens, g, p, scores[g, p], typemarker,
                              xleft, xright, rect_col, trans_split_gp)
        }
      } else {
        p_choice <- which.max(proposals[2, ])
        splits[g, p_choice] <- proposals[1, p_choice]
        scores[g, p_choice] <- proposals[2, p_choice]
        progress[g, p_choice] <- TRUE

        dens01_gp <- dens01(x[subsetter[, g], p_choice])

        trans_split_gp <- (splits[g, p_choice] - dens01_gp$min) / (dens01_gp$max - dens01_gp$min)

        rect_col <- "green"
        plot_targeted_split(dens01_gp$dens, g, p_choice, scores[g, p_choice], typemarker,
                            xleft, xright, rect_col, trans_split_gp)
      }
    }
  }

  if (G > 1) {
    equal_subsets <- matrix(nrow = G, ncol = G)
    is_a_duplicate <- rep(FALSE, G)
    for (g in 1:(G-1)) {
      for (h in (g + 1):G) {
        equal_subsets[g, h] <- all(subsetter[, g] == subsetter[, h])
        if (equal_subsets[g, h]) {
          print(paste0("Failed to distinguish between populations ",
                       g, " & ", h, "."))
          is_a_duplicate[h] <- TRUE
        }
      }
    }
    if (any(is_a_duplicate)) {
      to_be_deleted <- which(is_a_duplicate)

      subsetter <- subsetter[, -to_be_deleted]
      progress <- progress[-to_be_deleted, ]
      splits <- splits[-to_be_deleted, ]
      scores <- scores[-to_be_deleted, ]
      paused <- paused[-to_be_deleted, ]
      typemarker <- typemarker[-to_be_deleted, ]
      G <- G - sum(is_a_duplicate)
    }
  }

  for (g in 1:G) {
    for (p in 1:P) {
      if (is.na(splits[g, p])) {
        typemarker[g, p] <- 0
      }
    }
  }

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
    for (j in 1:ncol(x)) {
      below_cutoff[, j] <- x[, j] < min_val_cutoff[j]
    }
  }
  if (is.null(max_val_cutoff)) {
    above_cutoff <- array(FALSE, dim = dim(x))
  } else {
    above_cutoff <- array(dim = dim(x))
    for (j in 1:ncol(x)) {
      above_cutoff[, j] <- x[, j] < max_val_cutoff[j]
    }
  }
  # inside_cutoffs is a matrix of the same dimensions as the data
  # an observation-marker pair is TRUE if it is inside the two cutoffs
  inside_cutoffs <- !below_cutoff & !above_cutoff

  return(inside_cutoffs)
}

plot_targeted_split <- function(dens_gp, g, p, score, typemarker,
                                xleft, xright, rect_col, trans_split_gp) {
  plot(dens_gp,
       main = paste0("g = ", g, ", p = ", p,
                     ", score = ", round(score, 3)),
       sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p]),
       panel.first = graphics::rect(xleft, 0, xright, max(dens_gp$y),
                                    col = rect_col, border =  NA))
  graphics::abline(v = trans_split_gp)
}

propose_splits <- function(x, g, P, subsetter, progress, min_score, min_height){
  proposals <- matrix(nrow = 2, ncol = P)

  # loop over all variables to propose splits
  for (p in 1:P){
    # find variables that are not NA or TRUE in the progress matrix
    if (!is.na(progress[g, p]) & !progress[g, p]){
      # 0-1 scale this variable for the pathway's current subset
      scale_gp <- dens01(x[subsetter[, g], p])

      proposals[, p] <- find_valley(
        scale_gp$dens,
        score = TRUE,
        min_score = min_score,
        min_height = min_height)

      proposals[1, p] <- proposals[1, p] * (scale_gp$max - scale_gp$min) + scale_gp$min
    }
  }

  return(proposals)
}

dens01 <- function(x){
  min_x <- min(x)
  max_x <- max(x)
  dens_x <- stats::density((x - min_x) / (max_x - min_x))

  return(list(min = min_x, max = max_x, dens = dens_x))
}
