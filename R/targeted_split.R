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

  for (g in 1:G){
    subsetter[, g] <- apply(inside_cutoffs[, typemarker[g, ] != 0], 1, all)

    graphics::par(mfrow = c(2, 2))
    while (any(!progress[g, ] & !paused[g, ], na.rm = TRUE)) {
      proposals <- matrix(nrow = 2, ncol = P)
      for (p in 1:P){
        if (!is.na(progress[g, p]) & !progress[g, p]){
          x_gp <- x[subsetter[, g], p]
          min_gp <- min(x_gp)
          max_gp <- max(x_gp)
          dens_gp <- stats::density((x_gp - min_gp) / (max_gp - min_gp))

          proposals[, p] <- find_valley(
            dens_gp,
            score = TRUE,
            min_score = min_score,
            min_height = min_height)
          proposals[1, p] <- min_gp + (max_gp - min_gp) * proposals[1, p]
        }
      }

      if (all(is.na(proposals[1, ]))) {
        paused[g, !is.na(progress[g, ]) & !progress[g, ]] <- TRUE

        for (p in which(!is.na(progress[g, ]) & !progress[g, ])) {
          x_gp <- x[subsetter[, g], p]
          min_gp <- min(x_gp)
          max_gp <- max(x_gp)
          dens_gp <- stats::density((x_gp - min_gp) / (max_gp - min_gp))

          gmm_out <- split_gmm(x_gp)

          if (gmm_out$bic_two < gmm_out$bic_one) {
            splits[g, p] <- gmm_out$split
            scores[g, p] <- -Inf
            trans_split_gp <- (splits[g, p] - min_gp) / (max_gp - min_gp)
            if (typemarker[g, p] == -1) {
              plot(dens_gp,
                   main = paste0("g = ", g, ", p = ", p),
                   sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p_choice]),
                   panel.first = graphics::rect(0, 0, trans_split_gp, max(dens_gp$y),
                                                col = "lightblue", border =  NA))
              subsetter[, g] <- subsetter[, g] & x[, p] < splits[g, p]
            } else if (typemarker[g, p] == +1) {
              plot(dens_gp,
                   main = paste0("g = ", g, ", p = ", p),
                   sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p_choice]),
                   panel.first = graphics::rect(trans_split_gp, 0, 1, max(dens_gp$y),
                                                col = "lightblue", border =  NA))
              subsetter[, g] <- subsetter[, g] & x[, p] > splits[g, p]
            }
            graphics::abline(v = trans_split_gp)
          } else {
            plot(dens_gp,
                 main = paste0("g = ", g, ", p = ", p),
                 sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p]))
          }
        }
        next
      } else {
        p_choice <- which.max(proposals[2, ])
        splits[g, p_choice] <- proposals[1, p_choice]
        scores[g, p_choice] <- proposals[2, p_choice]
        progress[g, p_choice] <- TRUE

        x_gp <- x[subsetter[, g], p_choice]
        min_gp <- min(x_gp)
        max_gp <- max(x_gp)
        dens_gp <- stats::density((x_gp - min_gp) / (max_gp - min_gp))
        trans_split_gp <- (splits[g, p_choice] - min_gp) / (max_gp - min_gp)
        if (typemarker[g, p_choice] == -1) {
          plot(dens_gp,
               main = paste0("g = ", g, ", p = ", p_choice,
                             ", score = ", round(scores[g, p_choice], 3)),
               sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p_choice]),
               panel.first = graphics::rect(0, 0, trans_split_gp, max(dens_gp$y),
                                            col = "green", border =  NA))
          subsetter[, g] <- subsetter[, g] & x[, p_choice] < splits[g, p_choice]
        } else if (typemarker[g, p_choice] == +1) {
          plot(dens_gp,
               main = paste0("g = ", g, ", p = ", p_choice,
                             ", score = ", round(scores[g, p_choice], 3)),
               sub = paste0(rownames(typemarker)[g], ", ", colnames(typemarker)[p_choice]),
               panel.first = graphics::rect(trans_split_gp, 0, 1, max(dens_gp$y),
                                            col = "green", border =  NA))
          subsetter[, g] <- subsetter[, g] & x[, p_choice] > splits[g, p_choice]
        }
        graphics::abline(v = trans_split_gp)
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
