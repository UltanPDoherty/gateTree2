targeted_splits <- function(
  x,
  plusminus_table,
  min_depth = 1,
  min_height = 1
) {
  to_be_split <- apply(plusminus_table, 2, function(column) any(column != 0))

  obs_num <- nrow(x)
  var_num <- ncol(x)

  valleys <- boundaries <- splits <- rep(NA, var_num)
  types <- rep("unsplit", var_num)
  for (p in which(to_be_split)){
    scale01_p <- scale01(x[, p])
    dens01_p <- stats::density(scale01_p$y)

    valleys[p] <- find_valley(
      dens01_p,
      depth = TRUE,
      min_depth = min_depth,
      min_height = min_height
    )[1]

    valleys[p] <- unscale01(valleys[p], scale01_p$min, scale01_p$max)
    boundaries[p] <- find_boundary(x[, p])[1]

    splits[p] <- ifelse(!is.na(valleys[p]),
                        valleys[p], boundaries[p])

    if (!is.na(valleys[p])) {
      types[p] <- "valley"
    } else if (!is.na(boundaries[p])) {
      types[p] <- "boundary"
    } else {
      types[p] <- "failed"
    }
  }

  pop_num <- nrow(plusminus_table)
  pop_names <- rownames(plusminus_table)

  # create a +1/-1 matrix the same size as the data flowFrame
  is_positive <- t(apply(x, 1, function(x_row) x_row > splits))
  obs_plusminus <- (2 * is_positive) - 1

  subsetter <- matrix(NA, nrow = obs_num, ncol = pop_num,
                      dimnames = list(NULL, pop_names))
  nonzero_table <- plusminus_table != 0
  pop_signs <- list()
  for (j in seq_len(pop_num)){
    pop_signs[[j]] <- plusminus_table[j, nonzero_table[j, ]]
    for (i in seq_len(obs_num)) {
      subsetter[i, j] <- all(
        obs_plusminus[i, nonzero_table[j, ]] == pop_signs[[j]]
      )
    }
  }

  assigned <- apply(subsetter, 1, any)
  labels <- rep(0, obs_num)
  labels[assigned] <- apply(subsetter[assigned, ], 1, which)

  return(list(splits = splits,
              valleys = valleys,
              boundaries = boundaries,
              types = types,
              subsetter = subsetter,
              labels = labels))
}
