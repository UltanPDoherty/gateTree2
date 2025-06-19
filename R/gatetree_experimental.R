#' @title User-informed clustering decision tree.
#'
#' @description
#' Construct a semi-supervised decision tree to identify user-described
#' clusters.
#'
#' @param matrices Dataset in `matrix` or `data.frame` form.
#' @param plusminus
#' Table indicating whether each population (row) is positive (+1),
#' negative (-1), or neutral / unknown (0) for each variable (column).
#' @param min_depth
#' Minimum depth, as a percentage of the height of the global density maximum,
#' for a split to be returned by [find_valley].
#' @param min_diff
#' Minimum value of difference between one-component and two-component BIC
#' divided by 2 log(obs_num).
#' @param use_gmm Logical value.
#' @param min_kde_size
#' Minimum number of events required to search for a KDE valley.
#' @param min_gmm_size
#' Minimum number of events required to search for a GMM boundary.
#' @param min_cutoffs Minimum values for observations used when finding splits.
#' @param max_cutoffs Maximum values for observations used when finding splits.
#' @param seed Random seed for GMM fitting.
#' @param verbose Logical value.
#'
#' @return A `list` object:
#'
#' @export
#'
#' @examples
#'
#' iris_plusminus <- rbind(
#'   "setosa" = c(0, 0, -1, 0),
#'   "versicolor" = c(0, 0, +1, -1),
#'   "virginica" = c(0, 0, +1, +1)
#' )
#' colnames(iris_plusminus) <- colnames(iris_plusminus)[1:4]
#'
#' iris_gatetree <- gatetree(
#'   list(list(iris[, -5])),
#'   iris_plusminus,
#'   min_depth = 10,
#'   min_diff = 0.01,
#'   min_gmm_size = 50
#' )
gatetree2 <- function(
    matrices,
    plusminus,
    min_depth = 0.05,
    min_diff = 0.05,
    use_gmm = TRUE,
    min_kde_size = 100,
    min_gmm_size = 100,
    min_cutoffs = NULL,
    max_cutoffs = NULL,
    seed = NULL,
    verbose = TRUE) {
  pop_num <- nrow(plusminus)
  var_num <- ncol(plusminus)
  samp_num <- vapply(matrices, length, integer(1L))

  if (is.null(rownames(plusminus))) {
    rownames(plusminus) <- paste0("pop", seq_len(pop_num))
  }
  if (is.null(colnames(plusminus))) {
    colnames(plusminus) <- paste0("var", seq_len(var_num))
  }

  if (is.null(seed)) {
    seed <- sample(1:1e6, size = 1)
  }

  if (is.null(min_cutoffs)) {
    min_cutoffs <- rep(-Inf, var_num)
  }
  if (is.null(max_cutoffs)) {
    max_cutoffs <- rep(Inf, var_num)
  }
  if (is.null(names(min_cutoffs))) {
    names(min_cutoffs) <- colnames(plusminus)
  }
  if (is.null(names(max_cutoffs))) {
    names(max_cutoffs) <- colnames(plusminus)
  }

  this_call <- call(
    "ombc_gmm",
    "matrices" = substitute(matrices),
    "plusminus" = plusminus,
    "min_depth" = min_depth, "min_diff" = min_diff,
    "use_gmm" = use_gmm,
    "min_kde_size" = min_kde_size, "min_gmm_size" = min_gmm_size,
    "min_cutoffs" = min_cutoffs, "max_cutoffs" = max_cutoffs, "seed" = seed,
    "verbose" = verbose
  )

  var_named_nas <- rep(NA, var_num)
  var_named_0s <- rep(0, var_num)
  names(var_named_nas) <- names(var_named_0s) <- colnames(plusminus)

  batch_num <- length(matrices)
  batch_matrix_list <- list()
  for (b in seq_len(batch_num)) {
    batch_matrix_list[[b]] <- matrix(nrow = samp_num[b], ncol = var_num)
    rownames(batch_matrix_list[[b]]) <- names(matrices[[b]])
    colnames(batch_matrix_list[[b]]) <- colnames(plusminus)
  }
  names(batch_matrix_list) <- names(matrices)
  
  pop_list <- list()
  for (p in seq_len(pop_num)) {
    pop_list[[p]] <- list(
      "subsetter" = lapply(
        matrices, \(x) lapply(x, \(y) as.matrix(rep(TRUE, nrow(y))))
      ),
      "splits" = batch_matrix_list,
      "depths" = batch_matrix_list,
      "diffs" = batch_matrix_list,
      "pm_future" = plusminus[p, ],
      "pm_previous" = var_named_0s,
      "other_pops" = plusminus[-p, , drop = FALSE],
      "order" = var_named_nas,
      "method" = batch_matrix_list,
      "min_cutoffs" = min_cutoffs,
      "max_cutoffs" = max_cutoffs,
      "terminated" = FALSE
    )
    names(pop_list)[p] <- rownames(plusminus)[p]

    pop_list[[p]] <- recursive_gatetree2(
      pop_list[[p]], matrices,
      min_depth = min_depth, min_diff = min_diff,
      use_gmm = use_gmm,
      min_kde_size = min_kde_size, min_gmm_size = min_gmm_size,
      seed = seed
    )

    if (verbose) {
      message("Population ", p, ", ", names(pop_list)[p], ", complete.")
    }
  }

  list("output" = pop_list, "call" = this_call)
}


recursive_gatetree2 <- function(
    pop, matrices, min_depth, min_diff, use_gmm,
    min_kde_size, min_gmm_size, seed) {
  if (pop$terminated) {
    return(pop)
  }

  var_num <- length(pop$pm_future)
  batch_num <- length(matrices)
  samp_num <- vapply(matrices, length, integer(1L))
  split_num <- sum(pop$pm_previous != 0) + 1

  splittable_vars <- logical(var_num)
  for (v in seq_len(var_num)) {
    splittable_vars[v] <- pop$pm_future[v] != 0 && all(pop$other_pops[, v] != 0)
  }

  valleys <- depths <- list()
  for (b in seq_len(batch_num)) {
    valleys[[b]] <- depths[[b]] <- matrix(nrow = samp_num[b], ncol = var_num)
    for (s in seq_len(samp_num[b])) {
      for (v in seq_len(var_num)) {
        if (splittable_vars[v]) {
          x <- matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], v]
          x <- x[x > pop$min_cutoffs[v]]
          x <- x[x < pop$max_cutoffs[v]]
          new_valley <- find_valley(x, min_kde_size = min_kde_size)
          if (is.na(new_valley[2]) || new_valley[2] < min_depth) {
            valleys[[b]][s, v] <- depths[[b]][s, v] <- NA
          } else {
            valleys[[b]][s, v] <- new_valley[1]
            depths[[b]][s, v] <- new_valley[2]
          }
        } else {
          valleys[[b]][s, v] <- depths[[b]][s, v] <- NA
        }
      }
    }
  }

  rbind_valleys <- Reduce(rbind, valleys)
  rbind_depths <- Reduce(rbind, depths)

  mean_depths <- apply(rbind_depths, 2, sum, na.rm = TRUE) / sum(samp_num)
  if (any(mean_depths > 0)) {
    var_choice <- which.max(mean_depths)
  } else {
    var_choice <- NA
  }

  boundary_needed <- matrix(nrow = var_num, ncol = batch_num)
  if (!use_gmm) {
    boundary_needed[] <- FALSE
  } else if (is.na(var_choice)) {
    boundary_needed[splittable_vars, ] <- TRUE
    boundary_needed[!splittable_vars, ] <- FALSE
  } else {
    for (v in seq_len(var_num)) {
      if (v == var_choice) {
        for (b in seq_len(batch_num)) {
          if (all(is.na(depths[[b]][, v]))) {
            boundary_needed[v, b] <- TRUE
          } else {
            boundary_needed[v, b] <- FALSE
          }
        }
      } else {
        boundary_needed[] <- FALSE
      }
    }
  }

  boundaries <- diffs <- list()
  for (b in seq_len(batch_num)) {
    boundaries[[b]] <- diffs[[b]] <- matrix(nrow = samp_num[b], ncol = var_num)
    for (s in seq_len(samp_num[b])) {
      for (v in seq_len(var_num)) {
        if (boundary_needed[v, b]) {
          x <- matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], v]
          x <- x[x > pop$min_cutoffs[v]]
          x <- x[x < pop$max_cutoffs[v]]
          set.seed(seed)
          new_boundary <- find_boundary(
            x,
            min_gmm_size = min_gmm_size
          )
          if (is.na(new_boundary[2]) || new_boundary[2] < min_diff) {
            boundaries[[b]][s, v] <- diffs[[b]][s, v] <- NA
          } else {
            boundaries[[b]][s, v] <- new_boundary[1]
            diffs[[b]][s, v] <- new_boundary[2]
          }
        } else {
          boundaries[[b]][s, v] <- diffs[[b]][s, v] <- NA
        }
      }
    }
  }

  rbind_boundaries <- Reduce(rbind, boundaries)
  rbind_diffs <- Reduce(rbind, diffs)

  if (is.na(var_choice)) {
    mean_diffs <- apply(rbind_diffs, 2, sum, na.rm = TRUE) / sum(samp_num)
    if (any(mean_diffs > 0)) {
      var_choice <- which.max(mean_diffs)
    } else {
      var_choice <- NA
    }
  }

  if (!is.na(var_choice)) {
    batch_valley_means <- batch_boundary_means <- double(batch_num)
    for (b in seq_len(batch_num)) {
      batch_valley_means[b] <- stats::weighted.mean(
        valleys[[b]][, var_choice],
        depths[[b]][, var_choice],
        na.rm = TRUE
      )
      batch_boundary_means[b] <- stats::weighted.mean(
        boundaries[[b]][, var_choice],
        diffs[[b]][, var_choice],
        na.rm = TRUE
      )
    }

    study_valley_mean <- stats::weighted.mean(
      rbind_valleys[, var_choice],
      rbind_depths[, var_choice],
      na.rm = TRUE
    )
    study_boundary_mean <- stats::weighted.mean(
      rbind_boundaries[, var_choice],
      rbind_diffs[, var_choice],
      na.rm = TRUE
    )

    splits <- mechanisms <- list()
    for (b in seq_len(batch_num)) {
      splits[[b]] <- double(samp_num[b])
      mechanisms[[b]] <- character(samp_num[b])
      for (s in seq_len(samp_num[b])) {
        if (!is.na(valleys[[b]][s, var_choice])) {
          splits[[b]][s] <- valleys[[b]][s, var_choice]
          mechanisms[[b]][s] <- "valley"
        } else if (any(!is.na(valleys[[b]][, var_choice]))) {
          splits[[b]][s] <- batch_valley_means[b]
          mechanisms[[b]][s] <- "batch_valley"
        } else if (!is.na(boundaries[[b]][s, var_choice])) {
          splits[[b]][s] <- boundaries[[b]][s, var_choice]
          mechanisms[[b]][s] <- "boundary"
        } else if (any(!is.na(boundaries[[b]][, var_choice]))) {
          splits[[b]][s] <- batch_boundary_means[b]
          mechanisms[[b]][s] <- "batch_boundary"
        } else if (any(!is.na(rbind_depths[, var_choice]))) {
          splits[[b]][s] <- study_valley_mean
          mechanisms[[b]][s] <- "study_valley"

          message(paste0(
            "Study-level valley imputation has been conducted. Variable: ",
            var_choice, ". Batch: ", b, "."
          ))
        } else if (any(!is.na(rbind_diffs[, var_choice]))) {
          splits[[b]][s] <- study_boundary_mean
          mechanisms[[b]][s] <- "study_boundary"

          message(paste0(
            "Study-level boundary imputation has been conducted. Variable: ",
            var_choice, ". Batch: ", b, "."
          ))
        } else {
          stop("Unexpected scenario.")
        }
      }
    }

    for (b in seq_len(batch_num)) {
      for (s in seq_len(samp_num[b])) {
        if (pop$pm_future[var_choice] == +1) {
          x <-
            matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], var_choice]
          temp_subsetter <- x >= splits[[b]][s]
          pop$pm_previous[var_choice] <- +1
        } else if (pop$pm_future[var_choice] == -1) {
          x <-
            matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], var_choice]
          temp_subsetter <- x < splits[[b]][s]
          pop$pm_previous[var_choice] <- -1
        } else {
          stop("pop$pm_future[var_choice] should not be 0.")
        }

        pop$splits[[b]][s, var_choice] <- splits[[b]][s]
        if (mechanisms[[b]][s] == "valley") {
          pop$depths[[b]][s, var_choice] <- depths[[b]][s, var_choice]
        } else if (mechanisms[[b]][s] == "boundary") {
          pop$diffs[[b]][s, var_choice] <- diffs[[b]][s, var_choice]
        }

        pop$subsetter[[b]][[s]] <-
          cbind(pop$subsetter[[b]][[s]], pop$subsetter[[b]][[s]][, split_num])
        pop$subsetter[[b]][[s]][
          pop$subsetter[[b]][[s]][, split_num], split_num + 1
        ] <- temp_subsetter

        pop$method[[b]][s, var_choice] <- mechanisms[[b]][s]

        if (b == batch_num && s == samp_num[b]) {
          pop$order[var_choice] <- split_num

          same_path <- pop$other_pops[, var_choice] == pop$pm_future[var_choice]
          pop$other_pops <- pop$other_pops[same_path, , drop = FALSE]
          pop$pm_future[var_choice] <- 0
          pop$other_pops[, var_choice] <- 0
        }
      }
    }
  } else {
    pop$terminated <- TRUE
  }

  if (all(pop$pm_future == 0)) {
    pop$terminated <- TRUE
  }

  recursive_gatetree2(
    pop, matrices,
    min_depth = min_depth, min_diff = min_diff,
    use_gmm = use_gmm, min_kde_size = min_kde_size, min_gmm_size = min_gmm_size,
    seed = seed
  )
}
