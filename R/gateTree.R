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
#'   min_diff = 0.01
#' )
gatetree <- function(
    matrices,
    plusminus,
    min_depth = 100,
    min_diff = 0.05,
    use_gmm = TRUE,
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

  if (is.null(seed)) {
    seed <- sample(1:1e6, size = 1)
  }

  if (is.null(min_cutoffs)) {
    min_cutoffs <- rep(-Inf, var_num)
  }
  if (is.null(max_cutoffs)) {
    max_cutoffs <- rep(Inf, var_num)
  }

  this_call <- call(
    "ombc_gmm",
    "matrices" = substitute(matrices),
    "plusminus" = plusminus,
    "min_depth" = min_depth, "min_diff" = min_diff, "use_gmm" = use_gmm,
    "min_cutoffs" = min_cutoffs, "max_cutoffs" = max_cutoffs, "seed" = seed,
    "verbose" = verbose
  )

  pop_list <- list()
  for (p in seq_len(pop_num)) {
    pop_list[[p]] <- list(
      "subsetter" = lapply(
        matrices, \(x) lapply(x, \(y) as.matrix(rep(TRUE, nrow(y))))
      ),
      "splits" = lapply(samp_num, \(x) rep(list(rep(NA, var_num)), x)),
      "depths" = lapply(samp_num, \(x) rep(list(rep(NA, var_num)), x)),
      "diffs" = lapply(samp_num, \(x) rep(list(rep(NA, var_num)), x)),
      "pm_future" = plusminus[p, ],
      "pm_previous" = rep(0, var_num),
      "other_pops" = plusminus[-p, , drop = FALSE],
      "order" = rep(NA, var_num),
      "method" = rep(NA, var_num),
      "min_cutoffs" = min_cutoffs,
      "max_cutoffs" = max_cutoffs,
      "terminated" = FALSE
    )
    names(pop_list)[p] <- rownames(plusminus)[p]

    pop_list[[p]] <- recursive_gatetree(
      pop_list[[p]], matrices,
      min_depth = min_depth, min_diff = min_diff,
      use_gmm = use_gmm, seed = seed
    )

    if (verbose) {
      message("Population ", p, " complete.")
    }
  }

  list("output" = pop_list, "call" = this_call)
}

recursive_gatetree <- function(
    pop, matrices, min_depth, min_diff, use_gmm, seed) {
  if (pop$terminated) {
    return(pop)
  }

  var_num <- length(pop$pm_future)
  batch_num <- length(matrices)
  samp_num <- vapply(matrices, length, integer(1L))
  split_num <- sum(pop$pm_previous != 0) + 1

  valleys <- list()
  samp_depths <- list()
  samp_depth_checks <- list()
  for (b in seq_len(batch_num)) {
    valleys[[b]] <- rep(list(matrix(nrow = var_num, ncol = 2)), samp_num[b])
    for (s in seq_len(samp_num[b])) {
      for (v in seq_len(var_num)) {
        if (pop$pm_future[v] != 0 && all(pop$other_pops[, v] != 0)) {
          x <- matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], v]
          x <- x[x > pop$min_cutoffs[v]]
          x <- x[x < pop$max_cutoffs[v]]
          valleys[[b]][[s]][v, ] <- find_valley(x)
        } else {
          valleys[[b]][[s]][v, ] <- c(NA, NA)
        }
      }
    }
    samp_depths[[b]] <- vapply(valleys[[b]], \(x) x[, 2], double(var_num))
    samp_depth_checks[[b]] <- samp_depths[[b]] > min_depth
  }

  bind_depth_checks <- Reduce(cbind, samp_depth_checks)

  if (any(bind_depth_checks, na.rm = TRUE)) {
    choices <- make_choices(
      valleys, samp_depths, samp_depth_checks, matrices, pop$subsetter
    )
    var_choice <- choices$var
    valley_choices <- choices$splits
    choice_depths <- choices$scores

    for (b in seq_len(batch_num)) {
      for (s in seq_len(samp_num[b])) {
        if (pop$pm_future[var_choice] == +1) {
          x <-
            matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], var_choice]
          temp_subsetter <- x >= valley_choices[[b]][s]
          pop$pm_previous[var_choice] <- +1
        } else if (pop$pm_future[var_choice] == -1) {
          x <-
            matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], var_choice]
          temp_subsetter <- x < valley_choices[[b]][s]
          pop$pm_previous[var_choice] <- -1
        } else {
          browser()
          stop("pop$pm_future[var_choice] should not be 0.")
        }

        pop$splits[[b]][[s]][var_choice] <- valley_choices[[b]][s]
        pop$depths[[b]][[s]][var_choice] <- choice_depths[[b]][s]

        pop$subsetter[[b]][[s]] <-
          cbind(pop$subsetter[[b]][[s]], pop$subsetter[[b]][[s]][, split_num])
        pop$subsetter[[b]][[s]][
          pop$subsetter[[b]][[s]][, split_num], split_num + 1
        ] <- temp_subsetter

        if (b == batch_num && s == samp_num[b]) {
          pop$order[var_choice] <- split_num
          pop$method[var_choice] <- "valley"

          same_path <- pop$other_pops[, var_choice] == pop$pm_future[var_choice]
          pop$other_pops <- pop$other_pops[same_path, , drop = FALSE]
          pop$pm_future[var_choice] <- 0
          pop$other_pops[, var_choice] <- 0
        }
      }
    }
  } else if (use_gmm) {
    boundaries <- list()
    samp_diffs <- list()
    samp_diff_checks <- list()
    for (b in seq_len(batch_num)) {
      boundaries[[b]] <-
        rep(list(matrix(nrow = var_num, ncol = 2)), samp_num[b])
      for (s in seq_len(samp_num[b])) {
        for (v in seq_len(var_num)) {
          if (pop$pm_future[v] != 0 && all(pop$other_pops[, v] != 0)) {
            x <- matrices[[b]][[s]][pop$subsetter[[b]][[s]][, split_num], v]
            x <- x[x > pop$min_cutoffs[v]]
            x <- x[x < pop$max_cutoffs[v]]
            boundaries[[b]][[s]][v, ] <- find_boundary(x)
          } else {
            boundaries[[b]][[s]][v, ] <- c(NA, NA)
          }
        }
      }
      samp_diffs[[b]] <- vapply(boundaries[[b]], \(x) x[, 2], double(var_num))
      samp_diff_checks[[b]] <- samp_diffs[[b]] > min_diff
    }

    bind_diff_checks <- Reduce(cbind, samp_diff_checks)

    if (any(bind_diff_checks, na.rm = TRUE)) {
      choices <- make_choices(
        boundaries, samp_diffs, samp_diff_checks, matrices, pop$subsetter
      )
      var_choice <- choices$var
      boundary_choices <- choices$splits
      choice_diffs <- choices$scores

      for (b in seq_len(batch_num)) {
        for (s in seq_len(samp_num[b])) {
          if (pop$pm_future[var_choice] == +1) {
            x <- matrices[[b]][[s]][
              pop$subsetter[[b]][[s]][, split_num],
              var_choice
            ]
            temp_subsetter <- x >= boundary_choices[[b]][s]
            pop$pm_previous[var_choice] <- +1
          } else if (pop$pm_future[var_choice] == -1) {
            x <- matrices[[b]][[s]][
              pop$subsetter[[b]][[s]][, split_num],
              var_choice
            ]
            temp_subsetter <- x < boundary_choices[[b]][s]
            pop$pm_previous[var_choice] <- -1
          } else {
            browser()
            stop("pop$pm_future[var_choice] should not be 0.")
          }

          pop$splits[[b]][[s]][var_choice] <- boundary_choices[[b]][s]
          pop$diffs[[b]][[s]][var_choice] <- choice_diffs[[b]][s]

          pop$subsetter[[b]][[s]] <-
            cbind(pop$subsetter[[b]][[s]], pop$subsetter[[b]][[s]][, split_num])
          pop$subsetter[[b]][[s]][
            pop$subsetter[[b]][[s]][, split_num], split_num + 1
          ] <- temp_subsetter

          if (b == batch_num && s == samp_num[b]) {
            pop$order[var_choice] <- split_num
            pop$method[var_choice] <- "boundary"

            same_path <-
              pop$other_pops[, var_choice] == pop$pm_future[var_choice]
            pop$other_pops <- pop$other_pops[same_path, , drop = FALSE]
            pop$pm_future[var_choice] <- 0
            pop$other_pops[, var_choice] <- 0
          }
        }
      }
    } else {
      pop$terminated <- TRUE
    }
  } else {
    pop$terminated <- TRUE
  }

  if (all(pop$pm_future == 0)) {
    pop$terminated <- TRUE
  }

  recursive_gatetree(
    pop, matrices,
    min_depth = min_depth, min_diff = min_diff, use_gmm = use_gmm, seed = seed
  )
}

make_choices <- function(splits, scores, checks, matrices, subsetter) {
  var_num <- nrow(scores[[1]])
  batch_num <- length(matrices)
  samp_num <- vapply(matrices, length, integer(1L))

  split_vals <- list()
  score_sums <- list()
  score_means <- list()
  for (b in seq_len(batch_num)) {
    scores[[b]][!checks[[b]]] <- NA
    split_vals[[b]] <- vapply(splits[[b]], \(x) x[, 1], double(var_num))
    split_vals[[b]][!checks[[b]]] <- NA

    score_sums[[b]] <- as.matrix(apply(scores[[b]], 1, sum, na.rm = TRUE))
    score_means[[b]] <- score_sums[[b]] / samp_num[b]
  }
  bind_score_means <- Reduce(cbind, score_means)
  mean_score_means <- apply(
    bind_score_means, 1, \(x) sum(x * samp_num / sum(samp_num))
  )
  var_choice <- which.max(mean_score_means)

  choice_scores <- lapply(scores, \(x) x[var_choice, ])

  comparable_splits_batch <- list()
  split_means <- double(batch_num)
  for (b in seq_len(batch_num)) {
    comparable_splits_batch[[b]] <- compare_splits(
      split_vals[[b]], scores[[b]], matrices[[b]], subsetter[[b]], var_choice
    )
    split_means[b] <- stats::weighted.mean(
      split_vals[[b]][var_choice, ] * comparable_splits_batch[[b]],
      scores[[b]][var_choice, ] * comparable_splits_batch[[b]],
      na.rm = TRUE
    )
  }

  comparable_splits_study <- compare_splits(
    Reduce(cbind, split_vals), Reduce(cbind, scores),
    Reduce(append, matrices), Reduce(append, subsetter),
    var_choice
  )
  split_mean_study <- stats::weighted.mean(
    Reduce(cbind, split_vals)[var_choice, ] * comparable_splits_study,
    Reduce(cbind, scores)[var_choice, ] * comparable_splits_study,
    na.rm = TRUE
  )
  split_means[is.na(split_means)] <- split_mean_study

  split_choices <- list()
  for (b in seq_len(batch_num)) {
    split_choices[[b]] <- split_vals[[b]][var_choice, ]
    split_choices[[b]][is.na(split_choices[[b]])] <- split_means[b]
  }

  list("var" = var_choice, "splits" = split_choices, "scores" = choice_scores)
}

compare_splits <- function(
    split_vals, scores, matrices, subsetter, var_choice) {
  samp_num <- ncol(split_vals)
  subset_num <- ncol(subsetter[[1]])

  balanced_accuracy <- double(samp_num)
  comparable_splits <- logical(samp_num)
  if (any(!is.na(scores[var_choice, ]))) {
    max_s <- which.max(scores[var_choice, ])

    for (s in seq_along(matrices)) {
      if (!is.na(scores[var_choice, s])) {
        obs_num <- sum(subsetter[[s]][, subset_num])

        original_neg <- sum(
          matrices[[s]][subsetter[[s]][, subset_num], var_choice] <
            split_vals[var_choice, s]
        )
        original_pos <- obs_num - original_neg

        maximum_neg <- sum(
          matrices[[s]][subsetter[[s]][, subset_num], var_choice] <
            split_vals[var_choice, max_s]
        )
        maximum_pos <- obs_num - maximum_neg

        tp <- min(original_pos, maximum_pos)
        tn <- min(original_neg, maximum_neg)
        fp <- max(maximum_pos - original_pos, 0)
        fn <- max(maximum_neg - original_neg, 0)

        fp + fn / obs_num
        balanced_accuracy[s] <- ((tp / (tp + fp)) + (tn / (tn + fn))) / 2
        comparable_splits[s] <- balanced_accuracy[s] >= 0.75
      }
    }
  }

  comparable_splits
}
