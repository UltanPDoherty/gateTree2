#' @title User-informed clustering decision tree.
#'
#' @description
#' Construct a semi-supervised decision tree to identify user-described
#' clusters.
#'
#' @param samples Dataset in `matrix` or `data.frame` form.
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
#'   list(iris[, -5]),
#'   iris_plusminus,
#'   min_depth = 10,
#'   min_diff = 0.01
#' )
gatetree <- function(
    samples,
    plusminus,
    min_depth = 50,
    min_diff = 0.01,
    use_gmm = TRUE,
    min_cutoffs = NULL,
    max_cutoffs = NULL,
    verbose = TRUE) {
  pop_num <- nrow(plusminus)
  var_num <- ncol(plusminus)
  samp_num <- length(samples)

  if (is.null(rownames(plusminus))) {
    rownames(plusminus) <- paste0("pop", seq_len(pop_num))
  }
  
  if (is.null(min_cutoffs)) {
    min_cutoffs <- rep(-Inf, var_num)
  }
  if (is.null(max_cutoffs)) {
    max_cutoffs <- rep(Inf, var_num)
  }

  pop_list <- list()
  for (p in seq_len(pop_num)) {
    pop_list[[p]] <- list(
      "subsetter" = lapply(samples, \(x) as.matrix(rep(TRUE, nrow(x)))),
      "pm_future" = plusminus[p, ],
      "pm_previous" = rep(0, var_num),
      "other_pops" = plusminus[-p, , drop = FALSE],
      "splits" = rep(list(rep(NA, var_num)), samp_num),
      "depths" = rep(list(rep(NA, var_num)), samp_num),
      "diffs" = rep(list(rep(NA, var_num)), samp_num),
      "order" = rep(NA, var_num),
      "method" = rep(NA, var_num),
      "min_cutoffs" = min_cutoffs,
      "max_cutoffs" = max_cutoffs,
      "terminated" = FALSE
    )

    names(pop_list)[p] <- rownames(plusminus)[p]

    pop_list[[p]] <- recursive_gatetree(
      pop_list[[p]], samples,
      min_depth = min_depth, min_diff = min_diff,
      use_gmm = use_gmm
    )

    if (verbose) {
      message("Population ", p, " complete.")
    }
  }

  pop_list
}

recursive_gatetree <- function(pop, samples, min_depth, min_diff, use_gmm) {
  if (pop$terminated) {
    return(pop)
  }

  var_num <- length(pop$pm_future)
  samp_num <- length(samples)
  split_num <- sum(pop$pm_previous != 0) + 1

  valleys <- rep(list(matrix(nrow = var_num, ncol = 2)), samp_num)
  for (s in seq_len(samp_num)) {
    for (v in seq_len(var_num)) {
      if (pop$pm_future[v] != 0 && all(pop$other_pops[, v] != 0)) {
        x <- samples[[s]][pop$subsetter[[s]][, split_num], v]
        x <- x[x > pop$min_cutoffs[v]]
        x <- x[x < pop$max_cutoffs[v]]
        # valleys[[s]][v, ] <- find_valley(stats::density(x), min_depth)
        valleys[[s]][v, ] <- find_valley2(x)
      } else {
        valleys[[s]][v, ] <- c(NA, NA)
      }
    }
  }

  samp_depths <- vapply(valleys, \(x) x[, 2], double(var_num))
  samp_depth_checks <- samp_depths > min_depth

  if (any(samp_depth_checks, na.rm = TRUE)) {
    samp_depths_na <- samp_depths
    samp_depths_na[!samp_depth_checks] <- NA
    samp_depth_means <- apply(samp_depths_na, 1, mean, na.rm = TRUE)
    samp_count_depth_checks <- apply(samp_depth_checks, 1, sum, na.rm = TRUE)
    max_depth_mean <- max(samp_depth_means[which.max(samp_count_depth_checks)])
    var_choice <- which(samp_depth_means == max_depth_mean)

    samp_valleys <- vapply(valleys, \(x) x[, 1], double(var_num))
    samp_valleys_na <- samp_valleys
    samp_valleys_na[!samp_depth_checks] <- NA
    samp_valley_means <- vapply(
      seq_len(var_num),
      FUN.VALUE = double(1L),
      \(x) {
        stats::weighted.mean(
          samp_valleys_na[x, ], samp_depths_na[x, ],
          na.rm = TRUE
        )
      }
    )

    valley_choices <- samp_valleys_na[var_choice, ]
    valley_choices[is.na(valley_choices)] <- samp_valley_means[var_choice]

    for (s in seq_len(samp_num)) {
      if (pop$pm_future[var_choice] == +1) {
        x <- samples[[s]][pop$subsetter[[s]][, split_num], var_choice]
        temp_subsetter <- x >= valley_choices[s]
        pop$pm_previous[var_choice] <- +1
      } else if (pop$pm_future[var_choice] == -1) {
        x <- samples[[s]][pop$subsetter[[s]][, split_num], var_choice]
        temp_subsetter <- x < valley_choices[s]
        pop$pm_previous[var_choice] <- -1
      } else {
        stop("pop$pm_future[var_choice] should not be 0.")
      }

      pop$splits[[s]][var_choice] <- valley_choices[s]
      pop$depths[[s]][var_choice] <- samp_depths_na[var_choice, s]

      pop$subsetter[[s]] <-
        cbind(pop$subsetter[[s]], pop$subsetter[[s]][, split_num])
      pop$subsetter[[s]][pop$subsetter[[s]][, split_num], split_num + 1] <-
        temp_subsetter

      if (s == samp_num) {
        pop$order[var_choice] <- split_num
        pop$method[var_choice] <- "valley"

        same_path <- pop$other_pops[, var_choice] == pop$pm_future[var_choice]
        pop$other_pops <- pop$other_pops[same_path, , drop = FALSE]
        pop$pm_future[var_choice] <- 0
        pop$other_pops[, var_choice] <- 0
      }
    }
  } else if (use_gmm) {
    boundaries <- rep(list(matrix(nrow = var_num, ncol = 2)), samp_num)
    for (s in seq_len(samp_num)) {
      for (v in seq_len(var_num)) {
        if (pop$pm_future[v] != 0 && all(pop$other_pops[, v] != 0)) {
          x <- samples[[s]][pop$subsetter[[s]][, split_num], v]
          x <- x[x > pop$min_cutoffs[v]]
          x <- x[x < pop$max_cutoffs[v]]
          boundaries[[s]][v, ] <- find_boundary(x, TRUE)
        } else {
          boundaries[[s]][v, ] <- c(NA, NA)
        }
      }
    }

    samp_diffs <- vapply(boundaries, \(x) x[, 2], double(var_num))
    samp_diff_checks <- samp_diffs > min_diff

    if (any(samp_diff_checks, na.rm = TRUE)) {
      samp_diffs_na <- samp_diffs
      samp_diffs_na[!samp_diff_checks] <- NA
      samp_diff_means <- apply(samp_diffs_na, 1, mean, na.rm = TRUE)
      samp_count_diff_checks <- apply(samp_diff_checks, 1, sum, na.rm = TRUE)
      max_diff_mean <- max(samp_diff_means[which.max(samp_count_diff_checks)])
      var_choice <- which(samp_diff_means == max_diff_mean)

      samp_boundaries <- vapply(boundaries, \(x) x[, 1], double(var_num))
      samp_boundaries_na <- samp_boundaries
      samp_boundaries_na[!samp_diff_checks] <- NA
      samp_boundary_means <- vapply(
        seq_len(var_num),
        FUN.VALUE = double(1L),
        \(x) {
          stats::weighted.mean(
            samp_boundaries_na[x, ], samp_diffs_na[x, ],
            na.rm = TRUE
          )
        }
      )

      boundary_choices <- samp_boundaries_na[var_choice, ]
      boundary_choices[is.na(boundary_choices)] <-
        samp_boundary_means[var_choice]

      for (s in seq_len(samp_num)) {
        if (pop$pm_future[var_choice] == +1) {
          x <- samples[[s]][pop$subsetter[[s]][, split_num], var_choice]
          temp_subsetter <- x >= boundary_choices[s]
          pop$pm_previous[var_choice] <- +1
        } else if (pop$pm_future[var_choice] == -1) {
          x <- samples[[s]][pop$subsetter[[s]][, split_num], var_choice]
          temp_subsetter <- x < boundary_choices[s]
          pop$pm_previous[var_choice] <- -1
        } else {
          stop("pop$pm_future[var_choice] should not be 0.")
        }

        pop$splits[[s]][var_choice] <- boundary_choices[s]
        pop$diffs[[s]][var_choice] <- samp_diffs_na[var_choice, s]

        pop$subsetter[[s]] <-
          cbind(pop$subsetter[[s]], pop$subsetter[[s]][, split_num])
        pop$subsetter[[s]][pop$subsetter[[s]][, split_num], split_num + 1] <-
          temp_subsetter

        if (s == samp_num) {
          pop$order[var_choice] <- split_num
          pop$method[var_choice] <- "boundary"

          same_path <- pop$other_pops[, var_choice] == pop$pm_future[var_choice]
          pop$other_pops <- pop$other_pops[same_path, , drop = FALSE]
          pop$pm_future[var_choice] <- 0
          pop$other_pops[, var_choice] <- 0
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
    pop, samples,
    min_depth = min_depth, min_diff = min_diff, use_gmm = use_gmm
  )
}
