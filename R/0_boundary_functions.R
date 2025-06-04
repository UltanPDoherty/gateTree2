#' @title Find the optimal univariate GMM boundary.
#'
#' @description
#' Find the optimal two-component univariate GMM boundary and check if its
#' scaled BIC difference is greater than `min_scaled_bic_diff`.
#'
#' @inheritParams gatetree
#'
#' @return Vector consisting of the boundary and its scaled BIC difference. The
#' boundary will be NA if its scaled BIC difference is less than
#' `min_scaled_bic_diff`.
find_boundary <- function(x, min_scaled_bic_diff = 0) {
  obs_num <- length(x)
  sortx <- sort(x)

  leftx <- sortx[10]
  rightx <- sortx[obs_num - 11]
  xseq <- seq(leftx, rightx, length.out = 100)
  model_num <- length(xseq)

  prop <- double(2)
  mu <- double(2)
  sigma <- double(2)
  ll <- double(model_num)
  comp_pdf <- matrix(nrow = obs_num, ncol = 2)
  for (j in 1:model_num) {
    neg <- x <= xseq[j]
    prop[1] <- sum(neg) / obs_num
    prop[2] <- 1 - prop[1]
    mu[1] <- mean(x[neg])
    mu[2] <- mean(x[!neg])
    sigma[1] <- stats::sd(x[neg])
    sigma[2] <- stats::sd(x[!neg])

    comp_pdf[, 1] <- stats::dnorm(x, mu[1], sigma[1])
    comp_pdf[, 2] <- stats::dnorm(x, mu[2], sigma[2])

    ll[j] <- sum(log(comp_pdf %*% prop))
  }

  ll_max <- max(ll)
  ll_one <- sum(stats::dnorm(x, mean(x), stats::sd(x), log = TRUE))

  bic_one <- 2 * log(obs_num) - 2 * ll_one
  bic_two <- 4 * log(obs_num) - 2 * ll_max

  scaled_bic_diff <- (bic_one - bic_two) / (2 * log(obs_num))

  if (scaled_bic_diff > min_scaled_bic_diff) {
    boundary <- xseq[which.max(ll)]
  } else {
    boundary <- NA
  }

  return(c(boundary, scaled_bic_diff))
}

# ==============================================================================

#' @title Apply [find_boundary] across a set of variables.
#'
#' @description
#' Apply [find_boundary] to each variable in a set of splittable
#' variables for a subset of observations.
#'
#' @inheritParams propose_splits
#'
#' @return List:
#' * splits: each variable's boundary
#' * scores: the corresponding scaled BIC difference.
propose_boundaries <- function(
    x,
    min_scaled_bic_diff = 0,
    subsetter_g,
    splittable_vars_g = rep(TRUE, ncol(x))) {
  boundaries <- matrix(nrow = 2, ncol = ncol(x))

  # loop over all variables to propose splits
  for (p in which(splittable_vars_g)) {
    boundaries[, p] <- find_boundary(x[subsetter_g, p], min_scaled_bic_diff)
  }

  return(list("splits" = boundaries[1, ], "scores" = boundaries[2, ]))
}

# ==============================================================================

#' @title Find the optimal univariate GMM boundary.
#'
#' @description
#' Find the optimal two-component univariate GMM boundary and check if its
#' scaled BIC difference is greater than `min_scaled_bic_diff`.
#'
#' @inheritParams gatetree
#'
#' @return Vector consisting of the boundary and its scaled BIC difference. The
#' boundary will be NA if its scaled BIC difference is less than
#' `min_scaled_bic_diff`.
find_boundary2 <- function(x, min_scaled_bic_diff = 0) {
  obs_num <- length(x)

  gmm1 <- mclust::Mclust(x, G = 1)
  gmm2 <- mclust::Mclust(x, G = 2)

  dens1 <- mclust::dens(x, gmm1$modelName, gmm1$parameters)
  dens2 <- mclust::dens(x, gmm2$modelName, gmm2$parameters)

  noise1 <- dens1 < mclust::hypvol(x, TRUE)
  noise2 <- dens2 < mclust::hypvol(x, TRUE)

  gmm1n <- mclust::Mclust(x, G = 1, initialization = list(noise = noise1))
  gmm2n <- mclust::Mclust(x, G = 2, initialization = list(noise = noise2))

  bic <- c(gmm1$bic, gmm2$bic, gmm1n$bic, gmm2n$bic)
  model_id <- c("G1", "G2", "G1N", "G2N")
  bic_choice <- model_id[which.max(bic)]

  if (bic_choice == "G1") {
    boundary <- NA
    scaled_bic_diff <- 0
  } else if (bic_choice == "G2") {
    comp1 <- which.min(gmm2$parameters$mean)
    comp2 <- 3 - comp1

    upper_x <- x > gmm2$parameters$mean[comp1]
    lower_x <- x < gmm2$parameters$mean[comp2]

    max_comp1 <- max(x[lower_x & upper_x & gmm2$classification == comp1])
    min_comp2 <- max(x[lower_x & upper_x & gmm2$classification == comp2])

    boundary <- mean(c(max_comp1, min_comp2))

    scaled_bic_diff <- (bic[1] - bic[2]) / (2 * log(obs_num))
  } else if (bic_choice == "G1N") {
    upper_x <- x > gmm1n$parameters$mean
    lower_x <- x < gmm1n$parameters$mean

    upper_count <- sum(upper_x & gmm1n$classification == 0)
    lower_count <- sum(lower_x & gmm1n$classification == 0)

    if (upper_count > lower_count) {
      max_comp1 <- max(x[gmm2$classification == 1])
      min_noise <- min(x[upper_x & gmm2$classification == 0])

      boundary <- mean(c(max_comp1, min_noise))
    } else {
      min_comp1 <- min(x[gmm2$classification == 1])
      max_noise <- max(x[lower_x & gmm2$classification == 0])

      boundary <- mean(c(min_comp1, max_noise))
    }
    scaled_bic_diff <- (bic[1] - bic[3]) / (2 * log(obs_num))
  } else {
    comp1 <- which.min(gmm2n$parameters$mean)
    comp2 <- 3 - comp1

    upper_x <- x > gmm2n$parameters$mean[comp1]
    lower_x <- x < gmm2n$parameters$mean[comp2]

    max_comp1 <- max(x[lower_x & upper_x & gmm2n$classification == comp1])
    min_comp2 <- max(x[lower_x & upper_x & gmm2n$classification == comp2])

    boundary <- mean(c(max_comp1, min_comp2))

    scaled_bic_diff <- (bic[1] - bic[4]) / (2 * log(obs_num))
  }

  return(c(boundary, as.numeric(scaled_bic_diff)))
}

# ==============================================================================

#' @title Find the optimal univariate GMM boundary.
#'
#' @description
#' Find the optimal two-component univariate GMM boundary and check if its
#' scaled BIC difference is greater than `min_scaled_bic_diff`.
#'
#' @inheritParams gatetree
#' @param verbose message about `best_model`.
#'
#' @return Vector consisting of the boundary and its scaled BIC difference. The
#' boundary will be NA if its scaled BIC difference is less than
#' `min_scaled_bic_diff`.
#'
#' @import mclust
find_boundary3 <- function(x, min_scaled_bic_diff = 0, verbose = FALSE) {
  obs_num <- length(x)
  sortx <- sort(x)

  leftx <- sortx[10]
  rightx <- sortx[obs_num - 11]
  xseq <- seq(leftx, rightx, length.out = 100)
  model_num <- length(xseq)

  npar <- c(
    "gmm_gmm" = 4,
    "gmm_gmmn" = 5,
    "gmmn_gmm" = 5,
    "gmmn_gmmn" = 6,
    "gmm_unif" = 3,
    "unif_gmm" = 3,
    "gmmn_unif" = 4,
    "unif_gmmn" = 4,
    "unif_unif" = 2
  )

  best_model <- integer(model_num)
  best_bic <- double(model_num)
  for (j in 1:model_num) {
    neg <- x <= xseq[j]

    neg_gmm <- mclust::Mclust(x[neg], G = 1, modelNames = "V", verbose = FALSE)
    pos_gmm <- mclust::Mclust(x[!neg], G = 1, modelNames = "V", verbose = FALSE)

    neg_dens <- mclust::dens(x[neg], "V", neg_gmm$parameters)
    pos_dens <- mclust::dens(x[!neg], "V", pos_gmm$parameters)

    neg_const <- mclust::hypvol(x[neg], TRUE)
    pos_const <- mclust::hypvol(x[!neg], TRUE)

    neg_noise <- neg_dens < neg_const
    pos_noise <- pos_dens < pos_const

    if (sum(neg) - sum(neg_noise) > 5) {
      neg_gmmn <- mclust::Mclust(
        x[neg],
        G = 1, modelNames = "V",
        initialization = list(noise = neg_noise), verbose = FALSE
      )
      neg_gmmn_bic <- neg_gmmn$bic
    } else {
      neg_gmmn_bic <- -Inf
    }

    if (sum(!neg) - sum(pos_noise) > 5) {
      pos_gmmn <- mclust::Mclust(
        x[!neg],
        G = 1, modelNames = "V",
        initialization = list(noise = pos_noise), verbose = FALSE
      )
      pos_gmmn_bic <- pos_gmmn$bic
    } else {
      pos_gmmn_bic <- -Inf
    }

    neg_unif_ll <- sum(neg) * log(neg_const)
    pos_unif_ll <- sum(!neg) * log(pos_const)

    neg_unif_bic <- 2 * neg_unif_ll - log(sum(neg))
    pos_unif_bic <- 2 * pos_unif_ll - log(sum(!neg))

    bic <- c(
      "gmm_gmm" = neg_gmm$bic + pos_gmm$bic,
      "gmm_gmmn" = neg_gmm$bic + pos_gmmn_bic,
      "gmmn_gmm" = neg_gmmn_bic + pos_gmm$bic,
      "gmmn_gmmn" = neg_gmmn_bic + pos_gmmn_bic,
      "gmm_unif" = neg_gmm$bic + pos_unif_bic,
      "unif_gmm" = neg_unif_bic + pos_gmm$bic,
      "gmmn_unif" = neg_gmmn_bic + pos_unif_bic,
      "unif_gmmn" = neg_unif_bic + pos_gmmn_bic,
      "unif_unif" = neg_unif_bic + pos_unif_bic
    )

    best_bic[j] <- max(bic)
    best_model[j] <- which.max(bic)

    if (verbose) {
      message("*", appendLF = FALSE)
      if (j %% 10 == 0) {
        message("\tSplit No. ", j)
      }
    }
  }

  ll_one <- sum(stats::dnorm(x, mean(x), stats::sd(x), log = TRUE))
  bic_one <- 2 * ll_one - 2 * log(obs_num)

  choice <- which.max(best_bic)

  scaled_bic_diff <- (best_bic[choice] - bic_one) / (2 * log(obs_num))

  if (verbose) {
    message("Chosen model was ", names(npar)[best_model[choice]])
  }
  if (scaled_bic_diff > min_scaled_bic_diff) {
    boundary <- xseq[choice]
  } else {
    boundary <- NA
  }

  return(c(boundary, scaled_bic_diff))
}
