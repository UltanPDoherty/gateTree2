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
    "gmm_gmmn" = 6,
    "gmmn_gmm" = 6,
    "gmmn_gmmn" = 8,
    "gmm_unif" = 3,
    "unif_gmm" = 3,
    "gmmn_unif" = 5,
    "unif_gmmn" = 5,
    "unif_unif" = 2
  )

  best_model <- integer(model_num)
  best_bic <- double(model_num)
  for (j in 1:model_num) {
    neg <- x <= xseq[j]

    neg_prop <- sum(neg) / obs_num
    pos_prop <- 1 - neg_prop

    neg_dens <- stats::dnorm(x[neg], mean(x[neg]), sd(x[neg]))
    pos_dens <- stats::dnorm(x[!neg], mean(x[!neg]), sd(x[!neg]))

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

      neg_gmmn_dens <- stats::dnorm(
        x[neg], neg_gmmn$parameters$mean, sqrt(neg_gmmn$parameters$variance$sigmasq)
      )

      neg_gmmn_unif_prop <- sum(neg_gmmn$classification == 0) / sum(neg)
      neg_gmmn_norm_prop <- 1 - neg_gmmn_unif_prop
    } else {
      neg_gmmn_dens <- rep(0, sum(neg))

      neg_gmmn_unif_prop <- 0.5
      neg_gmmn_norm_prop <- 1 - neg_gmmn_unif_prop
    }

    if (sum(!neg) - sum(pos_noise) > 5) {
      pos_gmmn <- mclust::Mclust(
        x[!neg],
        G = 1, modelNames = "V",
        initialization = list(noise = pos_noise), verbose = FALSE
      )

      pos_gmmn_dens <- stats::dnorm(
        x[!neg], pos_gmmn$parameters$mean, sqrt(pos_gmmn$parameters$variance$sigmasq)
      )
      pos_gmmn_unif <- 1 / pos_gmmn$hypvol

      pos_gmmn_unif_prop <- sum(pos_gmmn$classification == 0) / sum(!neg)
      pos_gmmn_norm_prop <- 1 - pos_gmmn_unif_prop
    } else {
      pos_gmmn_dens <- rep(0, sum(!neg))

      pos_gmmn_unif_prop <- 0.5
      pos_gmmn_norm_prop <- 1 - pos_gmmn_unif_prop
    }

    neg_unif_ll <- sum(neg) * log(neg_const)
    pos_unif_ll <- sum(!neg) * log(pos_const)

    neg_unif_bic <- 2 * neg_unif_ll - log(sum(neg))
    pos_unif_bic <- 2 * pos_unif_ll - log(sum(!neg))

    ll <- c(
      "gmm_gmm" = sum(log(neg_dens)) + sum(log(pos_dens)),
      "gmm_gmmn" = sum(log(neg_dens)) +
          sum(log(pos_gmmn_norm_prop * pos_gmmn_dens + pos_gmmn_unif_prop * rep(pos_const, sum(!neg)))),
      "gmmn_gmm" =
        sum(log((neg_gmmn_norm_prop * neg_gmmn_dens + neg_gmmn_unif_prop * rep(neg_const, sum(neg))))) +
          sum(log(pos_dens)),
      "gmmn_gmmn" =
        sum(log((neg_gmmn_norm_prop * neg_gmmn_dens + neg_gmmn_unif_prop * rep(neg_const, sum(neg))))) +
          sum(log(pos_gmmn_norm_prop * pos_gmmn_dens + pos_gmmn_unif_prop * rep(pos_const, sum(!neg)))),
      "gmm_unif" = sum(log(neg_dens)) + sum(log(rep(pos_const, sum(!neg)))),
      "unif_gmm" = sum(log(rep(neg_const, sum(neg)))) + sum(log(pos_dens)),
      "gmmn_unif" =
        sum(log((neg_gmmn_norm_prop * neg_gmmn_dens + neg_gmmn_unif_prop * rep(neg_const, sum(neg))))) +
          sum(log(rep(pos_const, sum(!neg)))),
      "unif_gmmn" = sum(log(rep(neg_const, sum(neg)))) +
          sum(log(pos_gmmn_norm_prop * pos_gmmn_dens + pos_gmmn_unif_prop * rep(pos_const, sum(!neg)))),
      "unif_unif" = sum(log(rep(neg_const, sum(neg)))) + sum(log(rep(pos_const, sum(!neg))))
    )

    bic <- 2 * ll - npar * log(obs_num)

    best_bic[j] <- max(bic)
    best_model[j] <- which.max(bic)

    if (verbose) {
      message("*", appendLF = FALSE)
      if (j %% 10 == 0) {
        message("\tSplit No. ", j)
      }
    }
  }

  one_dens <- stats::dnorm(x, mean(x), stats::sd(x))
  one_const <- mclust::hypvol(x, TRUE)
  one_noise <- one_dens < one_const
  one_gmmn <- mclust::Mclust(
    x, 1, "V", initialization = list(noise = one_noise), verbose = FALSE
  )
  one_gmm_bic <- 2 * (sum(log(one_dens))) - 2 * log(obs_num)
  one_unif_bic <- 2 * (obs_num * log(one_const)) - log(obs_num)
  one_gmmn_bic <- one_gmmn$bic
  one_bic_vec <- c(
    "gmm" = one_gmm_bic, "unif" = one_unif_bic, "gmmn" = one_gmmn_bic
  )
  one_choice <- which.max(one_bic_vec)
  one_bic <- one_bic_vec[one_choice]

  choice <- which.max(best_bic)

  scaled_bic_diff <- (best_bic[choice] - one_bic) / (2 * log(obs_num))

  if (verbose) {
    message("Chosen one-component model was ", names(one_bic_vec)[one_choice])
    message("Chosen two-component model was ", names(npar)[best_model[choice]])
  }
  if (scaled_bic_diff > min_scaled_bic_diff) {
    boundary <- xseq[choice]
  } else {
    boundary <- NA
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
find_boundary4 <- function(x, min_scaled_bic_diff = 0, verbose = FALSE) {
  obs_num <- length(x)
  sortx <- sort(x)

  leftx <- sortx[10]
  rightx <- sortx[obs_num - 11]
  xseq <- seq(leftx, rightx, length.out = 100)
  model_num <- length(xseq)

  npar <- c(
    "gmm_gmm" = 4,
    "gmm_gmmn" = 6,
    "gmmn_gmm" = 6,
    "gmmn_gmmn" = 8,
    "gmm_unif" = 3,
    "unif_gmm" = 3,
    "gmmn_unif" = 5,
    "unif_gmmn" = 5,
    "unif_unif" = 2
  )

  best_model <- integer(model_num)
  best_bic <- double(model_num)
  for (j in 1:model_num) {
    neg <- x <= xseq[j]

    neg_prop <- sum(neg) / obs_num
    pos_prop <- 1 - neg_prop

    neg_dens <- stats::dnorm(x, mean(x[neg]), sd(x[neg]))
    pos_dens <- stats::dnorm(x, mean(x[!neg]), sd(x[!neg]))

    neg_const <- mclust::hypvol(x[neg], TRUE)
    pos_const <- mclust::hypvol(x[!neg], TRUE)

    neg_noise <- neg_dens[neg] < neg_const
    pos_noise <- pos_dens[!neg] < pos_const

    if (sum(neg) - sum(neg_noise) > 5) {
      neg_gmmn <- mclust::Mclust(
        x[neg],
        G = 1, modelNames = "V",
        initialization = list(noise = neg_noise), verbose = FALSE
      )

      neg_gmmn_dens <- stats::dnorm(
        x, neg_gmmn$parameters$mean, sqrt(neg_gmmn$parameters$variance$sigmasq)
      )

      neg_gmmn_unif_prop <- sum(neg_gmmn$classification == 0) / sum(neg)
      neg_gmmn_norm_prop <- 1 - neg_gmmn_unif_prop
    } else {
      neg_gmmn_dens <- 0

      neg_gmmn_unif_prop <- 0.5
      neg_gmmn_norm_prop <- 1 - neg_gmmn_unif_prop
    }

    if (sum(!neg) - sum(pos_noise) > 5) {
      pos_gmmn <- mclust::Mclust(
        x[!neg],
        G = 1, modelNames = "V",
        initialization = list(noise = pos_noise), verbose = FALSE
      )

      pos_gmmn_dens <- stats::dnorm(
        x, pos_gmmn$parameters$mean, sqrt(pos_gmmn$parameters$variance$sigmasq)
      )
      pos_gmmn_unif <- 1 / pos_gmmn$hypvol

      pos_gmmn_unif_prop <- sum(pos_gmmn$classification == 0) / sum(!neg)
      pos_gmmn_norm_prop <- 1 - pos_gmmn_unif_prop
    } else {
      pos_gmmn_dens <- 0

      pos_gmmn_unif_prop <- 0.5
      pos_gmmn_norm_prop <- 1 - pos_gmmn_unif_prop
    }

    neg_unif_ll <- sum(neg) * log(neg_const)
    pos_unif_ll <- sum(!neg) * log(pos_const)

    neg_unif_bic <- 2 * neg_unif_ll - log(sum(neg))
    pos_unif_bic <- 2 * pos_unif_ll - log(sum(!neg))

    ll <- c(
      "gmm_gmm" = sum(log(neg_prop * neg_dens + pos_prop * pos_dens)),
      "gmm_gmmn" = sum(log(
        neg_prop * neg_dens +
          pos_prop *
          (pos_gmmn_norm_prop * pos_gmmn_dens +
             pos_gmmn_unif_prop * pos_const)
      )),
      "gmmn_gmm" = sum(log(
        neg_prop *
          (neg_gmmn_norm_prop * neg_gmmn_dens +
             neg_gmmn_unif_prop * neg_const) +
        pos_prop * pos_dens
      )),
      "gmmn_gmmn" = sum(log(
        neg_prop *
          (neg_gmmn_norm_prop * neg_gmmn_dens +
             neg_gmmn_unif_prop * neg_const) +
          pos_prop *
          (pos_gmmn_norm_prop * pos_gmmn_dens +
             pos_gmmn_unif_prop * pos_const)
      )),
      "gmm_unif" = sum(log(
        neg_prop * neg_dens + pos_prop * rep(pos_const, obs_num)
      )),
      "unif_gmm" = sum(log(
         neg_prop * rep(neg_const, obs_num) + pos_prop * pos_dens
      )),
      "gmmn_unif" = sum(log(
        neg_prop *
          (neg_gmmn_norm_prop * neg_gmmn_dens +
             neg_gmmn_unif_prop * neg_const) +
          pos_prop * rep(pos_const, obs_num)
      )),
      "unif_gmmn" = sum(log(
        neg_prop * rep(neg_const, obs_num) +
          pos_prop *
          (pos_gmmn_norm_prop * pos_gmmn_dens +
             pos_gmmn_unif_prop * pos_const)
      )),
      "unif_unif" = sum(log(
        neg_prop * rep(neg_const, obs_num) + pos_prop * rep(pos_const, obs_num)
      ))
    )

    bic <- 2 * ll - npar * log(obs_num)

    best_bic[j] <- max(bic)
    best_model[j] <- which.max(bic)

    if (verbose) {
      message("*", appendLF = FALSE)
      if (j %% 10 == 0) {
        message("\tSplit No. ", j)
      }
    }
  }

  one_dens <- stats::dnorm(x, mean(x), stats::sd(x))
  one_const <- mclust::hypvol(x, TRUE)
  one_noise <- one_dens < one_const
  one_gmmn <- mclust::Mclust(
    x, 1, "V", initialization = list(noise = one_noise), verbose = FALSE
  )
  one_gmm_bic <- 2 * (sum(log(one_dens))) - 2 * log(obs_num)
  one_unif_bic <- 2 * (obs_num * log(one_const)) - log(obs_num)
  one_gmmn_bic <- one_gmmn$bic
  one_bic_vec <- c(
    "gmm" = one_gmm_bic, "unif" = one_unif_bic, "gmmn" = one_gmmn_bic
  )
  one_choice <- which.max(one_bic_vec)
  one_bic <- one_bic_vec[one_choice]

  choice <- which.max(best_bic)

  scaled_bic_diff <- (best_bic[choice] - one_bic) / (2 * log(obs_num))

  if (verbose) {
    message("Chosen one-component model was ", names(one_bic_vec)[one_choice])
    message("Chosen two-component model was ", names(npar)[best_model[choice]])
  }
  if (scaled_bic_diff > min_scaled_bic_diff) {
    boundary <- xseq[choice]
  } else {
    boundary <- NA
  }

  return(c(boundary, as.numeric(scaled_bic_diff)))
}
