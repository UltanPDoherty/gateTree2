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
#' @return Matrix in which the columns contain each variable's boundary and its
#' scaled BIC difference.
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

  return(boundaries)
}
