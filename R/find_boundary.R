find_boundary <- function(x, min_scaled_bic_diff = 0, return_all = FALSE) {
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

  if (return_all) {
    return(list(ll_vec = ll,
                ll_max = ll_max,
                ll_one = ll_one,
                bic_one = bic_one,
                bic_two = bic_two,
                boundary = xseq[which.max(ll)]))
  } else if (scaled_bic_diff > min_scaled_bic_diff) {
    return(c(xseq[which.max(ll)], scaled_bic_diff))
  } else {
    return(c(NA, NA))
  }

}
