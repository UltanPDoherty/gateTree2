split_gmm <- function(x){
  N <- length(x)
  sortx <- sort(x)

  leftx <- sortx[10]
  rightx <- sortx[N - 11]
  xseq <- seq(leftx, rightx, length.out = 100)
  J <- length(xseq)

  prop <- double(2)
  mu <- double(2)
  sigma <- double(2)
  label <- logical(2)
  ll <- double(J)
  comp_mixpdf <- matrix(nrow = N, ncol = 2)
  for (j in 1:J) {
    neg <- x <= xseq[j]
    prop[1] <- sum(neg) / N
    prop[2] <- 1 - prop[1]
    mu[1] <- mean(x[neg])
    mu[2] <- mean(x[!neg])
    sigma[1] <- stats::sd(x[neg])
    sigma[2] <- stats::sd(x[!neg])
    # for (k in 1:2) {
    #   comp_mixpdf[, k] <- prop[k] * dnorm(x, mu[k], sigma[k])
    # }
    comp_mixpdf[, 1] <- prop[1] * truncnorm::dtruncnorm(x, -Inf, xseq[j], mu[1], sigma[1])
    comp_mixpdf[, 2] <- prop[2] * truncnorm::dtruncnorm(x, xseq[j], +Inf, mu[2], sigma[2])
    ll[j] <- sum(log(rowSums(comp_mixpdf)))
  }

  ll_max = max(ll)
  ll_one <- sum(stats::dnorm(x, mean(x), stats::sd(x), log = TRUE))

  bic_one = 2 * log(N) - 2 * ll_one
  bic_two = 4 * log(N) - 2 * ll_max

  if (return_all) {
    return(list(ll_vec = ll,
                ll_max = ll_max,
                ll_one = ll_one,
                bic_one = bic_one,
                bic_two = bic_two,
                boundary = xseq[which.max(ll)]))
  } else if (bic_two < bic_one) {
    return(c(xseq[which.max(ll)], bic_two))
  } else {
    return(c(NA, NA))
  }

}
