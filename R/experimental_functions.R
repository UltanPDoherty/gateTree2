#' @importFrom mvtnorm dmvnorm
#' @importFrom abind abind
gmm_estep <- function(x, labels, exclude_index = NULL) {
  obs_num <- nrow(x)
  var_num <- ncol(x)
  label_num <- length(unique(labels))

  x_list <- lapply(split(x, labels), matrix, ncol = var_num)

  mu_list    <- lapply(x_list, colMeans)
  sigma_list <- lapply(x_list, stats::cov)
  counts_vec <- as.numeric(table(labels))

  include_index <- setdiff(seq_len(label_num), exclude_index)

  params <- list()
  params$mu <- Reduce(rbind, mu_list[include_index])
  params$sigma <- abind::abind(sigma_list[include_index], along = 3)
  params$prop <- counts_vec[include_index] / sum(counts_vec[include_index])

  clust_num <- length(include_index)

  lpdf <- vapply(
    1:clust_num,
    FUN.VALUE = double(obs_num),
    function(k) {
      mvtnorm::dmvnorm(x, log = TRUE,
                       mean = params$mu[k, ],
                       sigma = params$sigma[, , k])
    }
  )


  # for loop computes the block posterior probability matrix, postprob_block
  unnorm <- postprob <- matrix(nrow = obs_num, ncol = clust_num)
  log_maxes <- loglike_vec <- unnorm_sums <- vector(mode = "numeric",
                                                    length = obs_num)

  for (l in 1:obs_num) {
    # Add the log mixing proportions and then un-log this sum with exp.
    # Subtract lpdf_block row maxes to prevent exp mapping large values to Inf.
    log_maxes[l]      <- max(lpdf[l, ] + log(params$prop))
    unnorm[l, ] <- exp(lpdf[l, ] + log(params$prop) - log_maxes[l])

    # Normalise rows of block_unnorm to obtain postprob_block.
    unnorm_sums[l] <- sum(unnorm[l, ])
    postprob[l, ] <- unnorm[l, ] / unnorm_sums[l]

    loglike_vec[l] <- log(unnorm_sums[l]) + log_maxes[l]
  }

  loglike  <- sum(loglike_vec)

  return(list(loglike = loglike,
              postprob = postprob,
              lpdf = lpdf))
}

#===============================================================================

unassigned_labelling <- function(z, pdf, z_threshold, pdf_threshold) {

  map_labels <- apply(z, 1, which.max)

  zmap_values <- pdfmap_values <- c()
  low_zmap <- low_pdfmap <- c()
  map_labels_with0 <- c()
  for (i in seq_len(nrow(z))) {
    zmap_values[i] <- z[i, map_labels[i]]
    pdfmap_values[i] <- pdf[i, map_labels[i]]

    low_zmap[i] <- zmap_values[i] < z_threshold
    low_pdfmap[i] <- pdfmap_values[i] < pdf_threshold

    map_labels_with0[i] <- ifelse(low_zmap[i] | low_pdfmap[i], 0,
                                  map_labels[i])
  }

  return(data.frame(map = map_labels,
                    map0 = map_labels_with0,
                    z_vals = zmap_values,
                    pdf_vals = pdfmap_values,
                    low_zmap = low_zmap,
                    low_pdfmap = low_pdfmap))
}

#===============================================================================

pair_transform <- function(
    x, y, trans_type, angle = NULL,
    min_depth = 1, use_boundaries = TRUE
) {

  if (is.null(angle)) {
    angle <- seq(0, 90, by = 7.5)
  } else if (any(angle < 0)) {
    warning("angle values are expected to be in [0, 90].\n")
  }

  theta <- angle * pi / 180
  alpha <- tan(theta) / (1 + tan(theta))

  signs <- switch(trans_type,
                  "top_right"       = c(+1, +1),
                  "top_left"        = c(-1, +1),
                  "bottom_right"    = c(+1, -1),
                  "bottom_left"     = c(-1, -1),
                  "diag_increasing" = c(+1, -1),
                  "diag_decreasing" = c(+1, +1))

  new_vars <- list()
  for (i in seq_along(theta)) {
    new_vars[[i]] <- signs[1] * alpha[i] * x + signs[2] * (1 - alpha[i]) * y
  }
  new_vars <- as.data.frame(new_vars)

  proposals <- propose_splits(new_vars,
                              subsetter_g = rep(TRUE, nrow(new_vars)),
                              splittable_vars_g = rep(TRUE, ncol(new_vars)),
                              min_size = 100, min_depth = min_depth,
                              use_boundaries = use_boundaries)

  if (proposals$scenario != "nothing") {
    trans_choice <- which.max(proposals$matrix[2, ])

    cat(paste0("angle = ", signs[2] * angle[trans_choice], ".\n"))
    cat(paste0("Chosen split was a ", proposals$scenario,
               " at ", proposals$matrix[1, trans_choice],
               " with score ", proposals$matrix[2, trans_choice], ".\n"))

    if (angle[trans_choice] == 0) {
      warning("x variable was ignored.\n")
    } else if (angle[trans_choice] == 90) {
      warning("y variable was ignored.\n")
    }

    return(new_vars[, trans_choice])
  } else {
    stop("Failed to find a split after any of the proposed transformations.\n")
  }
}

