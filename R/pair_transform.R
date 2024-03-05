pair_transform <- function(x, y, trans_type = TRUE, angle = NULL) {

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
  for(i in seq_along(theta)) {
    new_vars[[i]] <- signs[1] * alpha[i] * x + signs[2] * (1 - alpha[i]) * y
  }
  new_vars <- as.data.frame(new_vars)

  proposals <- propose_splits(new_vars,
                              subsetter_g = rep(TRUE, nrow(new_vars)),
                              splittable_vars_g = rep(TRUE, ncol(new_vars)),
                              min_size = 100, min_depth = 1,
                              use_boundaries = TRUE)

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
  } else (
    stop("Failed to find a split after any of the proposed transformations.\n")
  )
}
