find_split_2d <- function(
    x, y, positive, mechanism, boundary_noise_comp = TRUE) {
  sign_mult <- ifelse(positive, -1, +1)
  
  theta <- seq(0, pi / 2, length.out = 10)[-c(1, 10)]
  alpha <- tan(theta) / (1 + tan(theta))
  split_mat <- matrix(nrow = length(theta), ncol = 4)
  for (i in seq_along(theta)) {
    z <- alpha[i] * iris[, 1] + sign_mult * (1 - alpha[i]) * iris[, 2]
    if (mechanism == "valley") {
      split_mat[i, 1:2] <- find_valley(z)
    } else if (mechanism == "boundary") {
      split_mat[i, 1:2] <- find_boundary(z)
    }
  }
  split_mat[, 3] <- - split_mat[, 1] / (1 - alpha)
  split_mat[, 4] <- alpha / (1 - alpha)
  
  max_ind <- which.max(split_mat[, 2])
  
  split_mat[max_ind, ]
}

find_valley_2d <- function(x, y, positive) {
  find_split_2d(x, y, positive, "valley")
}

find_boundary_2d <- function(x, y, positive, boundary_noise_comp = TRUE) {
  find_split_2d(x, y, positive, "boundary", boundary_noise_comp)
}
