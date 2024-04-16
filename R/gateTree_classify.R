gateTree_classify <- function(x, gateTree_output) {

  pop_num <- nrow(gateTree_output$splits)
  var_num <- ncol(x)
  obs_num <- nrow(x)

  new_subsetter <- matrix(TRUE, nrow = obs_num, ncol = pop_num)
  for (k in seq_len(pop_num)) {
    for (p in seq_len(var_num)) {
      if (gateTree_output$signs[k, p] == 1) {
        new_subsetter[, k] <- new_subsetter[, k] &
          x[, p] > gateTree_output$splits[k, p]
      } else if (gateTree_output$signs[k, p] == -1) {
        new_subsetter[, k] <- new_subsetter[, k] &
          x[, p] <= gateTree_output$splits[k, p]
      }
    }
  }

  new_labels <- rep(NA, obs_num)
  new_unassigned <- rowSums(new_subsetter) == 0
  new_labels[new_unassigned] <- 0
  new_labels[!new_unassigned] <- unlist(apply(
    new_subsetter[!new_unassigned, ],
    MARGIN = 1,
    which
  ))

  return(new_labels)
}
