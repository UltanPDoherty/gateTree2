#' @title Check if table populations are theoretically distinguishable.
#'
#' @description
#' If all splits were achieved, is there enough information in the user-provided
#' table to distinguish between all of the described populations.
#'
#' @inheritParams gatetree
#'
#' @return List:
#' * pass: Logical, whether the table passed the check.
#' * indistinct_pops: List of indistinguishable sub-tables.
#' @export
check_table <- function(plusminus_table) {

  pm_list <- list(plusminus_table)
  new_pm_list <- list()

  result <- rep("continue", length(pm_list))
  new_result <- c()

  while_count <- 0
  while (any(result == "continue") && while_count < 10) {
    while_count <- while_count + 1

    pop_group_label <- list()
    result <- rep(NA, length(pm_list))
    new_result <- c()

    for (i in seq_along(pm_list)) {
      basic_out <- check_table_inner(pm_list[[i]])
      pop_group_label[[i]] <- basic_out$pop_group_label
      result[i] <- basic_out$result

      new_pm_list <- append(new_pm_list,
                            split(pm_list[[i]], pop_group_label[[i]]))

      new_result <- append(new_result,
                           rep(result[i], length(unique(pop_group_label[[i]]))))
    }

    pm_list <- new_pm_list
    new_pm_list <- list()
  }

  final_check <- any(new_result == "fail")
  indistinct_pops <- lapply(pm_list[new_result == "fail"], rownames)

  if (final_check) {
    cat("The following groups of populations were indistinguishable:\n")
    lapply(indistinct_pops, function(x) cat(x, "\n"))
  } else {
    cat("All populations are theoretically distinguishable.")
  }

  return(list(pass = final_check, indistinct_pops = indistinct_pops))
}

#===============================================================================


#' @title Inner function used by [check_table].
#'
#' @description
#' Top-level check of a single plusminus_table sub-table.
#'
#' @inheritParams gatetree
#'
#' @return List
#' * pop_group_label: integer vector of population group labels
#' * result: `"success"`, `"continue"`, or `"fail"`.
check_table_inner <- function(plusminus_table) {

  pop_num <- nrow(plusminus_table)

  common_variables <- apply(plusminus_table, 2, function(x) all(x != 0))
  common_num <- sum(common_variables)

  pop_group_binary <- rep(0, pop_num)
  for (j in seq_len(common_num)) {
    pop_group_binary <- (
      pop_group_binary +
        2^j * (plusminus_table[, which(common_variables)[j]] == 1)
    )
  }

  pop_group_label <- as.numeric(as.factor(pop_group_binary))

  if (length(unique(pop_group_label)) == pop_num) {
    result <- "success"
  } else if (length(unique(pop_group_label)) > 1) {
    result <- "continue"
  } else {
    result <- "fail"
  }

  return(list(pop_group_label = pop_group_label, result = result))
}
