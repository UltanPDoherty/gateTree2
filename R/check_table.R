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
  values_check <- rep(NA, 2)
  values_check[1] <- all(!is.na(plusminus_table))
  values_check[2] <- all(as.matrix(plusminus_table) %in% c(-1, 0, +1))
  if (!values_check[1]) {
    cat(paste0("This table contains NAs.\n"))
  } else if (!values_check[2]) {
    cat(paste0("This table contains values other than -1, 0, & 1.\n"))
  }
  if (!all(values_check)) {
    return(list(
      pass = FALSE,
      indistinct_pops = NULL
    ))
  }

  check_table_inner(plusminus_table, initial = TRUE)

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

      new_pm_list <- append(
        new_pm_list,
        split(pm_list[[i]], pop_group_label[[i]])
      )

      new_result <- append(
        new_result,
        rep(result[i], length(unique(pop_group_label[[i]])))
      )
    }

    pm_list <- new_pm_list
    new_pm_list <- list()
  }

  final_check <- all(new_result == "success")
  indistinct_pops <- lapply(pm_list[new_result == "fail"], rownames)
  names(indistinct_pops) <- NULL

  indistinct_subtables <- pm_list[new_result == "fail"]
  names(indistinct_subtables) <- NULL

  if (final_check) {
    cat("All populations are theoretically distinguishable.")
  } else {
    cat("The following groups of populations are indistinguishable:\n")
    lapply(
      seq_along(indistinct_pops),
      function(x) {
        cat(paste0(
          x, ". ",
          paste0(indistinct_pops[[x]], collapse = ", "),
          ".\n"
        ))
      }
    )
    cat("\n")
    lapply(
      seq_along(indistinct_subtables),
      function(x) print(indistinct_subtables[[x]])
    )
    cat("\n")
  }

  return(list(
    pass = final_check,
    indistinct_pops = indistinct_pops,
    indistinct_subtables = indistinct_subtables
  ))
}

# ==============================================================================

#' @title Inner function used by [check_table].
#'
#' @description
#' Non-recursive / non-iterative check of a single `plusminus_table` sub-table.
#'
#' @inheritParams gatetree
#' @param initial Logical: whether this is the initial check before the loop.
#'
#' @return If `initial = FALSE`, then `check_table_inner` returns a list:
#' * pop_group_label: integer vector of population group labels
#' * result: `"success"`, `"continue"`, or `"fail"`.
#'
#' If `initial = TRUE`, `check_table_inner` returns `NULL`.
check_table_inner <- function(plusminus_table, initial = FALSE) {
  pop_num <- nrow(plusminus_table)

  common_variables <- apply(plusminus_table, 2, function(x) all(x != 0))
  common_num <- sum(common_variables)

  if (initial && common_num == 0) {
    cat(paste0(
      "gateTree requires at least one variable for which each of the ",
      "populations in the table are described as positive or negative.\n",
      "This table does not satisfy this requirement."
    ))
  }

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

  if (initial && common_num != 0 && result == "fail") {
    cat(paste0(
      "gateTree requires at least one variable for which each of the ",
      "populations in the table are described as positive or negative.\n",
      "This table does satisfy this requirement."
    ))
    cat(paste0(
      "However, for these variables, either all populations are defined to be ",
      "positive or all populations are defined to be negative.\n",
      "Therefore, it is not possible to differentiate between the populations."
    ))
  }

  if (initial) {
    return(NULL)
  } else {
    return(list(pop_group_label = pop_group_label, result = result))
  }
}
