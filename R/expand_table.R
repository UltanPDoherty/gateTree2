#' @title Add missing variables to a plusminus table.
#'
#' @description
#' Given a plusminus subtable which does not include all of the variables from
#' the data set, add in zero columns for the unused variables.
#'
#' @inheritParams gatetree
#' @param plusminus_subtable A `matrix` or `data.frame` whose column names are a
#' subset of the column names of `x`.
#'
#' @return A `data.frame` with the same number of columns as `x``.`
#' @export
expand_table <- function(x, plusminus_subtable) {
  pop_num <- nrow(plusminus_subtable)
  var_num <- ncol(x)

  plusminus_table <- as.data.frame(matrix(0, nrow = pop_num, ncol = var_num))
  colnames(plusminus_table) <- colnames(x)

  plusminus_table[, colnames(plusminus_subtable)] <- plusminus_subtable
  rownames(plusminus_table) <- rownames(plusminus_subtable)

  plusminus_table
}
