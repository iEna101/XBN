#' Compute odds ratio
#'
#' @param p Query
#'
#' @return Odds ratio
#' @keywords internal
#'
odds <- function(p) {

  # function to compute the odds

  p/(1 - p)

}
