#' Compare current best solution with best neighbour
#'
#' @param y_gbf Current best solution
#' @param neighbour_gbf Best neighbouring solution
#'
#' @return Updated best solution
#' @keywords internal
compare_gbf <- function(y_gbf, neighbour_gbf) {

  # function to compare the current best solution with best neighbour

  # y_gbf: current best solution
  # neighbour_gbf: top neighbouring solution

  gbf_sol <- data.frame()

  # if generalised Bayes factor of y_gbf is greater than or equal to
  # generalised Bayes factor of neighbour_gbf
  if (y_gbf$GBF >= neighbour_gbf[1,]$GBF) {
    # keep y_gbf as best solution
    updated_y <- y_gbf
  }
  else {
    # update current best to neighbour_gbf
    updated_y <- neighbour_gbf[1,]
  }

  # populate updated best solution
  gbf_sol <- plyr::rbind.fill(gbf_sol, updated_y)

  # remove duplicates
  gbf_sol <- gbf_sol[!duplicated(gbf_sol),]

  return(gbf_sol)

}
