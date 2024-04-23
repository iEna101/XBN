#' Compute the generalised Bayes factor score for each solution in the solution set (sol_set)
#'
#' @param sol_set Solution set
#' @param evidence_set Set of observed variable names
#' @param evidence_states Set of observed variable states
#' @param bn_grain Bayesian network as grain object
#'
#' @return Generalised Bayes factor scores for all variables in the solution set (sol_set)
#' @keywords internal

gbf_set = function(sol_set, evidence_set, evidence_states, bn_grain) {

  # function to compute the gbf scores for each solution in sol_set

  # sol_set = set of variables we want to explore to explain observed evidence
  # evidence_set = set of observed variables, for example: dyspnoea
  # evidence_states = observed states of evidence variables, for example: true
  # bn_grain = Bayesian network as grain object

  # create empty set to save generalised bayes factor (gbf) scores for each target
  gbf_results <- c()

  # set observed evidence
  e_1 <- gRain::setEvidence(bn_grain, nodes = evidence_set, states = evidence_states)

  # for row in sol_set:
  for (j in 1:nrow(sol_set)) {

    # extract node j from sol_set
    node_input <- as.character(sol_set[j,])

    # numerator: observed hypotheses
    # query, i.e., obtain the conditional distribution
    odds_1 <- odds(gRain::querygrain(e_1, nodes = node_input, type = "joint"))

    # denominator: all alternative hypotheses
    # query, i.e., obtain the conditional distribution
    odds_2 <- odds(gRain::querygrain(bn_grain, nodes = node_input, type = "joint"))

    # generalised bayes factor
    gbf <- odds_1 / odds_2
    gbf_df <- as.data.frame.table(gbf)

    gbf_results <- plyr::rbind.fill(gbf_results, gbf_df)
    # Reorder the columns, moving the specified column to the last position
    gbf_results <- gbf_results[, c(setdiff(names(gbf_results), 'Freq'), 'Freq')]
  }

  # rename `Freq` to `GBF`
  names(gbf_results)[names(gbf_results) == 'Freq'] <- 'GBF'

  # return the gbf scores for all variables in sol_set
  return(gbf_results)

}
