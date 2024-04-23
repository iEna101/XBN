#' Compute the initiliasation values for each target in the set
#'
#' @param target_set Set of variables to explore to explain the observed evidence
#' @param evidence_set Set of observed variable names
#' @param evidence_states Set of observed variable states
#' @param bn_grain Bayesian network as grain object
#'
#' @return Starting solutions for each target variable in target_set
#' @importFrom magrittr '%>%'
#' @export

init_gbf <- function(target_set, evidence_set, evidence_states, bn_grain) {

  # function to compute the starting values for each target in the set

  # create empty set to save gbf scores for each target
  gbf_results <- c()

  # set observed evidence
  e_1 <- gRain::setEvidence(bn_grain, nodes = evidence_set, states = evidence_states)

  # the `best pivot` rule requires each target variable to be in their most
  # optimal state.
  # for each starting solution set t:
  for (t in target_set) {

    # numerator: observed hypotheses
    # query, i.e., obtain the conditional distribution
    odds_1 <- odds(gRain::querygrain(e_1, nodes = t, type = 'joint'))

    # denominator: all alternative hypotheses
    # query, i.e., obtain the conditional distribution
    odds_2 <- odds(gRain::querygrain(bn_grain, nodes = t, type = 'joint'))

    # generalised bayes factor
    gbf <- odds_1 / odds_2
    gbf_df <- as.data.frame.table(gbf)
    best_gbf <- gbf_df[which.max(gbf_df$Freq),]

    gbf_results <- plyr::rbind.fill(gbf_results, best_gbf)
    # Reorder the columns, moving the specified column to the last position
    gbf_results <- gbf_results[, c(setdiff(names(gbf_results), 'Freq'), 'Freq')]

  }

  # rename `Freq` to `GBF`
  names(gbf_results)[names(gbf_results) == 'Freq'] <- 'GBF'

  # arrange the starting solutions in descending order
  gbf_results <- gbf_results %>%
    dplyr::arrange(dplyr::desc(GBF)) %>%

    # return the starting solutions for each target variable in target_set
    return(gbf_results)

}
