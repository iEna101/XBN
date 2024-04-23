#' Brute force search to obtain all explanations according to generalised Bayes factor
#'
#' @param target_set Set of variables to explore to explain the observed evidence
#' @param evidence_set Set of observed variable names
#' @param evidence_states Set of observed variable states
#' @param bn_grain Bayesian network as grain object
#'
#' @return Dataframe containing all explanations
#' @export
#'

mre_brute <- function(target_set, evidence_set, evidence_states, bn_grain) {

  # brute force function to obtain all explanations

  gbf_table <- c()

  # set observed evidence
  e_1 <- gRain::setEvidence(bn_grain, nodes = evidence_set, states = evidence_states)

  for (i in 1:length(target_set)) {
    X1 <- gtools::combinations(n = length(target_set), r = i, v = target_set, set = T,
                               repeats.allowed = F)
    node_record <- data.frame(X1)

    for (j in 1:nrow(node_record)) {

      node_input <- as.character(node_record[j,])

      # numerator: observed hypotheses
      # query, i.e., obtain the conditional distribution
      odds_1 <- odds(gRain::querygrain(e_1, nodes = node_input, type = "joint"))

      # denominator: all alternative hypotheses
      # query, i.e., obtain the conditional distribution
      odds_2 <- odds(gRain::querygrain(bn_grain, nodes = node_input, type = "joint"))

      # generalised bayes factor
      gbf <- (odds_1 / odds_2)

      gbf_df <- as.data.frame.table(gbf)
      gbf_table <- plyr::rbind.fill(gbf_table, gbf_df)

      # Reorder the columns, moving the specified column to the last position
      gbf_table <- gbf_table[, c(setdiff(names(gbf_table), 'Freq'), 'Freq')]

    }
  }

  # rename `Freq` to `GBF`
  names(gbf_table)[names(gbf_table) == 'Freq'] <- 'GBF'

  # arrange the starting solutions in descending order
  gbf_table <- gbf_table %>%
    dplyr::distinct() %>%
    tidyr::drop_na(GBF) %>%
    dplyr::arrange(dplyr::desc(GBF)) %>%
    dplyr::select(-GBF,dplyr::everything()) %>%
    dplyr::mutate(mre_size = rowSums(!is.na(.)) - 1)

  gbf_table <- gbf_table %>%
    dplyr::mutate_at(dplyr::vars(dplyr::all_of(target_set)), as.character)

  return(gbf_table)

}
