#' Return minimal set based on dominance relations
#'
#' @param best_set Set of most relevant explanations
#'
#' @return Minimal explanations
#' @keywords internal
minimal_exp <- function(best_set) {

  # function to return minimal set based on dominance relations

  # best_set = set of most relevant explanations

  not_minimal_df <- data.frame()
  minimal_df <- data.frame()

  minimal_set <- data.frame()

  n <- nrow(best_set)

  best_set_targets <- best_set
  best_set_targets[, c('GBF')] <- list(NULL)

  best_set_cols <- colnames(best_set)
  # best_set_cols[, c('mre_size')] <- list(NULL)

  for (i in 1:n) {

    # x
    x <- best_set[i, ]
    x_n <- best_set_targets[i, ]
    x_names <- names(which(sapply(x_n, function(x) any(!is.na(x)))))
    x_length <- length(which(sapply(x_n, function(x) any(!is.na(x))))) # non NA length

    for (j in 1:n) {

      if (i == j) next

      # y
      y <- best_set[j,]
      y_n <- best_set_targets[j, ]
      y_names <- names(which(sapply(y_n, function(x) any(!is.na(x)))))
      y_length <- length(which(sapply(y_n, function(x) any(!is.na(x))))) # non NA length

      # dom_relation: dominance relation
      # 1: remove solutions that are strongly dominated by another
      # 2: minimal
      # 3: not in above categories

      # check if x_names is a subset of y_names
      x_subset_y_names <- all(x_names %in% y_names)
      # check if x contains less variables than y
      x_subset_y_length <- (x_length < y_length)
      # check if x_names is a superset of y_names
      x_superset_y_names <- all(y_names %in% x_names)
      # check if x contains more variables than y
      x_superset_y_length <- (x_length > y_length)

      # logical vector showing matching variables in x and y
      x_match_y <- (x_n == y_n)
      # count of all matches (TRUE)
      x_match_y_true <- sum(x_n == y_n, na.rm = TRUE)
      # count of all NA
      x_match_y_na <- sum(is.na(x_n == y_n))
      # count of all non-matches (FALSE)
      # We use this to check if x is a subset/superset of y.
      # If there are more than 0 FALSE, x not a subset/superset of y.
      x_match_y_false <- (length(x_match_y) - x_match_y_true - x_match_y_na)

      if (x_subset_y_names && x_match_y_false == 0) {
        if (x_subset_y_length) {
          # x subset of y
          if (x$GBF >= y$GBF) {
            # x strongly dominates y. keep x
            dom_relation = '2'
          } else {
            # x does not strongly dominate y. Discard x
            dom_relation = '1'
          }
        } else {
          # x not subset/superset of y
          # x and y consist of the exactly same variables, however, some states differ
          # x not dominated by another explanation. Keep x
          dom_relation <- 2
        }
      } else if (x_superset_y_names && x_match_y_false == 0 && x_superset_y_length) {
        # x superset of y
        if (x$GBF > y$GBF) {
          # x weakly dominates y. Keep x
          dom_relation <- '2'
        } else {
          # x does not weakly dominate y. Discard x
          dom_relation <- '1'
        }
      } else {
        # x not subset/superset of y
        # x and y differ by variables + states.
        # x not dominated by another explanation. Keep x
        dom_relation <- '2'
      }

      if (dom_relation == '1') {
        not_minimal_df <- plyr::rbind.fill(not_minimal_df, x)
      } else if (dom_relation == '2') {
        minimal_df <- plyr::rbind.fill(minimal_df, x)
      }

      # remove duplicated
      not_minimal_df <- unique(not_minimal_df)
      minimal_df <- unique(minimal_df)

      if (nrow(not_minimal_df) > 0 && nrow(minimal_df) > 0) {
        minimal_set <- dplyr::anti_join(minimal_df, not_minimal_df, by = best_set_cols)
      } else {
        minimal_set <- minimal_df
      }

    }
  }

  return(minimal_set)

}
