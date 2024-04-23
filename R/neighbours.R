#' Title
#'
#' @param a1_neighbour_preglasso Add one neighbours
#' @param c1_neighbour_preglasso Change one neighbours
#' @param target_set Set of variables to explore to explain the observed evidence
#' @param target_set_length Number of variables in the target set
#' @param glasso_result_inv Result of gLasso
#'
#' @return List with two elements: all the neighbours, remaining neighbours after gLasso
#'
#' @keywords internal
neighbours <- function(a1_neighbour_preglasso, c1_neighbour_preglasso, target_set, target_set_length, glasso_result_inv) {

  preglasso_neighbours <- data.frame(neighbours = c(a1_neighbour_preglasso, c1_neighbour_preglasso))

  colnames(preglasso_neighbours) <- c('neighbours')
  preglasso_neighbours_rows <- nrow(preglasso_neighbours)

  all_preglasso_neighbours <- as.data.frame(matrix(NA_character_, nrow = preglasso_neighbours_rows, ncol = target_set_length))
  colnames(all_preglasso_neighbours) <- target_set

  glasso_subset <- cbind(all_preglasso_neighbours, glasso_inv = character(preglasso_neighbours_rows))

  neighbour_combos_split <- strsplit(preglasso_neighbours$neighbours, " ")

  for (i in seq_len(preglasso_neighbours_rows)) {
    combos_len <- length(neighbour_combos_split[[i]])
    combos_len_u <- length(unique(neighbour_combos_split[[i]]))

    rows <- neighbour_combos_split[[i]][seq_len(combos_len - 1)]
    cols <- neighbour_combos_split[[i]][combos_len]

    r_node <- stringi::stri_replace_last_regex(rows, "__.*$", "")
    r_state <- stringi::stri_replace_first_regex(rows, ".*__", "")

    c_node <- if (combos_len_u > 1) stringi::stri_replace_last_regex(cols, "__.*$", "") else NULL
    c_state <- if (combos_len_u > 1) stringi::stri_replace_first_regex(cols, ".*__", "") else NULL

    rc_node <- c(r_node, c_node)
    rc_state <- c(r_state, c_state)

    all_preglasso_neighbours[i, match(rc_node, colnames(all_preglasso_neighbours))] <- rc_state

    glasso_subset[i, ] <- all_preglasso_neighbours[i, ]
    glasso_subset[i, 'glasso_inv'] <- glasso_result_inv[paste(rows, collapse = " "), cols]
  }

  glasso_subset <- glasso_subset[glasso_subset[, "glasso_inv"] != 0, ]

  return(list(all_preglasso_neighbours = all_preglasso_neighbours,
              glasso_subset = glasso_subset))

}
