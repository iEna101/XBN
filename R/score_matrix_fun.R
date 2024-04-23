#' Create score matrix for gLasso
#'
#' @param df_target_neighbours Set of neighbours
#' @param my_gbfs Generalised Bayes factor scores
#'
#' @return Score matrix for gLasso
#' @keywords internal
score_matrix_fun <- function(df_target_neighbours, my_gbfs) {

  # create score_matrix -- this is used to compute sigma for glasso
  # elements in df_target_neighbours$neighbours == row and column names of score_matrix
  matrix_cols <- df_target_neighbours$neighbours
  # the number of rows in df_target_neighbours will determine the dimension of score_matrix
  matrix_dim <- nrow(df_target_neighbours)
  # create score_matrix filled with zero's
  score_matrix <- matrix(0, nrow = matrix_dim,
                         ncol = matrix_dim,
                         dimnames = list(matrix_cols, matrix_cols))

  # Preprocess row and column names
  row_names <- strsplit(rownames(score_matrix), " ")
  # col_names <- row_names

  # Extract node and state information from names
  row_info <- lapply(row_names, split_names)

  # Calculate row sums outside the loop
  row_sums <- rowSums(!is.na(my_gbfs))

  # Iterate over rows and columns
  for (r in seq_len(matrix_dim)) {

    node_r <- row_info[[r]]$nodes
    state_r <- row_info[[r]]$states
    row_names_r <- row_names[[r]]

    for (c in r:matrix_dim) {

      node_c <- row_info[[c]]$nodes
      state_c <- row_info[[c]]$states
      col_names_c <- row_names[[c]]

      common_rows <- !is.na(match(row_names_r, col_names_c))
      common_cols <- !is.na(match(col_names_c, row_names_r))

      state_set <- c(state_r, state_c)
      node_set <- c(node_r, node_c)

      if (any(common_rows) || any(common_cols)) {

        length_diff_cr <- length(setdiff(col_names_c, row_names_r))

        if (length_diff_cr == 0) {
          state_set <- state_r
          node_set <- node_r
        } else if (length_diff_cr == 1) {
          rc <- unique(c(row_names_r, col_names_c))
          node_set <- substr(rc, 1, regexpr("__", rc) - 1)
          state_set <- substr(rc, regexpr("__", rc) + 1, nchar(rc))
        }
      } else if (any(node_c %in% node_r)) {
        state_set <- NULL
        node_set <- NULL
      }

      length_node_set <- length(node_set)
      c_match <- colnames(my_gbfs) %in% node_set

      row_match2 <- FALSE

      if (any(c_match)) {
        subset_all_combos_gbf <- my_gbfs[row_sums - 1 == length_node_set, c("GBF", colnames(my_gbfs)[c_match])]
        subset_all_combos_gbf <- subset_all_combos_gbf[stats::complete.cases(subset_all_combos_gbf), ]

        if (length_node_set != 1 || !any(node_r == node_c) || !any(state_r == state_c)) {
          row_match <- t(apply(subset_all_combos_gbf[, node_set], 1,
                               function(x) x == state_set))
          row_match2 <- rowSums(row_match) == length_node_set
        } else {
          row_match <- subset_all_combos_gbf[, node_set] == state_set
          row_match2 <- sapply(row_match, function(x) all(x))
        }

        if (any(row_match2)) {
          score_matrix[r, c] <- subset_all_combos_gbf[row_match2, "GBF"][1]

        }
      }

    }

  }

  # Assign values symmetrically
  score_matrix[lower.tri(score_matrix)] <- t(score_matrix)[lower.tri(score_matrix)]

  return(score_matrix)

}
