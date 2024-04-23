#' Solve the most relevant explanation in Bayesian networks using a forward-gLasso search algorithm.
#'
#' @param target_set Set of variables to explore to explain the observed evidence
#' @param evidence_set Set of observed variable names
#' @param evidence_states Set of observed variable states
#' @param bn_grain Bayesian network as grain object
#' @param bn_rho L1 regularisation parameter for gLasso
#' @param score_scale Parameter to add very small random noise to score matrix for gLasso
#'
#' @return Set of most relevant explanations
#' @export
mre_fwd_glasso <- function(target_set, evidence_set, evidence_states, bn_grain,
                           bn_rho, score_scale = TRUE) {

  # Function to solve the most relevant explanation in Bayesian networks using
  # a Forward-gLasso search algorithm.

  # initialise the starting solution set, I
  # init.GBF based on best pivot rule
  i_set <- init_gbf(target_set = target_set,
                    evidence_set = evidence_set,
                    evidence_states = evidence_states,
                    bn_grain = bn_grain)

  # initialise the current best solution
  y_best <- data.frame()

  # create empty set to save all candidate solutions
  candidate_solutions <- data.frame()

  # obtain length of target set
  target_set_length <- length(target_set)

  # create empty set to save all combinations visited
  # certain combinations may be visited more than once
  all_combos_visit <- data.frame(matrix(NA, ncol = target_set_length))
  all_combos_visit_cnames <- colnames(all_combos_visit)

  # empty set for skipping
  empty_gbf <- data.frame(matrix(NA, ncol = ncol(i_set)))
  colnames(empty_gbf) <- colnames(i_set)

  # empty set for saving combination solutions
  my_gbfs <- c()

  # for each starting solution do:
  for (s.id in seq_len(nrow(i_set))) {

    # s: single starting solution in I
    s <- i_set[s.id,]
    # set y = s
    y <- s

    # the process will repeat until y no longer updates
    # repeat {}
    while (!is.infinite(y$GBF)) {

      # get variable names present in y
      # remove `GBF` column from y since we only want the names of the variables
      y_target_nogbf <- y[-length(y)]
      # obtain names of variables present in y. This will exclude those variable
      # names that are NA-cells
      y_target <- names(which(sapply(y_target_nogbf, function(x) any(!is.na(x)))))
      y_target <- as.data.frame(y_target)

      # get name + "_" + state of target node(s). This is used to set up score_matrix
      # create empty set to save all name + "_" + state combinations
      targets_name_state <- c()
      # convert y_target_nogbf list to vector
      y_target_vec <- y_target_nogbf[unlist(y_target)]
      for (col_name in colnames(y_target_vec)) {
        cell_value <- unlist(y_target_vec[col_name])
        # get name + "_" + state of target node
        target_name_state <- paste(names(cell_value), cell_value[1], sep = "__")
        # save all name + "_" + state of target nodes
        targets_name_state <- c(targets_name_state, target_name_state)
      }

      # concatenate elements in targets_name_state into single string
      target_name_state_s <- paste(targets_name_state, collapse = " ")

      # the (starting) solution is improved by either adding an additional feature
      # or by changing the state of an existing feature in the solution

      # `add one` node as neighbour:
      y_target_unlist <- unlist(y_target)
      # select nodes from target_set not present in y_target (y_best)
      a1_nodes <- setdiff(target_set, y_target_unlist)

      # create empty vector to save all `add one` neighbours
      a1_neighbours <- character()

      # add states for each node in a1_nodes
      for (i in a1_nodes) {
        # extract states from bn_grain for node i
        a1_state <- bn_grain$universe$levels[[i]]
        # concatenate node i + "_" + a1_state
        a1_node_state <- paste(i, a1_state, sep = "__")
        # append to a1_neighbours vector
        a1_neighbours <- c(a1_neighbours, a1_node_state)
      }

      # will be used to set up preglasso neighbours
      a1_neighbour_preglasso <- paste(target_name_state_s, a1_neighbours, sep = " ")

      # save target node(s) and `add one` neighbours into a data.frame
      df_target_neighbours <- data.frame(c(target_name_state_s, t(a1_neighbours)))
      colnames(df_target_neighbours) <- c('neighbours')

      # `change one` node as a neighbour
      # extract the states from bn_grain for all nodes
      c1_states <- lapply(names(y_target_vec), function(node) bn_grain$universe$levels[[node]])

      # create empty sets for `change one`.diff and .remain
      # we are `splitting` the `change one` neighbours into two separate elements
      # such that we can obtain the partial correlation coefficients in glasso
      c1_neighbours_diff <- c()
      c1_neighbours_remain <- c()

      # for each node in the target set
      for (i in seq_len(ncol(y_target_vec))) {

        # make a copy of y_target_vec
        c1_state_change <- y_target_vec
        # extract states not included in c1_state_change
        c1_state_diff <- setdiff(c1_states[[i]], c1_state_change[, i])

        # iterate through all states in c1_state_diff
        for (j in seq_along(c1_state_diff)) {
          c1_state_change[, i] <- c1_state_diff[j]
          # concatenate node + "_" + state
          c1_neighbour <- sprintf("%s__%s", names(c1_state_change), format(c1_state_change))

          # extract `change one`.diff
          c1_neighbours_split_diff = setdiff(c1_neighbour, targets_name_state)
          # extract `change one`.remain
          c1_neighbours_split_remain <- c1_neighbours_split_diff
          if (length(c1_neighbour) > 1) {
            c1_neighbours_split_remain <- intersect(c1_neighbour, targets_name_state)
          }

          # concatenate elements in c1_neighbours_split_diff into single string
          c1_neighbour_name_state_diff <- paste(c1_neighbours_split_diff, collapse = " ")
          # concatenate elements in c1_neighbours_split_remain into single string
          c1_neighbour_name_state_remain <- paste(c1_neighbours_split_remain, collapse = " ")

          # save `change one`.diff in a data.frame of its own
          c1_neighbours_diff <- c(c1_neighbours_diff, c1_neighbour_name_state_diff)
          # save `change one`remain in a data.frame of its own
          c1_neighbours_remain <- c(c1_neighbours_remain, c1_neighbour_name_state_remain)

          # save updated `change one` neighbours in a data.frame
          c1_neighbour_df <- data.frame(c(c1_neighbour_name_state_diff, c1_neighbour_name_state_remain))
          colnames(c1_neighbour_df) <- c('neighbours')
          # save `change one` neighbours into the same data.frame as target and `add one` neighbours
          df_target_neighbours <- rbind(df_target_neighbours, c1_neighbour_df)

        }

      }

      df_target_neighbours <- unique(df_target_neighbours)

      c1_neighbour_preglasso <- paste(c1_neighbours_remain, c1_neighbours_diff, sep = " ")

      # split neighbours and extract variable names
      neighbour_names <- unique(sapply(strsplit(df_target_neighbours$neighbours, " "),
                                       function(x) paste0(sub("__.*$", "", x), collapse = " ")))

      # generate data.frame with all combinations
      neighbour_names_df <- unique(data.frame(t(utils::combn(neighbour_names, 2))))

      # concatenate the elements of neighbour_names_df
      neighbour_names_all <- paste(neighbour_names_df[,1], neighbour_names_df[,2], sep = " ")

      neighbour_names_all <- c(neighbour_names_all, y_target_unlist, a1_nodes)
      names(neighbour_names_all) <- NULL

      # Split and obtain unique sorted neighbors
      neighbour_names_all <- stringr::str_split(neighbour_names_all, " ")
      neighbour_names_unique <- lapply(neighbour_names_all, function(vec) unique(stringr::str_sort(vec)))
      neighbour_names_unique <- unique(neighbour_names_unique)

      all_combos <- data.frame(matrix(NA, ncol = target_set_length, nrow = length(neighbour_names_unique)))

      # Calculate the number of unique neighbor combinations
      num_combos <- length(neighbour_names_unique)

      for (i in seq_len(num_combos)) {
        n <- length(neighbour_names_unique[[i]])
        all_combos[i, 1:n] <- neighbour_names_unique[[i]]
      }

      all_combos <- setdiff(all_combos, all_combos_visit)

      # update combinations already visited
      all_combos_visit <- unique(rbind(all_combos_visit, all_combos))

      all_combos <- as.data.frame(all_combos)
      # compute gbf score for combinations not yet visited
      filtered_combos_gbf <- empty_gbf
      if (nrow(all_combos) > 0) {
        filtered_combos_gbf <- gbf_set(sol_set = all_combos,
                                       evidence_set = evidence_set,
                                       evidence_states = evidence_states,
                                       bn_grain = bn_grain)
      }

      filtered_combos_gbf <- unique(filtered_combos_gbf)
      my_gbfs <- unique(plyr::rbind.fill(my_gbfs, filtered_combos_gbf))
      my_gbfs <- my_gbfs[!apply(is.na(my_gbfs), 1, all), ]

      # generalised bayes factor matrix
      score_matrix <- score_matrix_fun(df_target_neighbours, my_gbfs)

      # Replace Inf values with the chosen Bayes factor value
      score_matrix[is.infinite(score_matrix)] <- max(my_gbfs$GBF[!is.infinite(my_gbfs$GBF)]) + 1

      # Replace NaN values with zero
      score_matrix[is.na(score_matrix)] <- stats::rnorm(1, mean = 0, sd = 0.00001)

      if (score_scale == TRUE) {
        # add very small random noise to zero diagonal values
        diag_indices <- which(diag(score_matrix) == 0)
        noisy_values <- stats::rnorm(length(diag_indices), mean = 0, sd = 0.00001)

        score_matrix[diag_indices, diag_indices] <- noisy_values

        score_matrix <- scale(score_matrix)
      }

      # obtain covariance matrix, sigma
      sigma <- stats::cov(score_matrix)

      # apply glasso
      glasso_result <- glassoFast::glassoFast(sigma, rho = bn_rho)
      glasso_result_inv <- glasso_result$wi

      # rename rows and columns
      rownames(glasso_result_inv) <- rownames(score_matrix)
      colnames(glasso_result_inv) <- colnames(score_matrix)

      neighbour_result <- neighbours(a1_neighbour_preglasso, c1_neighbour_preglasso, target_set, target_set_length, glasso_result_inv)
      glasso_subset <- neighbour_result$glasso_subset
      all_preglasso_neighbours <- neighbour_result$all_preglasso_neighbours

      df_neighbours_preglasso <- dplyr::left_join(all_preglasso_neighbours, my_gbfs, by = target_set)
      # shrink neighbours using results from glasso_subset
      df_neighbours_update <- dplyr::semi_join(df_neighbours_preglasso, glasso_subset, by = target_set)

      df_neighbours_update <- df_neighbours_update[!is.na(df_neighbours_update$GBF), ]

      # check minimal explanations, adding candidate solution here to see if any
      # of the neighbours are dominated by a candidate solution. This will
      # ensure minimal explanations, we later remove the candidate solutions
      # form this list
      df_neighbours_update <- dplyr::bind_rows(dplyr::filter(df_neighbours_update, !is.nan(GBF)), candidate_solutions)
      df_neighbours_update <- minimal_exp(df_neighbours_update)

      # check if neighbouring set is in candidate solutions. If so, remove set from
      # neighbouring set
      if (nrow(candidate_solutions) >= 1) {

        df_neighbours_update <- dplyr::anti_join(df_neighbours_update, candidate_solutions, by = target_set)
      }

      df_neighbours_update <- df_neighbours_update[order(-df_neighbours_update$GBF), ]

      # compare y_best with neighbours
      if (nrow(df_neighbours_update) >= 1) {
        y_update <- compare_gbf(y, df_neighbours_update[1,])
      } else {
        y_update = y
      }

      # save previous y for stop criteria
      y_old <- y
      # update y
      y <- y_update

      # save all candidate solutions visited
      candidate_solutions <- unique(plyr::rbind.fill(candidate_solutions, y, y_old))

      # end repeat or while loop
      # Check if stop rule is satisfied or there are no missing values
      if (setequal(y, y_old) || sum(is.na(y)) == 0) {
        break
      }

    }

    # update y_best
    y_best <- unique(plyr::rbind.fill(y_best, y))

  }

  # arrange solution set according to generalised bayes factor score
  mre_set <- y_best %>%
    dplyr::distinct(GBF, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(GBF)) %>%
    dplyr::select(-GBF, dplyr::everything()) %>%
    dplyr::mutate(mre_size = rowSums(!is.na(.)) - 1)

  mre_set <- mre_set %>%
    dplyr::mutate_at(dplyr::vars(dplyr::all_of(target_set)), as.character)

  return(mre_set)

}
