#' Solve the most relevant explanation in Bayesian networks using a forward search algorithm.
#'
#' @param target_set Set of variables to explore to explain the observed evidence
#' @param evidence_set Set of observed variable names
#' @param evidence_states Set of observed variable states
#' @param bn_grain Bayesian network as grain object
#'
#' @return Set of most relevant explanations
#' @export


mre_fwd <- function(target_set, evidence_set, evidence_states, bn_grain) {

  # Function to solve the most relevant explanation in Bayesian networks using
  # a Forward search algorithm.

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

  all_neighbours_calc_gbf <- c()

  # empty set for skipping
  empty_gbf <- data.frame(matrix(NA, ncol = ncol(i_set)))
  colnames(empty_gbf) <- colnames(i_set)

  # for each starting solution do:
  for (s_id in 1:nrow(i_set)) {

    # s: single starting solution in I
    s <- i_set[s_id,]
    # set y = s
    y <- s

    # the process will repeat until y no longer updates
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

      # save target + `add one`
      a1_neighb <- paste(paste(names(y_target_vec), collapse = " "), a1_nodes, sep = " ")

      # create empty set to save all `add one` neighbours
      a1_neighbours <- c()

      # add states for each node in a1_nodes
      for (i in a1_nodes) {
        # extract states from bn_grain for node i
        a1_state <- bn_grain$universe$levels[[i]]
        # concatenate node i + "_" + a1_state
        a1_node_state <- paste(i, a1_state, sep = "__")
        # concatenate target + a1_node_state
        a1_neighbour <- paste(target_name_state_s, a1_node_state, sep = " ")
        # save `add one` neighbours
        a1_neighbours <- c(a1_neighbours, a1_neighbour)
      }

      # `change one` node as a neighbour
      # create empty set to save all `change one` neighbours
      c1_neighbours <- c()

      # for each node in the target set
      for (i in 1:ncol(y_target_vec)) {
        # extract states from bn_grain for node i
        c1_state <- bn_grain$universe$levels[[names(y_target_vec[i])]]
        # make copy of y_target_vec
        c1_state_change <- y_target_vec
        # extract states not included in y_target_vec
        c1_state_diff <- setdiff(c1_state, y_target_vec[,i])

        # iterate through all states in c1_state_diff
        for (j in 1:length(c1_state_diff)) {

          c1_state_change[,i] <- c1_state_diff[j]
          # extract node name
          c1_nodes <- names(c1_state_change)
          # concatenate nodes in c1_nodes -- this is used to get gbf of neighbours
          c1_nodes_concat <- paste(c1_nodes, collapse = " ")
          # concatenate node + "_" + state
          c1_neighbour1 <- paste(c1_nodes, format(c1_state_change), sep = "__")
          # concatenate elements in c1_neighbour
          c1_neighbour <- paste(c1_neighbour1, collapse = " ")
          # save `change one` neighbours
          c1_neighbours <- rbind(c1_neighbours, c1_neighbour)
        }
      }

      # save `add one` and `change one` neighbours into a data.frame
      df_neighbours <- data.frame(c(t(a1_neighbours),
                                    t(c1_neighbours)))
      colnames(df_neighbours) <- c('neighbours')

      # keep only unique rows in df_neighbours
      df_neighbours <- unique(df_neighbours)

      # get nodes for generalised bayes factor calculations
      # create data.frame from empty matrix with,
      # number of columns  = number of target variables
      # number of rows = number of neighbours
      all_neighbours_calc <- data.frame(matrix(NA,
                                               ncol = target_set_length,
                                               nrow = (length(a1_nodes)) + length(c1_neighbour)))

      df_neighbours_calc <- data.frame(c(t(a1_neighb), c1_nodes_concat))
      colnames(df_neighbours_calc) <- c("neighbours")

      for (i in 1:nrow(all_neighbours_calc)) {
        n <- length(strsplit(df_neighbours_calc[i,], " ")[[1]])
        for (j in 1:n) {
          all_neighbours_calc[i,j] <- strsplit(df_neighbours_calc[i,], " ")[[1]][j]
        }
      }

      # compute gbf score for all neighbours
      all_neighbours_calc <- setdiff(all_neighbours_calc, all_combos_visit)

      # update combinations already visited
      all_combos_visit <- unique(rbind(all_combos_visit, all_neighbours_calc))

      # compute gbf score for combinations not yet visited
      filtered_combos_gbf <- empty_gbf
      all_neighbours_calc <- as.data.frame(all_neighbours_calc)

      if (nrow(all_neighbours_calc) > 0) {
        filtered_combos_gbf <- gbf_set(sol_set = all_neighbours_calc,
                                       evidence_set = evidence_set,
                                       evidence_states = evidence_states,
                                       bn_grain = bn_grain)
      }

      filtered_combos_gbf <- unique(filtered_combos_gbf)
      all_neighbours_calc_gbf <- unique(plyr::rbind.fill(all_neighbours_calc_gbf, filtered_combos_gbf))
      all_neighbours_calc_gbf <- all_neighbours_calc_gbf[!apply(is.na(all_neighbours_calc_gbf), 1, all), ]


      neighbour_subset <- data.frame(matrix(NA,
                                            ncol = target_set_length,
                                            nrow = (length(df_neighbours))))
      colnames(neighbour_subset) <- target_set

      for (i in 1:nrow(df_neighbours)) {

        # extract row
        neighbour_i <- strsplit(df_neighbours[i,], " ")[[1]]

        # extract nodes
        neighbour_i_n <- stringi::stri_replace_last_regex(neighbour_i, "__.*$", "")

        # extract states
        neighbour_i_s <- stringi::stri_replace_first_regex(neighbour_i, ".*__", "")

        neighbour_subset[i, neighbour_i_n] <- neighbour_i_s
      }


      # select only those values from all_neighbours_calc_gbf that occur in df_neighbours
      neighbour_subset_gbf <- dplyr::semi_join(all_neighbours_calc_gbf, neighbour_subset, by = target_set)

      neighbour_subset_gbf <- neighbour_subset_gbf[!is.na(neighbour_subset_gbf$GBF), ]

      # check minimal explanations, adding candidate solution here to see if any
      # of the neighbours are dominated by a candidate solution. This will
      # ensure minimal explanations, we later remove the candidate solutions
      # form this list
      all_neighbours_gbf <- dplyr::bind_rows(dplyr::filter(neighbour_subset_gbf, !is.nan(GBF)), candidate_solutions)
      all_neighbours_gbf <- minimal_exp(all_neighbours_gbf)

      # check if neighbouring set is in candidate solutions. If so, remove set from
      # neighbouring set
      if (nrow(candidate_solutions) >= 1) {
        all_neighbours_gbf <- dplyr::anti_join(all_neighbours_gbf, candidate_solutions, by = target_set)
      }

      # arrange neighbours according to generalised bayes factor
      all_neighbours_gbf <- all_neighbours_gbf[order(-all_neighbours_gbf$GBF), ]

      # compare y_best with neighbours
      if (nrow(all_neighbours_gbf) >= 1) {
        y_update <- compare_gbf(y, all_neighbours_gbf[1,])
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
