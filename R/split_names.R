#' Split names into nodes and states
#'
#' @param names
#'
#' @return A list with node names and states
#' @keywords internal
split_names <- function(names) {

  # function to split names into nodes and states

  node_parts <- gsub("__.*$", "", names)
  state_parts <- gsub("^[^_]*__", "", names)

  list(nodes = node_parts, states = state_parts)
}
