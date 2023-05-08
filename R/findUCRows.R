#' Find Unique and Common Rows
#'
#' This is a helper function that is designed to find those rows that are
#' unique for a given dataframe or that are shared among a list of dataframes.
#' This function looks at a provided column and if the values for that column
#' are unique in a certain dataframe, it stores those rows as unique for that
#' dataframe; on the other hand, if it finds common values among the provided
#' dataframes, it stores them in another list. Thus, it only needs two
#' arguments: the list of dataframes and the column in which it has to look
#' those values. It returns a list with the unique and common rows that can
#' be used by other functions.
#'
#' @param dataframes List of dataframes.
#' @param column Column in which we have to look the unique/common values.
#' @return A list with the unique and common rows.
#' @importFrom dplyr filter
#' @export


findUCRows <- function(dataframes, column) {
  # We create lists for storing the unique and common rows:
  unique_rows <- list()
  common_rows <- list()

  # Then, we iterate through the list of dataframes:
  for (i in seq_along(dataframes)) {
    # We look for the unique values:
    unique_rows[[i]] <- dataframes[[i]] |>
      filter(!.data[[column]] %in% unlist(sapply(dataframes[-i], `[[`, column)))

    # And then the shared values:
    common_rows[[i]] <- dataframes[[i]] |>
      filter(.data[[column]] %in% Reduce(intersect, lapply(dataframes, `[[`, column)))
  }

  # Finally, we return a list with the unique and common rows:
  return(list(unique_rows = unique_rows, common_rows = common_rows))
}
