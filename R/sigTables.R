#' Table of Significant Biomarkers
#'
#' This function returns a table with the number of significant biomarkers of a
#' provided top table. It accepts as many levels of significance as desired
#' by the user. In addition, it can work with both adjusted and raw p-values
#' as long as the user provides the name of their columns in the top table.
#' It has five arguments: the name of the toptable, the desired significance
#' levels and the names of the estimate, raw p-value and adjusted p-value columns.
#' It returns a kable-formatted table with the desired information.
#'
#' @param toptable The top table with the necessary information.
#' @param sigs Vector with the desired significance levels (e.g., c(0.01, 0.05)).
#' @param estimate Name of the 'estimate' column in the top table.
#' @param pval Name of the 'raw p-value' column in the top table.
#' @param adj Name of the 'adjusted p-value' column in the top table.
#' @param visit0 The name of the first visit (which contains baseline info) (e.g., "Visit 0").
#' @param bioCols Range of columns that contain the biomarkers (e.g., 5:500).
#' @param PDF Set to TRUE if the Tables are going to be exported to a PDF (FALSE by default).
#' @return A kable-formatted table with the number of significant biomarkers by sig. level.
#' @importFrom dplyr summarize
#' @import kableExtra
#' @export

sigTables <- function(toptable, sigs, estimate, pval, adj, PDF = FALSE) {
  # We create a dataframe to store all the information:
  fdf <- data.frame()

  # We iterate through each significance level to extract the data:
  for (sig in sigs) {
    df <- toptable %>%
      summarize(significance = sig,
                up = sum(toptable[estimate] > 0
                         & toptable[pval] < sig, na.rm = TRUE),
                upAdj = sum(toptable[estimate] > 0
                         & toptable[adj] < sig, na.rm = TRUE),
                dw = sum(toptable[estimate] < 0
                         & toptable[pval] < sig, na.rm = TRUE),
                dwAdj = sum(toptable[estimate] < 0
                         & toptable[adj] < sig, na.rm = TRUE))

    # We transform the tibble to dataframe to operate with the rows:
    df <- as.data.frame(df)

    # Now we fuse the new data with the existing data:
    fdf <- rbind(fdf, df)

  }

  # We format the table depending on the "PDF" argument:
  if (!isTRUE(PDF)){
    sigTable <- kable(fdf, row.names = F,
                      col.names = c("Significance level", "P.value",
                                    "Adj. P. Value", "P.value", "Adj. P. Value"),
                      caption = "Significantly changed metabolites") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") %>%
      add_header_above(c("", "Up-regulated" = 2, "Down-regulated" = 2))
  } else {
    sigTable <- kable(fdf, row.names = F,
                      col.names = c("Significance level", "P.value",
                                    "Adj. P. Value", "P.value", "Adj. P. Value"),
                      caption = "$Significantly changed metabolites$") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") %>%
      add_header_above(c("", "Up-regulated" = 2, "Down-regulated" = 2))
  }


  # And we return the final dataframe:
  return(sigTable)
}
