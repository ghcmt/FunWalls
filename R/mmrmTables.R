#' MMRM Table
#'
#' This function will create a Kable-formatted table with the most relevant
#' biomarkers in a given model. It has four arguments: a model that has to be
#' provided; a dataframe in which we have the data that we used to run the model;
#' the range of columns in which we have the biomarker information; and an 'excel'
#' argument that expects the name of the desired file and which is null by default
#' (i.e., it will not create an Excel file with the information of the table
#' if we don't specify it).
#'
#' @param model MMRM model that we have generated before
#' @param df The dataframe that we used to create the model.
#' @param bioCols Range of columns that contain the biomarkers (e.g., 5:500).
#' @param excel Name of the Excel file which will contain the table data. Default NULL.
#' @return A Kable-formatted table with the relevant information of the model.
#' @importFrom dplyr arrange
#' @importFrom tidyr pivot_wider
#' @import kableExtra writexl
#' @export


mmrmTable <- function(model, df, bioCols, visits, excel = NULL) {
  # First, we create a dataframe with the provided model:
  data2table <- as.data.frame(do.call(rbind, model))

  # We add the biomarkers into a new column:
  data2table$Biomarker <- colnames(df)[bioCols]

  # Then, we do the same with the visit:
  data2table$Covar <- rep(visits, each = length(bioCols))

  View(data2table)

  # And pivot from long to wide:
  end.table <- pivot_wider(data = data2table,
                           id_cols = Biomarker,
                           names_from = Covar,
                           values_from = c("Estimate", "Pr...t.."))

  # Next, we adjust the p-value with the Benjamini-Hochberg Procedure:
  end.table$Adjp <- p.adjust(end.table[2], method = "BH")

  # Now we can change the name of the columns:
  colnames(end.table) <- c("Biomarker", "Estimate", "p.Value", "Adj.p.Value")
  end.table <- end.table[, c("Biomarker", "Estimate", "p.Value", "Adj.p.Value")]

  # We save the table in a new dataframe:
  end.table <- as.data.frame(end.table)
  rownames(end.table) <- end.table$Biomarker

  # We order by the adjusted p-value:
  end.table <- end.table %>% arrange(Adj.p.Value)

  # And we save to a new Excel file if desired:
  if (!is.null(excel)) {
    write_xlsx(end.table, path = excel)
  }

  # Next, we calculate the number of rows with Adjusted and Raw p-values under
  # the standard significance threshold of 0.05:
  adjRows <- sum(end.table$Adj.p.Value < 0.05)
  pvalRows <- sum(end.table$Adj.p.Value < 0.05)

  # We color the rows with significant differences. We use if-else to determine
  # the number of rows to be colored, due to the fact that if we don't account
  # for this number, we might encounter a subscript error when we try to
  # show only the first rows of the table:
  if(adjRows > 1 & adjRows < 10) {
    color.me <- which(end.table$Adj.p.Value < 0.05)[1:adjRows]
  } else {
    color.me <- which(end.table$Adj.p.Value < 0.05)[1:10]
  }

  # We do the same thing with the raw p-values:
  if (pvalRows > 1 & pvalRows < 10) {
    color.me2 <- which(end.table$p.Value < 0.05)[1:pvalRows]
  } else {
    color.me2 <- which(end.table$p.Value < 0.05)[1:10]
  }

    # And then we format the table depending on the number of significant values:
  if (adjRows != 0 & pvalRows != 0) {
    tab <- kable(end.table[1:10, ], row.names = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") %>%
      row_spec(color.me2, bold = TRUE, color = "orange") %>%
      row_spec(color.me, bold = TRUE, color = "red") %>%
      footnote(general = "Ordered by the adjusted p-value.")
  } else if (adjRows == 0 & pvalRows != 0) {
    tab <- kable(end.table[1:10, ], row.names = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") %>%
      row_spec(color.me2, bold = TRUE, color = "orange") %>%
      footnote(general = "Ordered by the adjusted p-value.")
  } else {
    tab <- kable(end.table[1:10, ], row.names = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") %>%
      footnote(general = "Ordered by the adjusted p-value.")
  }

  # Finally, we return the Kable-formatted table with the relevant information
  # of the model:
  return(tab)


}

#firstmb[1]
#lastmb[1]
#mmrmTable(model = run.model, df = dfpos, visits = c("V2", "V3"), bioCols = firstmb[1]:lastmb[1])
