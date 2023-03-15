#' MMRM Table
#'
#' This function will create a Kable-formatted table with the most relevant
#' biomarkers in a given model for each visit. It has eight arguments: a model
#' that has to be provided; a dataframe in which we have the data that we used
#' to run the model; the range of columns in which we have the biomarker
#' information; and an 'excel'argument that expects the name of the desired file
#' and which is null by default (i.e., it will not create an Excel file with the
#' information of the table if we don't specify it). Essential to run it
#' with "results = 'asis'", as it contains Kable tables.
#'
#' @param model MMRM model that we have generated before
#' @param df The dataframe that we used to create the model.
#' @param bioCols Range of columns that contain the biomarkers (e.g., 5:500).
#' @param visits Vector with the visits to be analyzed.
#' @param excel Name of the desired Excel file for the output. Without ".xlsx".
#' @param classFile Name and path of the Excel file (with .xlsx) with class
#' information for the biomarkers.
#' @param classBio Name of the Biomarker column in the class file.
#' @param classGroup Name of the desired Group column in the class file.
#' @importFrom dplyr arrange
#' @importFrom tidyr pivot_wider
#' @import kableExtra writexl readxl
#' @return A toptable with all the relevant information.
#' @export


mmrmTable <- function(model, df, bioCols, visits, excel = NULL, classFile = NULL,
                      classBio = NULL, classGroup = NULL) {
  # Dataframe with all the selected elements:
  data2table <- as.data.frame(do.call(rbind, run.model))

  # We add the metabolites into a new column:
  data2table$Biomarker <- rep(colnames(df)[bioCols], each = length(visits))

  #Then, we do the same with the visits:
  data2table$Covar <- rep(visits, lastmb[1] - firstmb[1] + 1)

  # And pivot from long to wide:
  end.table <- pivot_wider(data = data2table,
                           id_cols = Biomarker,
                           names_from = Covar,
                           values_from = c("Estimate", "Pr...t.."))

  for(visit in visits) {
    name_raw <- paste0("Pr...t.._", visit)
    name_adj <- paste0("Adj.p.Value_", visit)
    end.table[name_adj] <- p.adjust(end.table[[name_raw]], method = "BH")
    end.table <- end.table |>
      relocate(all_of(c(name_raw, name_adj)), .after = paste0("Estimate_", visit))
  }


  # Now we can change the name of the columns:
  names(end.table) <- sub("Pr...t..", "p.Value", colnames(end.table))

  # And we change the rownames to the Biomarker:
  end.table <- as.data.frame(end.table)
  rownames(end.table) <- end.table$Biomarker

  # If we want Class information, we merge our dataframe with the Class dataframe:
  if(!is.null(classFile)) {
    # We read the file:
    bioclass <- read_xlsx(path = classFile)

    # We get the relevant information:
    bioclass <- bioclass |>
      rename(Biomarker = all_of(classBio),
             Class = all_of(classGroup)) |>
      select(Biomarker, Class)

    # And we merge the two dataframes:
    end.table <- left_join(end.table, bioclass)
  }

  # Finally, we generate the tables for each Visit:
  for(visit in visits) {
    name_raw <- paste0("p.Value_", visit)
    name_adj <- paste0("Adj.p.Value_", visit)

    # We arrange the table to order it by ascending adjusted p-value for that visit:
    end.table <- end.table |>
      arrange(.data[[name_adj]], .data[[name_raw]])


    if (!is.null(excel)) {
      # We check if the "results" folder already exists; if not, we create it:
      ifelse(!dir.exists("results/"), dir.create("results/"), FALSE)

      # And we store the results in there:
      write_xlsx(end.table, path = paste0("results/", excel, visit, ".xlsx"))
    }

    # Now, we prepare the data for the table:
    color_raw <- which(end.table[name_raw] < 0.05)
    color_adj <- which(end.table[name_adj] < 0.05)

    # As we will only show the first ten rows, we discard indexes beyond that:
    color_raw <- color_raw[color_raw <= 10]
    color_adj <- color_adj[color_adj <= 10]

    # And we configure the table:
    t <- kable(end.table[1:10, ], row.names = F, caption = paste0("Most significally changed metabolites after the treatment in ", visit, ".")) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") %>%
      row_spec(color_raw, bold = TRUE, color = "orange") %>%
      row_spec(color_adj, bold = TRUE, color = "red") %>%
      footnote(paste0("Ordered by ascendent adjusted p-value of ", visit, "."))
    cat("\n \n")
    cat(paste0("The following table shows the 10 most significantly changed metabolites in ", visit, ". Full table can be found in the `", excel, visit, ".xlsx` file at the `Results` folder."))
    cat("\n \n")


    # Finally, we print the final table:
    print(t)

    # And we change the rownames to the Biomarker for future analysis:
    end.table <- as.data.frame(end.table)
    rownames(end.table) <- end.table$Biomarker
  }

  return(end.table)
}

