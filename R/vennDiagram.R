#' Venn diagrams
#'
#' This function creates a Venn Diagram to compare two or more top-tables. It also
#' generates a table with the shared biomarkers between the top tables,
#' as well as tables for the unique biomarkers of each of the top-tables. IT IS
#' ESSENTIAL TO RUN THE CHUNK WITH 'results = "asis"' IN ORDER TO VIEW THE
#' KABLE TABLES. It has eight arguments: the list of top-tables that are going to
#' be compared, a vector with the names of the "estimate" and "adjusted P-value"
#' columns, the name of the biomarker column and, also, the desired estimate and
#' adjusted P-value thresholds. In addition, it also has an 'excel' argument to
#' save those tables that are too large to print. This argument has to be filled
#' without ".xlsx" and has to provide some information about the top-tables
#' (e.g., "Gender" or "MMSE"). Finally, it has a set_names argument to provide
#' the desired names for each of the top-tables (NULL by default)
#'
#'
#' @param endtables List of top-tables. If there is only one top-table, the input
#' should be "list(toptable)".
#' @param estimate Vector with the names of the "estimate" column in the top-tables.
#' There is no need for a vector if the name is shared between top-tables.
#' @param adj Vector with the names of the "adjusted p-value" column in the top-tables.
#' There is no need for a vector if the name is shared between top-tables.
#' @param biomarker Name of the "biomarker" column in the top-tables.
#' @param estimateVal Absolute estimate value threshold for up or down-regulation.
#' @param adjVal Adjusted p-value threshold.
#' @param excel Name of the Excel file to store results if tables are large (without .xlsx)
#' @param set_names Vector with the desired names for the top-tables
#' (i.e., c("Top1", "Top2", "Top3"))
#' @importFrom ggvenn ggvenn
#' @importFrom dplyr filter
#' @importFrom openxlsx write.xlsx
#' @import ggplot2
#' @import kableExtra
#' @import writexl
#' @export


vennDiagram <- function(endtables, estimate, adj, biomarker, estimateVal,
                         adjVal, excel, set_names = NULL) {
  # We use mapply to apply the filter function to each endtable and its corresponding estimate and adj columns:
  up <- mapply(function(x, e, a) filter(x, x[e] > estimateVal & x[a] < adjVal),
               endtables, estimate, adj, SIMPLIFY = FALSE)

  # We do the same for the downregulated metabolites:
  down <- mapply(function(x, e, a) filter(x, x[e] < -estimateVal & x[a] < adjVal),
                 endtables, estimate, adj, SIMPLIFY = FALSE)

  # And we create the diagrams. First, for upregulated metabolites:
  vlistup <- setNames(lapply(up, `[[`, biomarker), set_names)

  # Generation of Venn Diagram for upregulated metabolites:
  p_up <- ggvenn(vlistup,
                 fill_color = c("#0073C2FF", "#EFC000FF", "#b5050e", "#CD534CFF",
                                "#00A17EFF", "#C853BFFF", "#EAC117FF", "#A67B5BFF"),
                 stroke_size = 0.5,
                 set_name_color = "black",
                 set_name_size = 4,
                 text_color = "black",
                 text_size = 4,
                 show_percentage = FALSE) +
    ggtitle("Venn Diagram for upregulated metabolites")

  # Now we print the Diagram:
  print(p_up)

  # And we find the unique and common Biomarkers for each table by calling
  # the findUCRows function:
  results <- findUCRows(up, biomarker)

  # Then, we focus on unique upregulated Biomarkers for each table:
  for(i in 1:length(up)) {
    df <- names(vlistup[i])
    if (nrow(results$unique_rows[[i]]) > 0 & nrow(results$unique_rows[[i]]) <= 10) {
      cat("\n \n")
      print(kable(results$unique_rows[[i]], row.names = F,
                  caption = paste0("Unique upregulated metabolites in ", df, ".")) %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                            full_width = TRUE, position = "center"))
      cat("\n \n")
    } else if (nrow(results$unique_rows[[i]]) > 10) {
      write_xlsx(results$unique_rows[[i]], paste0("results/", excel, "_", df,
                                                  "_unique_upregulated.xlsx"))
      cat("\n \n")
      print(kable(results$unique_rows[[i]][1:10,], row.names = F,
                  caption = paste0("Some of the unique upregulated metabolites in ", df, ".")) %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                            full_width = TRUE, position = "center"))
      cat("\n \n")
      cat(paste0("Full table can be found in the file `", excel, "_", df,
                 "_unique_upregulated.xlsx` at the `results` folder."))
      cat("\n \n")
    } else {
      cat("\n \n")
      cat(paste0("There are no unique differentially upregulated metabolites in ", df, "."))
      cat("\n \n")
    }
  }

  # And then for common biomarkers:
  if (nrow(results$common_rows[[1]]) > 0 & nrow(results$common_rows[[1]]) <= 10) {
    cat("\n \n")
    print(kable(results$common_rows[[1]], row.names = F,
                caption ="Common upregulated metabolites") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                          full_width = TRUE, position = "center"))
    cat("\n \n")
  } else if (nrow(results$common_rows[[1]]) > 10) {
    write.xlsx(x = results$common_rows,
               file = paste0("results/", excel, "_common_upregulated.xlsx"),
               sheetName = names(vlistup))
    cat("\n \n")
    print(kable(results$common_rows[[1]][1:10,], row.names = F,
                caption = paste0("Some of the shared upregulated metabolites (values from ",
                                 names(vlistup[1]), ").")) %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                          full_width = TRUE, position = "center"))
    cat("\n \n")
    cat(paste0("Full table can be found in the file `", excel,
               "_common_upregulated.xlsx` at the `results` folder."))
    cat("\n \n")
  } else {
    cat("\n \n")
    cat("There are no shared upregulated metabolites in these tables.")
    cat("\n \n")
  }

  # And, now, for downregulated biomarkers:
  vlistdw <- setNames(lapply(down, `[[`, biomarker), set_names)

  # Generation of Venn Diagram for upregulated metabolites:
  p_dw <- ggvenn(vlistdw,
                 fill_color = c("#0073C2FF", "#EFC000FF", "#b5050e", "#CD534CFF",
                                "#00A17EFF", "#C853BFFF", "#EAC117FF", "#A67B5BFF"),
                 stroke_size = 0.5,
                 set_name_color = "black",
                 set_name_size = 4,
                 text_color = "black",
                 text_size = 4,
                 show_percentage = FALSE) +
    ggtitle("Venn Diagram for downregulated metabolites")

  # Now we print the Diagram:
  print(p_dw)

  # And we create the tables. First, for unique downregulated biomarkers:
  results <- findUCRows(down, biomarker)

  for(i in 1:length(down)) {
    df <- names(vlistdw[i])
    if (nrow(results$unique_rows[[i]]) > 0 & nrow(results$unique_rows[[i]]) <= 10) {
      cat("\n \n")
      print(kable(results$unique_rows[[i]], row.names = F,
                  caption = paste0("Unique downregulated metabolites in ", df)) %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                            full_width = TRUE, position = "center"))
      cat("\n \n")
    } else if (nrow(results$unique_rows[[i]]) > 10) {
      write_xlsx(results$unique_rows[[i]], paste0("results/", excel, "_", df,
                                                  "_unique_downregulated.xlsx"))
      cat("\n \n")
      print(kable(results$unique_rows[[i]][1:10,], row.names = F,
                  caption = paste0("Some of the unique downregulated metabolites in ", df, ".")) %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                            full_width = TRUE, position = "center"))
      cat("\n \n")
      cat(paste0("Full table can be found in the file `", excel, "_", df,
                 "_unique_downregulated.xlsx` at the `results` folder."))
      cat("\n \n")
    } else {
      cat("\n \n")
      cat(paste0("There are no unique differentially downregulated metabolites in ", df, "."))
      cat("\n \n")
    }
  }

  # And then for common downregulated biomarkers:
  if (nrow(results$common_rows[[1]]) > 0 & nrow(results$common_rows[[1]]) <= 10) {
    cat("\n \n")
    print(kable(results$common_rows[[1]], row.names = F,
                caption ="Common downregulated metabolites") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                          full_width = TRUE, position = "center"))
    cat("\n \n")
  } else if (nrow(results$common_rows[[1]]) > 10) {
    write.xlsx(x = results$common_rows,
               file = paste0("results/", excel, "_common_downregulated.xlsx"),
               sheetName = names(vlistup))
    cat("\n \n")
    print(kable(results$common_rows[[1]][1:10,], row.names = F,
                caption = paste0("Some of the shared downregulated metabolites (values from ",
                                 names(vlistdw[1]), ").")) %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                          full_width = TRUE, position = "center"))
    cat("\n \n")
    cat(paste0("Full table can be found in the file `", excel,
               "_common_downregulated.xlsx` at the `results` folder."))
    cat("\n \n")
  } else {
    cat("\n \n")
    cat("There are no shared downregulated metabolites in these tables.")
    cat("\n \n")
  }

}
