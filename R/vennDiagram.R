#' Venn diagrams
#'
#' This function creates a Venn Diagram to compare two top-tables. It also
#' generates a table with the shared biomarkers between the two top tables,
#' as well as two tables for the unique biomarkers of each top-tables. IT IS
#' ESSENTIAL TO RUN THE CHUNK WITH 'results = "asis"' IN ORDER TO VIEW THE
#' KABLE TABLES. It has nine arguments: the two top-tables that are going to
#' be compared, the name of the estimate, adjusted P-value and biomarker
#' columns and, finally, the desired estimate and adjusted P-value thresholds.
#'
#'
#' @param endtable1 The first top-table.
#' @param endtable2 The second top-table.
#' @param estimate1 Name of the "estimate" column in the first top-table.
#' @param estimate2 Name of the "estimate" column in the second top-table.
#' @param adj1 Name of the "adjusted p-value" column in the first top-table.
#' @param adj2 Name of the "adjusted p-value" column in the second top-table.
#' @param biomarker Name of the "biomarker" column in both top-tables.
#' @param estimateVal Absolute estimate value threshold for up or down-regulation.
#' @param adjVal Adjusted p-value threshold.
#' @importFrom ggvenn ggvenn
#' @importFrom dplyr filter
#' @import ggplot2
#' @import kableExtra
#' @export


vennDiagram <- function(endtable1, endtable2, estimate1, estimate2, adj1,
                        adj2, biomarker, estimateVal, adjVal) {
  # We use dplyr's filter to extract the upregulated metabolites for each visit:
  up1 <- endtable1 %>% filter(endtable1[estimate1] > estimateVal &
                                 endtable1[adj1] < adjVal)
  up2 <- endtable2 %>% filter(endtable2[estimate2] > estimateVal &
                                 endtable2[adj2] < adjVal)

  # We do the same for the downregulated metabolites:
  down1 <- endtable1 %>% filter(endtable1[estimate1] < -estimateVal &
                                 endtable1[adj1] < adjVal)
  down2 <- endtable2 %>% filter(endtable2[estimate2] < -estimateVal &
                                 endtable2[adj2] < adjVal)

  # And we create the diagrams. First, for upregulated metabolites:
  vlistup <- list(first_up = up1[[biomarker]], second_up = up2[[biomarker]])

  # Generation of Venn Diagram for upregulated metabolites:
  p_up <- ggvenn(vlistup,
               fill_color = c("#0073C2FF", "#EFC000FF"),
               stroke_size = 0.5,
               set_name_color = "black",
               set_name_size = 4,
               text_color = "black",
               text_size = 4,
               show_percentage = FALSE
  ) + ggtitle("Venn Diagram for upregulated metabolites")

  print(p_up)

  # We store the different and shared up-regulated metabolites in new variables:
  dif1 <- setdiff(up1[[biomarker]], up2[[biomarker]])
  dif2 <- setdiff(up2[[biomarker]], up1[[biomarker]])
  common <- intersect(up1[[biomarker]], up2[[biomarker]])


  # And we create the tables:
  venn_up_1 <- endtable1[dif1, ]
  cat("\n \n")
  print(kable(venn_up_1, row.names = F,
        caption = "Upregulated metabolites of first top-table") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center"))
  cat("\n \n")

  venn_up_2 <- endtable2[dif2, ]
  cat("\n \n")
  print(kable(venn_up_2, row.names = F,
        caption = "Upregulated metabolites of second top-table") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center"))
  cat("\n \n")

  venn_up_com <- endtable1[common, ]
  cat("\n \n")
  print(kable(venn_up_com, row.names = F,
        caption = "Common upregulated metabolites between the two top-tables") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center"))
  cat("\n \n")

  # Now, for downregulated metabolites:
  vlistdw <- list(first_dw = down1[[biomarker]], second_dw = down2[[biomarker]])

  # Generation of Venn Diagram for downregulated metabolites:
  p_dw <- ggvenn(vlistdw,
                 fill_color = c("#0073C2FF", "#EFC000FF"),
                 stroke_size = 0.5,
                 set_name_color = "black",
                 set_name_size = 4,
                 text_color = "black",
                 text_size = 4,
                 show_percentage = FALSE
  ) + ggtitle("Venn Diagram for downregulated metabolites")

  print(p_dw)

  # We store the different and shared down-regulated metabolites in new variables:
  dif1 <- setdiff(down1[[biomarker]], down2[[biomarker]])
  dif2 <- setdiff(down2[[biomarker]], down1[[biomarker]])
  common <- intersect(down1[[biomarker]], down2[[biomarker]])

  # And we create the tables:
  venn_dw_1 <- endtable1[dif1, ]
  cat("\n \n")
  print(kable(venn_dw_1, row.names = F,
        caption = "Downregulated metabolites of first top-table") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center"))
  cat("\n \n")

  venn_up_2 <- endtable2[dif2, ]
  cat("\n \n")
  print(kable(venn_up_2, row.names = F,
        caption = "Downregulated metabolites of second top-table") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center"))
  cat("\n \n")

  venn_up_com <- endtable1[common, ]
  cat("\n \n")
  print(kable(venn_up_com, row.names = F,
        caption = "Common downregulated metabolites between the two top-tables") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center"))
  cat("\n \n")

}
