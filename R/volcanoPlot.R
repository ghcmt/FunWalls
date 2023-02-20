#' Volcano plots
#'
#' This function will generate a Volcano plot from a provided top table. It will
#' highlight the differentially expressed biomarkers that are down or up-regulated.
#' The user can provide the estimate threshold from which a biomarker will be
#' up or down-regulated, as well as the desired adjusted p-value threshold.
#' The function has six arguments: the top-table with the data, the name
#' of the estimate column in that top-table, the name of the adjusted p-value
#' in the top-table, the name of the Biomarker column (i.e., the column that
#' contains the name or identification of each biomarker) and then the desired
#' thresholds for the estimate and adjusted p-value. It returns an interactive
#' ggplotly object that allows the user to focus on the relevant biomarkers, as the basic information of the highlighted
#' biomarkers will be displayed on mouseover.
#'
#' @param endtable Top-Table with the biomarker's data.
#' @param estimate Name of the "estimate" column in the top-table.
#' @param adj Name of the "adjusted p-value" column in the top-table.
#' @param biomarker Name of the "biomarker" column in the top-table.
#' @param estimateVal Absolute estimate value threshold for up or down-regulation.
#' @param adjVal Adjusted p-value threshold.
#' @return An interactive ggplotly object with the Volcano plot.
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr mutate
#' @export


volcanoPlot <- function(endtable, estimate, adj, biomarker, estimateVal, adjVal) {
  # We create a general label for all metabolites:
  endtable$diffexpressed <- "NO"

  # Then, we change this label if metabolites are upregulated or downregulated.
  # The threshold for change is set at 0.2. Metabolites not affected will
  # mantain the "NO" label:
  endtable$diffexpressed[endtable[estimate] >= estimateVal &
                           endtable[adj] < adjVal] <- "UP"
  endtable$diffexpressed[endtable[estimate] <= -estimateVal &
                           endtable[adj] < adjVal] <- "DOWN"

  # Next, we assign colors to the differentially expressed labels:
  mycolors <- c("blue", "red", "grey")
  names(mycolors) <- c("DOWN", "UP", "NO")

  # And we only label the differentially expressed metabolites with the name of
  # the biomarker. This would be useful for the interactive plot:
  end.table <- end.table %>% mutate(
    metabolite = case_when(
      diffexpressed != "NO" ~ Biomarker))

  # Finally, we proceed with the Volcano Plot:
  vplot <- ggplot(data = endtable, aes(x = .data[[estimate]],
                                    y = -log10(.data[[adj]]),
                                    col = diffexpressed,
                                    label = Biomarker)) +
    geom_point() +
    scale_colour_manual(values = mycolors) +
    theme_minimal() +
    ggtitle(paste("Volcano plot \nAdj.Pvalue <", adjVal,
                  "& abs(Estimate) >=", estimateVal)) +
    labs(x = "Estimate")

  # And we return the plotly object:
  return(ggplotly(vplot))
}
