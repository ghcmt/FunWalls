#' Bar Plot of Metabolite Class by Significance
#'
#' This functions is a variation of metTables. Instead of a table with
#' the different N and percentatge of the biomarker Classes, it returns
#' a bar plot with that information: the percentage at the Y axis and the
#' count (n) at a given significance level specified inside the bar.
#' It has five arguments: the toptable, the names of the 'class', 'estimate'
#' and 'adjusted p-value' columns in that toptable and the significance
#' thresholds desired by the user. It returns a bar plot with the percentage of
#' biomarkers up-regulated or down-regulated depending on the significance level.
#' (if not, it returns a brief message communicating the lack of biomarkers at
#' those levels).
#'
#' @param toptable The top table with the necessary information.
#' @param class Name of the 'class' column in the top table.
#' @param estimate Name of the 'estimate' column in the top table.
#' @param adj Name of the 'adjusted p-value' column in the top table.
#' @param sigs Vector with the desired adjusted p-values thresholds.
#' @return A bar plot with the percentage of biomarkers up-regulated or down-
#' regulated depending on the significance level.
#' @importFrom tidyr pivot_wider
#' @import kableExtra
#' @rawNamespace import(dplyr, except = group_rows)
#' @export



metBarPlot <- function(toptable, class, estimate, adj, sigs) {
  # We create a dataframe to store all the information:
  fdf <- data.frame()

  # We calculate the number of biomarkers for each class present in the toptable:
  total <- toptable |>
    count(.data[[class]]) |>
    rename(total_count = n)

  # We iterate through each significance level to extract the data:
  for (sig in sigs) {
    # First, we get upregulated biomarkers:
    up <- toptable |>
      filter(toptable[estimate] > 0 & toptable[adj] < sig) |>
      group_by(.data[[class]]) |>
      summarize(significance = sig,
                diff = "Upregulated",
                n = n())

    # Next, we proceed with downregulated biomarkers:
    dw <- toptable |>
      filter(toptable[estimate] < 0 & toptable[adj] < sig) |>
      group_by(.data[[class]]) |>
      summarize(significance = sig,
                diff = "Downregulated",
                n = n())

    # We fuse the upregulated and downregulated rows:
    df <- bind_rows(up, dw)

    # And we calculate the percentage of each class with respect to the
    # total number by class that we got before:
    df <- df |>
      left_join(total, by = class) |>
      mutate(percent = n / total_count * 100, .before = total_count)

    # Then, we format the new column:
    df$percent <- as.numeric(format(df$percent, digits = 2))

    # And we fuse the new data with the existing data:
    fdf <- rbind(fdf, df)

  }

  # Now, if we have some biomarkers differentially expressed at the provided
  # significance thresholds, we start to prepare the information for the
  # final table:
  if(nrow(fdf) != 0) {
    # First, we get how many significance levels are in our final dataframe:
    sigdf <- unique(fdf$significance)

    # We do the same thing with the "diff" column to see if there are
    # downregulated and/or upregulated biomarkers:
    regdf <- unique(fdf$diff)

    # We reorder the vector to make sure that the order is the same than
    # in the final dataframe:
    regdf <- sort(regdf)

    # Now, we use complete to add rows with n = 0 if there are not upregulated
    # or downregulated biomarkers at a given significance level:
    fdf <- fdf %>%
      group_by(significance, Class, total_count) %>%
      complete(diff = distinct(., diff)$diff, fill = list(n = 0)) |>
      mutate(percent = ifelse(is.na(percent), 0, percent))

    # Next, we assess how many significance levels there are in the table to
    # label the plot accordingly:
    unique_levels <- unique(fdf$significance)
    sig_labels <- paste0("Adjusted P value = ", c(unique_levels))
    names(sig_labels) <- unique_levels

    # And we generate the plot:
    p <- ggplot(fdf, aes(x = Class, y = percent, fill = diff)) +
      geom_bar(stat="identity", size = 0.5, color = "black") +
      ylab("Percentage (%)") +
      geom_text(aes(label = ifelse(n > 0, n, NA)), position = position_stack(vjust = 0.5), size = 3.3) +
      ggtitle("Percentage (and n) of biomarker classes depending on significance") +
      theme(axis.text.x = element_text(angle = 45,  size = 12, hjust = 1, vjust = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "bottom", legend.title = element_blank()) +
      scale_fill_manual(values = c("#00b8a9", "#ff6b6b")) +
      facet_wrap(~significance, nrow = 2, labeller = as_labeller(sig_labels))

    # Finally, we return the plot:
    return(p)

    # If there are no rows in the dataframe, we return a message:
  } else {
    return(cat("There are no metabolites up-regulated or down-regulated at the
               provided significance values"))
  }

}
