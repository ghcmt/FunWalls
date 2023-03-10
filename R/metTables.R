#' Table of Metabolite Class by Significance
#'
#' This functions returns a table with the different classes of biomarkers
#' in a given top-table and the number of total counts and the counts for
#' each significance level (provided by the user). It also calculates the
#' percentages of counts for each significance threshold with respect to the
#' total count for that class and divides the biomarkers in upregulated
#' and downregulated. Thus, it provides a general view of which classes
#' are being upregulated or downregulated. It has five arguments: the toptable,
#' the names of the 'class', 'estimate' and 'adjusted p-value' columns in that
#' toptable and the significance thresholds desired by the user. It returns a
#' table with all the information if there are biomarkers upregulated or down-
#' regulated below the provided significance thresholds (if not, it returns
#' a brief message communicating the lack of biomarkers at those levels).
#'
#' @param toptable The top table with the necessary information.
#' @param class Name of the 'class' column in the top table.
#' @param estimate Name of the 'estimate' column in the top table.
#' @param adj Name of the 'adjusted p-value' column in the top table.
#' @param sigs Vector with the desired adjusted p-values thresholds.
#' @return A kable-formatted table with the number of differentially expressed biomarkers
#' by sig. level and class.
#' @importFrom tidyr pivot_wider
#' @import kableExtra
#' @rawNamespace import(dplyr, except = group_rows)
#' @export


metTables <- function(toptable, class, estimate, adj, sigs) {
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

    # Then, we use pivot_wider to ease the understanding of the data:
    fdf <- fdf |>
      pivot_wider(names_from = c(diff, significance),
                  values_from = n:percent,
                  values_fill = 0,
                  names_vary = "slowest",
                  names_glue = "{.value} {diff} at sig = {significance}")

    # We add the percentage symbol to the percentages columns:
    fdf <- fdf |> mutate(across(contains("percent"), ~paste0(., "%")))

    # And we prepare the headers for the table:
    regheader <- c("", "", rep(setNames(rep(2, length(regdf)), regdf), length(sigdf)))
    sigheader <- c("", "", setNames(rep(2*length(regdf), length(sigdf)), paste("Adj. p-Val =", sigdf)))


    # Finally, we create a vector with the names of the columns:
    colstab <- c("Class", "Total Count", rep(c("n", "Percentage"), length(sigdf) * length(regdf)))

    # Now, we format the table:
    metTable <- kable(fdf, row.names = F, col.names = colstab) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") |>
      add_header_above(regheader) |>
      add_header_above(sigheader)

    # And we return it:
    return(metTable)

    # If there are no rows in our dataframe, we return a message explaining the
    # situation
  } else {
    return(cat("There are no metabolites up-regulated or down-regulated at the provided significance values"))
  }

}
