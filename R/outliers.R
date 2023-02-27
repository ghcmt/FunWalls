#' Outliers
#'
#' This function is designed to find all the outliers that are present in
#' a given dataframe. It iterates through every numeric column to get these
#' outliers values, which can be described as data points that differ
#' significantly from other observations. The function defines an outlier
#' as any value that is 1.5 times the interquantile range (IQR)  below the lower
#' quartile or 1.5IQR above the upper quartile. The function generates a boxplot
#' for each variable and highlights the outliers in red. In addition, the user
#' has the option of labelling the outliers by setting the argument label as
#' TRUE (it is FALSE by default), which individually labels each outlier with an
#' unique name and it also adds a subtitle to the plot with the outliers' values.
#' Finally, the function stores all the outliers detected in the dataframe in a
#' new dataframe, which it is also the object that it returns.
#'
#' @param df A given dataframe
#' @param label Individual label for each outlier (by default it is FALSE).
#' @return A dataframe with all the outliers found.
#' @importFrom dplyr filter mutate select
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export

outliers <- function(df, label = FALSE) {
  # First, we get the numeric variables of the dataframe:
  num_vars <- df |> select(where(is.numeric)) |> names()

  # Next, we create a dataframe to store the outliers:
  outlier_df <- data.frame()

  # And we create a custom theme for our boxplots:
  plot_theme <- theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.line = element_line(linewidth = 1),
                      panel.background = element_rect(fill = 'white'))

  # Then, we start the iteration through all the numeric variables:
  for(var in num_vars) {
    # Generation of interquantile range for each variable:
    q1 <- quantile(df[[var]], probs = 0.25, na.rm = TRUE)
    q3 <- quantile(df[[var]], probs = 0.75, na.rm = TRUE)
    iqr <- q3 - q1

    # We set the limits of non-outlier data:
    outlier_min <- q1 - 1.5 * iqr
    outlier_max <- q3 + 1.5 * iqr

    # And we proceed with all the data points that are beyond these limits:
    outlier_rows <- df |>
      filter(.[[var]] < outlier_min | .[[var]] > outlier_max) |>
      # Creation of a unique label for each outlier, so it eases its
      # later identification:
      mutate(outlier_label = paste0(var, "_outlier_", row_number()))

    # Now, we begin with the plots. First, if there are outliers and we want
    # labels:
    if(nrow(outlier_rows) > 0 & isTRUE(label)) {
      plot <- ggplot(df, aes(x = var, y = .data[[var]])) +
        geom_boxplot(outlier.color = "red") +
        geom_text_repel(data = outlier_rows, aes(label = outlier_label)) +
        ggtitle(paste("Outliers in variable", var),
                subtitle = paste("Outliers values:", paste(outlier_rows[[var]],
                                                           collapse = ", "))) +
        plot_theme

      # If there are outliers but we don't want labels:
    } else if(nrow(outlier_rows) > 0) {
      plot <- ggplot(df, aes(x = var, y = .data[[var]])) +
        geom_boxplot(outlier.color = "red") +
        ggtitle(paste("Outliers in variable", var)) +
        plot_theme

      # And if there are no outliers in the variable:
    } else {
      plot <- ggplot(df, aes(y = .data[[var]])) +
        geom_boxplot() +
        ggtitle(paste("No outliers in variable", var)) +
        plot_theme
    }

    # Finally, we add the outliers to our newly created dataframe:
    outlier_df <- rbind(outlier_df, outlier_rows)

    # And we print the boxplot for each variable:
    print(plot)
  }

  # We return a dataframe with all the outliers:
  return(outlier_df)
}
