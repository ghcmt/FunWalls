#' Correlation Maps
#'
#' This function generates a customized correlation plot from a provided
#' dataframe. It will extract the numeric variables from the dataframe and it
#' will return a plot, highlighting positive correlations in red and negative
#' correlations in blue. In addition, we can provide a desired significance
#' threshold and it will only color those correlations that are under the
#' targeted level of significance. It only needs two arguments: a dataframe and
#' a significance threshold (NULL by default).
#'
#' @param df A provided dataframe.
#' @param sig Significance threshold (e.g., 0.05). Null by default.
#' @importFrom Hmisc rcorr
#' @importFrom dplyr select
#' @import corrplot
#' @return A correlation plot.
#' @export


corrMap <- function(df, sig = NULL) {
  # First, we only select the numerical variables, as the function won't
  # work with categorical variables:
  corrdf <- df |> select(where(is.numeric))

  # Now, we create the correlation matrix:
  corr <- rcorr(as.matrix(corrdf))

  # We calculate the signifance matrix with a set confidence level of 0.95:
  p.mat <- cor.mtest(corrdf, conf.level = 0.95)

  # And we plot the correlation matrix with our desired parameters:
  p <- corrplot(corr$r, p.mat = p.mat$p, col = rev(col(50)), cl.pos = "b",
                method = "color", number.cex = 0.8, diag = F,
                insig = "blank", tl.col = "black", order = "hclust",
                sig.level = sig, addrect = 1, rect.lwd = 1)$corrPos -> p1
  text(p1$x, p1$y, round(p1$corr, 2), cex = 0.75)

  # Finally, we return the plot:
  return(p)
}

