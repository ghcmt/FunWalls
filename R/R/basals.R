#' Baseline Levels
#'
#' This function returns a dataframe with the baseline levels of the selected
#' biomarkers, adding new columns with those baseline levels with the name of
#' the metabolite and the suffix "_basal". It has five arguments: the provided
#' dataframe, the name of the column in which the patient information is stored
#' (in order to group by this variable, as each patient has different baseline
#' levels for each biomarker), the column in which the visit information is saved,
#' the name of the first visit (e.g., "FPE1 ANTES") and, finally, the range of the
#' biomarkers columns (e.g., if biomarkers levels begin at column 5 and finish at
#' 35, this parameter should be 5:35).
#'
#' @param df The dataframe with the expression levels of the biomarkers.
#' @param patientCol The column name with the information of the patients (e.g., "patients").
#' @param visitCol The column name with the information of the visits (e.g., "visits").
#' @param visit0 The name of the first visit (which contains baseline info) (e.g., "Visit 0").
#' @param bioCols Range of columns that contain the biomarkers (e.g., 5:500).
#' @return A dataframe with the original columns and the baseline expression columns.
#' @importFrom dplyr arrange group_by mutate
#' @export

basals <- function(df, patientCol, visitCol, visit0, bioCols) {
  # We pass to symbol the two columns provided, as it eases the evaluation
  # with dplyr:
  patientCol <- sym(patientCol)
  visitCol <- sym(visitCol)

  # Now we generate the dataframe with the baseline levels for each biomarker:
  df <- df %>%
    # We group by the patient column:
    group_by(!!patientCol) %>%
    # And we create a new column for each biomarker:
    mutate(across(
      .cols = all_of(bioCols-1),
      # If there is a visit0, we get that value; else, the patient will have
      # NA as the baseline level for that biomarker:
      ~ifelse(visit0 %in% !!visitCol, .[!!visitCol == visit0], NA),
      # We add the label "_basal" to the name of the column:
      .names = "{.col}_basal"))

  # And we return a dataframe:
  return(df)
}
