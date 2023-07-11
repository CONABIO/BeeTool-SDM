#' Sample occurrences stratified by regular grid of given resolution
#'
#' @param occs a `SpatVector` object of occurrences
#' @param grid_res resolution for the regular grid in degree units
#'
#' @returns a `SpatVector` point dataset with the selected occurrences
#' @export
sample_occurrences <- function(occs, grid_res = 1) {
  aoi <- terra::rast(occs)
  terra::res(aoi) <- grid_res

  aoi <- terra::extend(aoi, terra::ext(aoi) + grid_res)

  occs_samples <- terra::spatSample(occs, size=1, "random", strata=aoi)

  return(occs_samples)
}


#' Save `SpatVector` points object to csv with x y columns
#' @param x `SpatVector` object to save.
#' @param file 	file or connection to write to.
#' @param na string used for missing values. Defaults to NA. Missing values
#'           will never be quoted; strings with the same value as na will always
#'           be quoted.
#' @param ... parameters forwarded to a `readr::write_csv` function
#' @export
write_points <- function(x, file, na = "", ...) {
  df <- as.data.frame(x, geom="XY")

  readr::write_csv(df, file, na=na, ...)
}

#' Create a report for a `fuzzySim::corSelect` output
#' @param results `corSelect` result list
#' @param output_file path to save the result report
#' @return None.
#' @export
create_report <- function(my_list, output_file) {
  # Extract elements from the list
  high_correlations <- my_list$high.correlations
  bivariate_significance <- my_list$bivariate.significance
  excluded_vars <- my_list$excluded.vars
  selected_vars <- my_list$selected.vars
  selected_var_cols <- my_list$selected.var.cols
  strongest_remaining_corr <- my_list$strongest.remaining.corr
  remaining_multicollinearity <- my_list$remaining.multicollinearity

  # Open a connection to the output file
  file_conn <- file(output_file, "w")

  # Write report content to the file in Markdown format
  cat("# Report\n\n", file = file_conn)

  cat("## High Correlations\n\n", file = file_conn)
  writeLines(knitr::kable(high_correlations, "pipe"), con = file_conn)
  cat("\n", file = file_conn)

  cat("## Bivariate Significance\n\n", file = file_conn)
  writeLines(knitr::kable(bivariate_significance, "pipe"), con = file_conn)
  cat("\n", file = file_conn)

  cat("## Excluded Variables\n\n", file = file_conn)
  writeLines(excluded_vars, con = file_conn)
  cat("\n", file = file_conn)

  cat("## Selected Variables\n\n", file = file_conn)
  writeLines(selected_vars, con = file_conn)
  cat("\n", file = file_conn)

  cat("## Selected Variable Columns\n\n", file = file_conn)
  writeLines(as.character(selected_var_cols), con = file_conn)
  cat("\n", file = file_conn)

  cat("## Strongest Remaining Correlation\n\n", file = file_conn)
  cat(strongest_remaining_corr, file = file_conn)
  cat("\n\n", file = file_conn)

  cat("## Remaining Multicollinearity\n\n", file = file_conn)
  writeLines(knitr::kable(remaining_multicollinearity, "pipe"), con = file_conn)
  cat("\n", file = file_conn)

  # Close the file connection
  close(file_conn)
}
