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
  coords <- terra::crds(x)
  df <- dplyr::bind_cols(
    as.data.frame(x),
    coords)

  readr::write_csv(df, file, na=na, ...)
}