library(terra)

library(readr)
library(dplyr)
library(magrittr)

options(java.parameters = "-Xmx8000m")

# Load config
config <- config::get()

# sdm_info_data <- "./entrega/sdm_models_info.csv" %>%
#   read_csv()

future_conditions <- config %$%
  future_projection_config %$%
  rasters %>%
  purrr::map(rast)
print(length(future_conditions))

model_names <- config %$%
  future_projection_config %$%
  rasters %>%
  stringr::str_split('_') %>%
  purrr::map_chr(\(x) paste(x[5], x[6], sep = "-"))

# aoi <- config %$%
#   future_projection_config %$%
#   aoi %>%
#   vect

## Read args ----
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please enter a single parameter (input folder).\n", call. = FALSE)
} else if (length(args) == 1 && fs::dir_exists(args[1])) {
  model_output_folder <- args[1]
  cat("Projecting model for folder ", model_output_folder, "\n")
} else {
  stop("Single parameter is needed (input folder).\n", call. = FALSE)
}

output_folder <- fs::path_join(c(model_output_folder, 'prediction'))
if (!fs::dir_exists(output_folder)) {
  fs::dir_create(output_folder)
}

# Binary models processing ----
# BASEPATH <- fs::as_fs_path("./entrega/")
RESULT_CSV <- "enmeval_results.csv"
MODEL_FILE <- "Maxent_models.Rds"
AOI <- "region_of_interest.shp"

# binary_sp_smd <- sdm_info_data %>%
#   filter(!is_binary) %>%
#   pull(spCode) %>%
#   unique()
# result_path <- fs::path(BASEPATH, binary_sp_smd)

# enm_result <- fs::path(result_path[2], RESULT_CSV) %>%
#   read_csv()
# enm_model <- fs::path(result_path[2], MODEL_FILE)
# aoi <- fs::path(result_path[2], AOI) %>% vect()

# enm_result <- fs::path('entrega/BOMEPH/', RESULT_CSV) %>%
#   read_csv()
# enm_model <- fs::path('entrega/BOMEPH/', MODEL_FILE)
# aoi <- fs::path('entrega/BOMEPH/', AOI) %>% vect()

enm_result <- fs::path(model_output_folder, RESULT_CSV) %>%
  read_csv()
enm_model <- fs::path(model_output_folder, MODEL_FILE)
aoi <- fs::path(model_output_folder, AOI) %>% vect()

future_conditions <- future_conditions %>%
  purrr::map(crop, aoi)

model_bestAICc <- enm_result %>%
  filter(delta.AICc == 0) %>%
  pull(tune.args) %>%
  as.character()
fitted_models <- enm_model %>% read_rds()

maxent_model <- fitted_models[[model_bestAICc[1]]]
bio_vars <- names(maxent_model@presence)
selected_vars_id <- bio_vars %>%
  stringr::str_split('_') %>%
  purrr::map_chr(-1) %>%
  stringr::str_c("_", .) %>%
  stringr::str_c(collapse = "|")

threshold <- maxent_model@presence %>%
  dismo::predict(maxent_model, .) %>%
  stats::quantile(c(0.1))

select_masks <- future_conditions %>%
  purrr::map(names) %>%
  purrr::map(stringr::str_ends, pattern=selected_vars_id)

future_conditions <- future_conditions %>%
  purrr::map2(select_masks, \(x, y) x[[y]])

future_conditions %>%
  purrr::map(terra::set.names, bio_vars)

print("Proyectando a futuro")
projections <- future_conditions %>%
  purrr::map(\(x) dismo::predict(maxent_model, x))


# projections <- list()

# for (i in length(future_conditions)) {
#   projections[i] <- dismo::predict(maxent_model, future_conditions[[i]])
# }

for (i in seq_along(future_conditions)) {
  fn <- paste0(model_names[i],'.tif')
  print(paste("Escribiendo raster ", fn))
  terra::writeRaster(projections[[i]]>=threshold, fs::path_join(c(output_folder, fn)), overwrite=TRUE)
}

# # TODO ----
# [ ] Decidir que aoi tomar
# [ ] Salvar raster a archivo

# # Draft ----
# fitted_models <- "/home/jbarrios/Projects/BeeTool-Bee_DM/done/single_model/BOMHUN/Maxent_models.Rds" %>%
#   read_rds()
# result_models <- "/home/jbarrios/Projects/BeeTool-Bee_DM/done/single_model/BOMHUN/enmeval_results.csv" %>%
#   read_csv()
#
# # fitted_models <- "/home/jbarrios/Projects/BeeTool-Bee_DM/done/ANDALI/Maxent_models.Rds" %>%
# #   read_rds()
# # result_models <- "/home/jbarrios/Projects/BeeTool-Bee_DM/done/ANDALI/enmeval_results.csv" %>%
# #   read_csv()
#
# future_raster <- "./data/future_worldclim/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2021-2040.tif" %>%
#   rast()
#
# aoi <- "./data/aoi/norteamerica.shp" %>% vect()
#
# model_bestAICc <- result_models %>%
#   filter(delta.AICc == 0) %>%
#   pull(tune.args) %>%
#   as.character()
#
# maxent_model <- fitted_models[[model_bestAICc[1]]]
# bio_vars <- names(maxent_model@presence)
# selected_vars_id <- bio_vars %>%
#   stringr::str_split('_') %>%
#   purrr::map_chr(-1) %>%
#   stringr::str_c("_", .) %>%
#   stringr::str_c(collapse = "|")
#
# select_mask <- names(future_raster) %>%
#   stringr::str_ends(selected_vars_id)
#
# future_raster <- future_raster[[select_mask]] %>%
#   crop(aoi)
# names(future_raster) <- bio_vars
#
# future_model <- predict(maxent_model, future_raster)
