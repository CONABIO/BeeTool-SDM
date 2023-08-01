#' BeeTool SDM script

# Libraries ----
library(terra)

library(fs)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(fuzzySim)
library(ENMeval)

library(magrittr)


box::use(./sdm/utils)

# Environment config ----
set.seed(1049)

# Load config
config <- config::get()


# Preprocessing ----
## Create base output dir ----
base_output <- config$base_output_path %>%
  as_fs_path()

if(!dir_exists(base_output)) {
  dir_create(base_output)
}

## Read args ----
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please enter a single parameter (input file).\n", call. = FALSE)
} else if (length(args) == 1) {
  cat("Processing model for file", args[1], "\n")
} else {
  stop("Single parameter is needed (input file).\n", call. = FALSE)
}

clean_occ_dir <- args[1]
# clean_occ_dir <- "./output/BOMIMP"

## Loading preprocessed files
crs_wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

sampled_occs_data <- path_join(c(clean_occ_dir, "data_clean.csv")) %>%
  read_csv() %>%
  as.data.frame() %>%
  vect(geom=c("x", "y"), crs=crs_wgs84)
small_sample <- if(nrow(sampled_occs_data) > 10) FALSE else TRUE

if (config$regional$use_regional_cutoff) {
  regions_of_interest <- path_join(c(clean_occ_dir, "region_of_interest.shp")) %>%
    vect()
}

## Load covars ----
if (config$covariables$is_worldclim) {
  covar_file_list <- str_c(
    config$covariables$path,
    "/wc2.1_2.5m_",
    config$covariables$variables,
    ".tif")
}

covar_rasters <- rast(covar_file_list)
if (!is.null(regions_of_interest)) {
  covar_rasters <- crop(covar_rasters, regions_of_interest, mask=TRUE)
}

sampled_occs_data_covar <- terra::extract(
  covar_rasters,
  sampled_occs_data,
  bind=TRUE)

## Tal vez sea necesario que se quiten los puntos que tengan valores NA

utils$write_points(
  sampled_occs_data_covar,
  path_join(c(clean_occ_dir, "data_clean_covar.csv"))
)

# Variable selection ----
# Add presence variable
presence_col <- "presence"
sampled_occs_data_covar[,presence_col] <- 1

covar_names <- names(covar_rasters)

covar_selection <- corSelect(
  data = sampled_occs_data_covar,
  sp.cols = presence_col,
  var.cols = covar_names
)

selected_vars <- covar_selection$selected.vars
utils$create_report(covar_selection, path_join(c(clean_occ_dir,"selection_variables_report.md")))

selected_covar_rasters <- covar_rasters[[selected_vars]]

# Create train/test set ----
sampled_row <- sample.int(
  nrow(sampled_occs_data_covar),
  size = floor(0.7*nrow(sampled_occs_data_covar))
)

sampled_occs_data_covar$train = FALSE
sampled_occs_data_covar[sampled_row, 'train'] = TRUE

utils$write_points(
  sampled_occs_data_covar,
  path_join(c(clean_occ_dir, "data_clean_covar.csv"))
)

#Pseudo-absent data
N_bg_points <- if(small_sample) 500 else 5000

background_points <- spatSample(
  selected_covar_rasters,
  N_bg_points,
  "random",
  na.rm=TRUE,
  as.points=TRUE
)

sampled_row <- sample.int(
  nrow(background_points),
  size = floor(0.7*nrow(background_points))
)

background_points$train = FALSE
background_points[sampled_row, 'train'] = TRUE

utils$write_points(
  background_points,
  path_join(c(clean_occ_dir, "background_data.csv"))
)

if (small_sample) {
  occs_train <- as.data.frame(sampled_occs_data_covar, geom="XY") %>%
    select(x,y) %>%
    rename("lon"=x, "lat"=y)

  bg_train <- as.data.frame(background_points, geom="XY") %>%
    select(x,y) %>%
    rename("lon"=x, "lat"=y)
} else {
  idx_occs_train <- which(sampled_occs_data_covar$train==TRUE)
  occs_train <- sampled_occs_data_covar[idx_occs_train, c(selected_vars, "presence")]

  occs_train <- as.data.frame(occs_train, geom="XY") %>%
    select(x,y) %>%
    rename("lon"=x, "lat"=y)


  idx_bg_train <- which(background_points$train == TRUE)
  bg_train <- background_points[idx_bg_train, selected_vars]
  bg_train <- as.data.frame(bg_train, geom="XY") %>%
    select(x,y) %>%
    rename("lon"=x, "lat"=y)
}

env <- raster::stack(selected_covar_rasters)

tune_args <- list(
  fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
  rm = seq(0.5, 4, 0.5)
)

# ENMeval ----
if (small_sample) {
  sp_models <- ENMevaluate(occs_train, env, bg_train, tune.args = tune_args,
                           partitions = "jackknife",
                           bin.output = TRUE,
                           parallel = TRUE, numCores = parallel::detectCores()-2,
                           algorithm="maxent.jar")
} else {
  sp_models <- ENMevaluate(occs_train, env, bg_train, tune.args = tune_args,
                           partitions = "randomkfold", partition.settings	= list(kfolds=4),
                           bin.output = TRUE,
                           parallel = TRUE, numCores = parallel::detectCores()-2,
                           algorithm="maxent.jar")
}

saveRDS(sp_models@models,file=path_join(c(clean_occ_dir, "Maxent_models.Rds")))
resultados_enmeval <- sp_models@results
write_csv(resultados_enmeval,
          file=path_join(c(clean_occ_dir, "enmeval_results.csv")))

model_bestAICc <- resultados_enmeval %>%
  filter(delta.AICc == 0) %>%
  pull(tune.args) %>%
  as.character()

var_imporance_best <-  eval.variable.importance(sp_models)[model_bestAICc]
# TODO: save importance

map(
  model_bestAICc,
  utils$save_raster_with_settings,
  predictions=sp_models@predictions,
  output_path=clean_occ_dir,
  prefix="ENM_prediction_M_raw_"
)

# Binary predictions ----
# thresholds considered
probs <- c(0, 0.05, 0.1)

selected_predictions <- rast(sp_models@predictions[[model_bestAICc]])

cuts_data_frame <- selected_predictions %>%
  terra::extract(sampled_occs_data) %>%
  select(any_of(model_bestAICc)) %>%
  map(quantile, probs=probs, na.rm=TRUE) %>%
  map(\(ts) as.data.frame(t(as.matrix(ts)))) %>%
  list_rbind(names_to = "tune_setting") %>%
  janitor::clean_names() %>%
  tidyr::pivot_longer(!tune_setting, names_to='threshold_name') %>%
  rename(min_presence_value=value)

output_binary_path <- path_join(c(clean_occ_dir, "binary_output"))

if(!dir_exists(output_binary_path)) {
  dir_create(output_binary_path)
}

pmap(cuts_data_frame,
     utils$save_binary_prediction,
     predictions=sp_models@predictions,
     output_path=output_binary_path,
     .progress=TRUE)
