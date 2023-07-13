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

## Load data ----
if (config$regional$use_regional_cutoff) {
  regional_data <- vect(config$regional$shapefile_region_path)
}

## Read args ----
# args = commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
#   stop("Please enter a single parameter (input file).\n", call. = FALSE)
# } else if (length(args) == 1) {
#   print(paste("Processing model for file ", args[1]))
# } else {
#   stop("Single parameter is needed (input file).\n", call. = FALSE)
# }
#
# inputDataFile <- args[1]

input_data_file <- "./data/BOMHUN.csv"

output_folder <- input_data_file %>%
  path_file() %>%
  path_ext_remove()

if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

### Load occurrence data ----
crs_wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
occs_data <- read_csv(input_data_file)

# Coordinates columns are X and Y
occs_data <- vect(
  as.matrix(select(occs_data, X, Y)),
  crs=crs_wgs84,
  atts=as.data.frame(select(occs_data, -X, -Y))
  )

### Create records datasets ----
# Sample records
sampled_occs_data <- utils$sample_occurrences(
  occs_data,
  grid_res = config$spatial_resolution)

utils$write_points(
  sampled_occs_data,
  path_join(c(output_folder, "data_clean.csv"))
  )

## Select regional data of interest ----
regions_of_interest <- NULL
if (config$regional$use_regional_cutoff) {
  regions_of_interest <- extract(regional_data, sampled_occs_data)
  writeVector(
    regions_of_interest,
    path_join(c(output_folder, "region_of_interest.shp")),
    overwrite=TRUE)
}

## Load covariables ----
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
  path_join(c(output_folder, "data_clean_covar.csv"))
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
utils$create_report(covar_selection, "selection_variables_report.md")

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
  path_join(c(output_folder, "data_clean_covar.csv"))
  )

#Pseudo-absent data

background_points <- spatSample(
  selected_covar_rasters,
  5000,
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
  path_join(c(output_folder, "background_data.csv"))
)

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

env <- raster::stack(selected_covar_rasters)

tune_args <- list(
  fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
  rm = seq(0.5, 4, 0.5)
)

# ENMeval ----
sp_models <- ENMevaluate(occs_train, env, bg_train, tune.args = tune_args,
                         partitions = "randomkfold", partition.settings	= list(kfolds=4),
                         bin.output = TRUE,
                         parallel = TRUE, numCores = parallel::detectCores()-2,
                         algorithm="maxent.jar")


saveRDS(sp_models@models,file=path_join(c(output_folder, "Maxent_models.Rds")))
resultados_enmeval <- sp_models@results
write_csv(resultados_enmeval,
          file=path_join(c(output_folder, "enmeval_results.csv")))

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
  output_path=output_folder,
  prefix="ENM_prediction_M_raw_"
)

# Binary predictions ----
# thresholds considered
probs <- c(0, 0.05, 0.1)

selected_predictions <- rast(sp_models@predictions[[model_bestAICc]])

cuts_data_frame <- selected_predictions %>%
  terra::extract(sampled_occs_data) %>%
  select(any_of(model_bestAICc)) %>%
  map(quantile, probs=probs) %>%
  map(\(ts) as.data.frame(t(as.matrix(ts)))) %>%
  list_rbind(names_to = "tune_settings") %>%
  janitor::clean_names() %>%
  tidyr::pivot_longer(!tune_settings, names_to='level')

output_binary_path <- path_join(c(output_folder, binary))

if(!dir_exists(output_binary_path)) {
  dir_create(output_binary_path)
}



# old ----
# dir.create(file.path(outputFolder, "Outputs_todos"))
# writeRaster(sp.models_p, file = file.path(outputFolder, "Outputs_todos/", paste0(outputFolder)),
#             suffix='names',
#             format = "GTiff",
#             bylayer=TRUE,
#             overwrite= TRUE)
#
# delta_aic <- which(resultados_enmeval$delta.AICc == 0)
# modelsAIC0 <- resultados_enmeval %>%
#   mutate(index = rownames(resultados_enmeval)) %>%
#   filter(delta.AICc == 0) %>%
#   select(index, settings) %>%
#   mutate(index = as.numeric(index), settings = as.character(settings))
#
# aic.opt <- sp.models@models[[which(sp.models@results$delta.AICc==0)]]
# importa <- var.importance(aic.opt)
# write.csv( importa,
#            file = file.path(outputFolder, "varImportance.csv"),
#            row.names = FALSE)
#
#
# ####ENMTest####
# #source("funciones_LAE.R")
# #Threslhold independent
#
# #AUC
# aucCalculator <- function(prediction, occs, bgPoints) {
#   data <- rbind(occs, setNames(bgPoints, names(occs)))
#   labels <- c(rep(1, nrow(occs)),
#               rep(0, nrow(bgPoints)))
#   scores <- raster::extract(prediction, data)
#   pred <- ROCR::prediction(scores, labels)
#   # perf <- performance(pred, "tpr", "fpr")
#   auc <- performance(pred, "auc")@y.values[[1]]
#   return(auc)
# }
#
# aucStatistcs <- function(model, models, env, occs, bgPoints) {
#   result <- apply(model, 1, function(x, models, env, occs, bgPoints){
#     choicedModel <- models[[as.integer(x["index"])]]
#     prediction <- dismo::predict(choicedModel, env)
#     auc <- aucCalculator(prediction, occs, bgPoints)
#     return(c(x["settings"], auc))
#   },
#   models = models,
#   env = env,
#   occs = occs,
#   bgPoints = bgPoints)
#
#   result <- data.frame(
#     matrix(unlist(result), nrow = nrow(model), byrow = TRUE),
#     stringsAsFactors = FALSE
#   )
#
#   names(result) <- c("settings", "AUC")
#
#   result <- result %>% mutate(AUC = as.numeric(AUC))
#
#   return(result)
# }
#
#
# # Testing background
# bg.df.test <- bg.df %>%
#   dplyr::filter(isTrain == 0) %>%
#   dplyr::select(x, y)
#
# resultsAUC <- aucStatistcs(modelsAIC0, sp.models@models, env, occsValidacion, bg.df.test)
# write.csv(resultsAUC,
#           file = file.path(outputFolder, "data_auc.csv"),
#           row.names = FALSE)
#
# #### Projections ####
# # predict choicemodel over current climate variables
# predictAndSave <- function(model, models, data, prefix, occs) {
#   choicedModel <- models[[as.integer(model["index"])]]
#   predictions <- dismo::predict(choicedModel, data)
#   raster::writeRaster(predictions,
#                       file.path(outputFolder, paste0(prefix,
#                                                      "log_",
#                                                      model["settings"],
#                                                      "_",
#                                                      outputFolder,
#                                                      "_",
#                                                      ".tif")),
#                       overwrite = TRUE)
#
#   # Threshold prection using minimum traning (min) and 10 percentil (q10) values
#   occsValues <- raster::extract(predictions, occs)
#
#   minValOcc <- min(occsValues, na.rm = TRUE)
#   raster::writeRaster(reclassify(predictions,
#                                  c(-Inf, minValOcc, 0, minValOcc, Inf, 1)),
#                       file.path(outputFolder, paste0(prefix,
#                                                      "bin_min_",
#                                                      model["settings"],
#                                                      "_",
#                                                      outputFolder,
#                                                      "_",
#                                                      ".tif")),
#                       overwrite = TRUE)
#
#   q10ValOcc <- quantile(occsValues, 0.1, na.rm = TRUE)
#   raster::writeRaster(reclassify(predictions,
#                                  c(-Inf, q10ValOcc, 0, q10ValOcc, Inf, 1)),
#                       file.path(outputFolder, paste0(prefix,
#                                                      "bin_q10_",
#                                                      model["settings"],
#                                                      "_",
#                                                      outputFolder,
#                                                      "_",
#                                                      ".tif")),
#                       overwrite = TRUE)
# }
#
#
# # log
# apply(modelsAIC0, 1, predictAndSave,
#       models = sp.models@models, data = env, prefix = "ENM_",
#       occs = occsCalibracion)
