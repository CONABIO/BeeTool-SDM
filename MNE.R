#' BeeTool SDM script

# Libraries ----
library(magrittr)
library(fs)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(terra)

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
}

## Load covariables ----
if (config$covariables$is_worldclim) {
  covar_file_list <- str_c(
    config$covariables$path,
    "/wc2.1_30s_",
    config$covariables$variables,
    ".tif")
}

covar_rasters <- rast(covar_file_list)
if (!is.null(regions_of_interest)) {
  covar_rasters <- crop(covar_rasters, ext(regions_of_interest))
  gc()
}

sampled_occs_data_covar <- terra::extract(
  covar_rasters,
  sampled_occs_data,
  bind=TRUE)


#### VARIABLES + PRESENCIAS####
covarData <- raster::extract(enviromentalVariables, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>%
  dplyr::select_(.dots=names(enviromentalVariables)) %>%
  complete.cases
covarData <- covarData[completeDataCases, ]

####SELECCION DE VARIABLES####
speciesCol <- match("Presence", names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)

select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.txt"))

selectedVariables <- enviromentalVariables[[select_var]]
selectedVariablesAOI <- enviromentalVariablesAOI[[select_var]]

# Selects the M of the species, base on Olson´s ecoregions
# Download: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
# Intersects the occurrence data with polygons
ecoregionsOfInterest <- sp::over(occsData, regionalizacion) %>%
  filter(!is.na(ECO_ID))

idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
polygonsOfInterest <- regionalizacion[regionalizacion$ECO_ID %in% idsEcoRegions, ]
writeOGR(polygonsOfInterest, layer = 'ecoregionsOI', outputFolder, driver="ESRI Shapefile")

# Mask present rasters with ecoregions of interest
selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
env <- raster::mask(selectedVariablesCrop,
                    polygonsOfInterest) #Species variables delimited by M

dir.create(file.path(outputFolder, "Presente"))
writeRaster(env,
            file.path(outputFolder, "Presente/.tif"),
            bylayer = T, suffix='names',
            overwrite = TRUE)

#### Calibration ####
# Divides your data into trainining and test data sets. 70/30 %
sampleDataPoints <- sample.int(
  nrow(covarData),
  size = floor(0.7*nrow(covarData))
)

selectedValues <- rep(0, nrow(covarData)) %>% inset(sampleDataPoints, 1)

covarData$isTrain <- selectedValues
write.csv(cbind(covarData@data, coordinates(covarData)), file = file.path(outputFolder,
                                                                          paste0(outputFolder,
                                                                                 "_",
                                                                                 "presencias",
                                                                                 ".csv")),
          row.names = FALSE)
# MAXENT calibration
# We used ENMeval package to estimate optimal model complexity (Muscarrella et al. 2014)
# Modeling process, first separate the calibration and validation data
occsCalibracion <- covarData %>%
  as.data.frame() %>%
  dplyr::filter(isTrain == 1) %>%
  dplyr::select(Long, Lat)

write.csv(occsCalibracion, file = file.path(outputFolder,paste0(outputFolder,
                                                                "_","Calibracion",".csv")),
          row.names = FALSE)

occsValidacion <- covarData %>%
  as.data.frame() %>%
  dplyr::filter(isTrain == 0) %>%
  dplyr::select(Long, Lat)

write.csv(occsValidacion, file = file.path(outputFolder,paste0(outputFolder,
                                                               "_","Validacion",".csv")),
          row.names = FALSE)

# Background
bg.df <- dismo::randomPoints(env[[1]], n = 10000) %>% as.data.frame()

#Divide backgeound into train and test
sample.bg <- sample.int(
  nrow(bg.df),
  size = floor(0.7*nrow(bg.df))
)
selectedValues.bg <- rep(0, nrow(bg.df)) %>% inset(sample.bg, 1)

bg.df$isTrain <- selectedValues.bg

sp::coordinates(bg.df) <- c("x", "y")
sp::proj4string(bg.df) <- crs.wgs84
bg.dfbio <- raster::extract(enviromentalVariables, bg.df)

bg.df<-as.data.frame(bg.df)
bg.dfbio <- cbind(bg.df, bg.dfbio) %>% as.data.frame()

write.csv(bg.dfbio, file = file.path(outputFolder, paste0(outputFolder,
                                                          "_",
                                                          "background_points",
                                                          ".csv")))

#training background
bg.df.cal <- bg.df %>%
  dplyr::filter(isTrain == 1) %>%
  dplyr::select(x, y)
write.csv(bg.df.cal, file = file.path(outputFolder,paste0(outputFolder,
                                                          "_","back_calibracion",".csv")),
          row.names = FALSE)

#testing back
bg.df.val <- bg.df %>%
  dplyr::filter(isTrain == 0) %>%
  dplyr::select(x, y)
write.csv(bg.df.val, file = file.path(outputFolder,paste0(outputFolder,
                                                          "_","back_validacion",".csv")),
          row.names = FALSE)


# ENMeval
sp.models <- ENMevaluate(occsCalibracion, env, bg.df.cal, RMvalues = seq(0.5, 4, 0.5),
                         fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                         method = "randomkfold", kfolds = 4, bin.output = TRUE,
                         parallel = TRUE, numCores = parallel::detectCores()-1,
                         updateProgress = TRUE)

resultados_enmeval <- sp.models@results

saveRDS(sp.models@models,
        file=file.path(outputFolder,"Maxent_models.Rds"))

write.csv(resultados_enmeval,
          file = file.path(outputFolder, "enmeval_results.csv"),
          row.names = FALSE)
sp.models_p<-sp.models@predictions

dir.create(file.path(outputFolder, "Outputs_todos"))
writeRaster(sp.models_p, file = file.path(outputFolder, "Outputs_todos/", paste0(outputFolder)),
            suffix='names',
            format = "GTiff",
            bylayer=TRUE,
            overwrite= TRUE)

# delta_aic <- which(resultados_enmeval$delta.AICc == 0)
modelsAIC0 <- resultados_enmeval %>%
  mutate(index = rownames(resultados_enmeval)) %>%
  filter(delta.AICc == 0) %>%
  select(index, settings) %>%
  mutate(index = as.numeric(index), settings = as.character(settings))

aic.opt <- sp.models@models[[which(sp.models@results$delta.AICc==0)]]
importa <- var.importance(aic.opt)
write.csv( importa,
           file = file.path(outputFolder, "varImportance.csv"),
           row.names = FALSE)

# save species niche (raw output) model over raster
saveRasterWithSettings <- function(models, predictions, prefix) {
  raster::writeRaster(predictions[[models["settings"]]],
                      file.path(outputFolder, paste0(prefix,
                                                     models["settings"],
                                                     ".tif")),
                      overwrite = TRUE)
}

apply(modelsAIC0, 1, saveRasterWithSettings,
      predictions = sp.models@predictions, prefix = "ENM_prediction_M_raw_")


####ENMTest####
#source("funciones_LAE.R")
#Threslhold independent

#AUC
aucCalculator <- function(prediction, occs, bgPoints) {
  data <- rbind(occs, setNames(bgPoints, names(occs)))
  labels <- c(rep(1, nrow(occs)),
              rep(0, nrow(bgPoints)))
  scores <- raster::extract(prediction, data)
  pred <- ROCR::prediction(scores, labels)
  # perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  return(auc)
}

aucStatistcs <- function(model, models, env, occs, bgPoints) {
  result <- apply(model, 1, function(x, models, env, occs, bgPoints){
    choicedModel <- models[[as.integer(x["index"])]]
    prediction <- dismo::predict(choicedModel, env)
    auc <- aucCalculator(prediction, occs, bgPoints)
    return(c(x["settings"], auc))
  },
  models = models,
  env = env,
  occs = occs,
  bgPoints = bgPoints)

  result <- data.frame(
    matrix(unlist(result), nrow = nrow(model), byrow = TRUE),
    stringsAsFactors = FALSE
  )

  names(result) <- c("settings", "AUC")

  result <- result %>% mutate(AUC = as.numeric(AUC))

  return(result)
}


# Testing background
bg.df.test <- bg.df %>%
  dplyr::filter(isTrain == 0) %>%
  dplyr::select(x, y)

resultsAUC <- aucStatistcs(modelsAIC0, sp.models@models, env, occsValidacion, bg.df.test)
write.csv(resultsAUC,
          file = file.path(outputFolder, "data_auc.csv"),
          row.names = FALSE)

#### Projections ####
# predict choicemodel over current climate variables
predictAndSave <- function(model, models, data, prefix, occs) {
  choicedModel <- models[[as.integer(model["index"])]]
  predictions <- dismo::predict(choicedModel, data)
  raster::writeRaster(predictions,
                      file.path(outputFolder, paste0(prefix,
                                                     "log_",
                                                     model["settings"],
                                                     "_",
                                                     outputFolder,
                                                     "_",
                                                     ".tif")),
                      overwrite = TRUE)

  #Threshold prection using minimum traning (min) and 10 percentil (q10) values
  occsValues <- raster::extract(predictions, occs)

  minValOcc <- min(occsValues, na.rm = TRUE)
  raster::writeRaster(reclassify(predictions,
                                 c(-Inf, minValOcc, 0, minValOcc, Inf, 1)),
                      file.path(outputFolder, paste0(prefix,
                                                     "bin_min_",
                                                     model["settings"],
                                                     "_",
                                                     outputFolder,
                                                     "_",
                                                     ".tif")),
                      overwrite = TRUE)

  q10ValOcc <- quantile(occsValues, 0.1, na.rm = TRUE)
  raster::writeRaster(reclassify(predictions,
                                 c(-Inf, q10ValOcc, 0, q10ValOcc, Inf, 1)),
                      file.path(outputFolder, paste0(prefix,
                                                     "bin_q10_",
                                                     model["settings"],
                                                     "_",
                                                     outputFolder,
                                                     "_",
                                                     ".tif")),
                      overwrite = TRUE)
}


# log
apply(modelsAIC0, 1, predictAndSave,
      models = sp.models@models, data = env, prefix = "ENM_",
      occs = occsCalibracion)
