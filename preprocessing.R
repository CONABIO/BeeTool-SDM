#' BeeTool SDM script

# Libraries ----
library(terra)

library(fs)
library(readr)
library(dplyr)
library(purrr)
library(stringr)

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

## Load data ----
if (config$regional$use_regional_cutoff) {
  regional_data <- vect(config$regional$shapefile_region_path)
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

input_data_file <- args[1]

output_folder <- input_data_file %>%
  path_file() %>%
  path_ext_remove()

output_folder <- path_join(c(base_output, output_folder))

if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

### Load occurrence data ----
crs_wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
occs_data <- read_csv(input_data_file)

# Coordinates columns are X and Y
occs_data <- occs_data %>%
  vect(geom=c("x", "y"), crs=crs_wgs84)
sampled_occs_data <- occs_data
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
  regions_of_interest <- extract(regional_data, occs_data)
  writeVector(
    regions_of_interest,
    path_join(c(output_folder, "region_of_interest.shp")),
    overwrite=TRUE)
}