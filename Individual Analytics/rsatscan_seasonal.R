# ---- 0: Load Packages ----

require(rsatscan)
require(sf)
require(tidyverse)

# ---- 1: Prepare Data ----

# Read casfile and geofiles
casfile = read.csv(here::here("ss_case.csv"))
geofile = read.csv(here::here("ss_geo.csv"))

# Subset to year of choice
# casfile <- casfile %>%
  # filter(year(ymd_hms(date)) == 2022) # 2023 Data Only
# geofile <- geofile %>%
  # filter(X %in% casfile$X)

# Create storage directory
dir.create("ss", showWarnings = FALSE)
ss_folder = here::here("ss")
sslocation = "C:/Program Files/SaTScan"

# ---- 2: Create Monte Carlo Sample SatScan Function ----

ss_run = function(casfile, geofile, season) {
  # 1.2: Get Date Ranges
  maxdate <- casfile$date %>%
    as.Date(format = "%Y-%m-%d %H:%M:%S") %>%
    max() %>%
    format("%Y/%m/%d") %>%
    gsub("/0", "/", .) 
  maxdate = paste("EndDate=", maxdate, sep="")
  print(maxdate)
  
  mindate <- casfile$date %>%
    as.Date(format = "%Y-%m-%d %H:%M:%S") %>%
    min() %>%
    format("%Y/%m/%d") %>%
    gsub("/0", "/", .)
  mindate = paste("StartDate=", mindate, sep="")
  print(mindate)
  
  # 1.3: Write Input Files
  write.cas(casfile, ss_folder, "CasFile")  # Case File
  write.geo(geofile, ss_folder, "GeoFile")  # Geo File
  
  # 1.4: Assign Input Files
  casfile_location = paste0(ss_folder, "/CasFile.cas")
  geofile_location = paste0(ss_folder, "/GeoFile.geo")
  
  # 2: Set Parameters
  # Reset Parameter File
  invisible(ss.options(reset=TRUE))
  
  # 2.1: Set Input Parameters
  ss.options(list(CaseFile =        casfile_location, 
                  CoordinatesFile = geofile_location,
                  CoordinatesType = 1      # Latitude/Longitude 
  ))
  
  # 2.2: Set Analysis Parameters
  ss.options(list(
    AnalysisType = 3, # 3 = Retrospective Space-Time, 7 = Seasonal Temporal
    ModelType = 2,    # 1 = Bernoulli, 2 = Space-Time Permutation
    ScanAreas = 3     # Both High & Low Rates Detect
  ))
  ss.options(list(# Try none, generic, and month
    PrecisionCaseTimes = 3,    # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
    TimeAggregationUnits = 3,  # 3 = Day
    TimeAggregationLength = 7, # 7 days = 1wk
    MinimumTemporalClusterSize = 7,  # Minimum size of 7 days
    MaxTemporalSize = 14,  # Maximum size of 14 days
    MaxSpatialSizeInPopulationAtRisk = 10, # 10% alt (50% default)
    CriteriaForReportingSecondaryClusters=0,
    MonteCarloReps = 999 # 99 rep minimum
  ))
  ss.options(c(mindate, maxdate))
  ss.options(c("UseDistanceFromCenterOption=y",
               "MaxSpatialSizeInDistanceFromCenter=0.125"))
  
  # 2.3: Set Output Parameters
  ss.options(list(ResultsFile = "SSResult.txt",
                  OutputGoogleEarthKML = "n",
                  OutputShapefiles = "y",
                  MostLikelyClusterEachCentroidASCII = "y", # Output Hotspot Raster
                  MostLikelyClusterEachCentroidDBase = "n",
                  MostLikelyClusterCaseInfoEachCentroidASCII = "y",
                  MostLikelyClusterCaseInfoEachCentroidDBase = "n"
  ))
  
  # Optional: Set Advanced Options Parameters
  ss.options(list(NumberParallelProcesses = 0))
  
  # 2.4: Write Parameter File
  write.ss.prm(ss_folder, "Parameters")
  
  # 3: Run SatScan
  satscan = satscan(ss_folder, 
                    "Parameters", 
                    sslocation = sslocation, 
                    ssbatchfilename = "SaTScanBatch64",
                    verbose = TRUE)
  
  # Save the result with iteration in the filename
  save(satscan, file = paste0(ss_folder, "/satscan_", season,"_999mcreps", ".rda"))
}

# ---- 3: Run Models for Each Season ----

# 3.0: Function to determine the season based on the date
get_season <- function(date) {
  month <- month(date)
  ifelse(month %in% c(3, 4, 5), "Spring",
         ifelse(month %in% c(6, 7, 8), "Summer",
                ifelse(month %in% c(9, 10, 11), "Fall", "Winter")))
}
get_season_vectorized <- Vectorize(get_season)

spring_casfile <- casfile %>%
  filter(get_season_vectorized(Date) == "Spring")

# 3.1: Run Seasonal SatScan Analysis Loop
seasons <- c("Spring", "Summer", "Fall", "Winter")
casfile$Date <- as.Date(casfile$date, format = "%Y-%m-%d %H:%M:%S")
for (season in seasons) {
  print(paste0("Processing ", season, " Points"))
  # Filter data for the season
  season_casfile <- casfile %>%
    filter(get_season_vectorized(Date) == season)
  
  # Use the row indices of season_casfile to subset geofile
  season_geofile <- geofile %>%
    filter(X %in% season_casfile$X)
  
  # Determine your max sample size by inspecting number of points per season
  season_counts <- casfile %>%
    mutate(Season = get_season_vectorized(Date)) %>%
    group_by(Season) %>%
    summarise(Count = n())
  print(season_counts)
  
  # Run SatScan for each season
  ss_run(season_casfile, season_geofile, season)
}


# ---- 3: Save Shapefiles ----

# Function to process files in a given folder
process_folder <- function(folder) {
  # List all .rda files in the folder
  rda_files <- list.files(folder, pattern = "\\.rda$", full.names = TRUE)
  
  # Apply the process_file function to each file with the season name
  lapply(rda_files, process_file)
}

# Function to process each file
process_file <- function(file_path) {
  # Load the .rda file
  load(file_path)
  
  # Check if 'satscan$shapeclust' exists
  if (exists("satscan") && !is.null(satscan$shapeclust)) {
    # Extract the file name without extension
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Create the path for the new shapefile with season in the name
    shp_path <- file.path(output_folder, paste0(file_name, ".shp"))
    
    # Write to shapefile
    st_write(satscan$shapeclust, shp_path, overwrite = TRUE)
    cat("Written shapefile for:", shp_path, "\n")
  } else {
    cat("No 'satscan$shapeclust' found in:", file_name, "\n")
  }
}

# Set the path to the folder where you want to save the shapefiles
output_folder <- here::here('ss', 'shp')

# Create the output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Process each season within the season folder
process_folder(ss_folder)

# ---- 4: Analyze Consistency ----

analyze_shapefile_overlap <- function(folder_path) {
  # List all shapefiles in the folder
  shapefile_paths <- list.files(folder_path, pattern = "\\.shp$", full.names = TRUE)
  
  # Read shapefiles
  shapefiles <- lapply(shapefile_paths, st_read)
  
  # Initialize a dataframe to store results
  overlap_results_df <- data.frame(
    Shapefile1 = character(),
    Shapefile2 = character(),
    OverlapArea = numeric(),
    PercentageOverlap = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Pairwise Intersection and Area Calculation
  combn(length(shapefiles), 2, function(idx) {
    intersect_shapefiles <- st_intersection(shapefiles[[idx[1]]], shapefiles[[idx[2]]])
    overlap_area <- sum(st_area(intersect_shapefiles))
    shapefile1_area <- sum(st_area(shapefiles[[idx[1]]]))
    
    # Calculate percentage overlap
    percentage_overlap <- (overlap_area / shapefile1_area) * 100
    
    # Add to dataframe
    overlap_results_df <<- rbind(overlap_results_df, data.frame(
      Shapefile1 = basename(shapefile_paths[idx[1]]),
      Shapefile2 = basename(shapefile_paths[idx[2]]),
      OverlapArea = as.numeric(overlap_area),
      PercentageOverlap = as.numeric(percentage_overlap)
    ))
  }, simplify = FALSE)
  
  # Return the dataframe
  return(overlap_results_df)
}

overlap_results <- analyze_shapefile_overlap(output_folder)

overlap_summary <- overlap_results %>%
  group_by(Shapefile1) %>%
  summarize(
    MeanPercentage = mean(PercentageOverlap, na.rm = TRUE),
    MedianPercentage = median(PercentageOverlap, na.rm = TRUE),
    SDPercentage = sd(PercentageOverlap, na.rm = TRUE),
    MinPercentage = min(PercentageOverlap, na.rm = TRUE),
    MaxPercentage = max(PercentageOverlap, na.rm = TRUE),
    Q1Percentage = quantile(PercentageOverlap, 0.25, na.rm = TRUE),
    Q3Percentage = quantile(PercentageOverlap, 0.75, na.rm = TRUE)
  )
