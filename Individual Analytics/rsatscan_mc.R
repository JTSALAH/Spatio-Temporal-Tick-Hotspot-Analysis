# ---- 0: Load Packages ----

  require(rsatscan)
  require(sf)
  require(tidyverse)

# ---- 1: Prepare Data ----
casfile = read.csv(here::here("ss_case.csv"))
geofile = read.csv(here::here("ss_geo.csv"))

# Create storage directory
dir.create("ss", showWarnings = FALSE)
ss_folder = here::here("ss")
sslocation = "C:/Program Files/SaTScan"

# ---- 2: Creation Sample SatScan Function ----
sample_ss = function(casfile, geofile, sample_size=100, iteration=1) {
  # ---- 1: Prepare Sample ----
  # Create a directory for model storage
  output_dir <- paste0(ss_folder, "/MonteCarloOutputs")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Sample Input Files
  sample_indices <- sample(1:nrow(casfile), sample_size)
  casfile_sampled = casfile[sample_indices,]
  geofile_sampled = geofile[sample_indices,]
  
  # Get Date Ranges
  maxdate <- casfile_sampled$date %>%
    as.Date(format = "%m/%d/%Y %H:%M:%S") %>%
    max() %>%
    format("%Y/%m/%d") %>%
    gsub("/0", "/", .) 
  maxdate = paste("EndDate=", maxdate, sep="")
  print(maxdate)
  
  mindate <- casfile_sampled$date %>%
    as.Date(format = "%m/%d/%Y %H:%M:%S") %>%
    min() %>%
    format("%Y/%m/%d") %>%
    gsub("/0", "/", .)
  mindate = paste("StartDate=", mindate, sep="")
  print(mindate)
  
  # Write Input Files
  write.cas(casfile_sampled, ss_folder, "CasFile")  # Case File
  write.geo(geofile_sampled, ss_folder, "GeoFile")  # Geo File
  
  # Assign Input Files
  casfile_location = paste0(ss_folder, "/CasFile.cas")
  geofile_location = paste0(ss_folder, "/GeoFile.geo")
  
  # ---- 2: Set Parameters ----
  # Reset Parameter File
  invisible(ss.options(reset=TRUE))
  
  # Set Input Parameters
  ss.options(list(CaseFile =        casfile_location, 
                  CoordinatesFile = geofile_location,
                  CoordinatesType = 1      # Latitude/Longitude 
  ))
  
  # Set Analysis Parameters
  ss.options(list(
    AnalysisType = 3, # 3 = Retrospective Space-Time, 7 = Seasonal Temporal
    ModelType = 2,    # 1 = Bernoulli, 2 = Space-Time Permutation
    ScanAreas = 3     # Both High & Low Rates Detect
  ))
  ss.options(list(# Try none, generic, and month
    PrecisionCaseTimes = 3,    # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
    TimeAggregationUnits = 3,  # 3 = Day
    TimeAggregationLength = 7, # 7 days = 1wk
    MinimumTemporalClusterSize = 14,  # Minimum size of 14 days
    MaxTemporalSize = 30,  # 1mo
    MaxSpatialSizeInPopulationAtRisk = 10, # 10% alt (50% default)
    CriteriaForReportingSecondaryClusters=5,
    MonteCarloReps = 0 # 99 rep minimum
  ))
  ss.options(c(mindate, maxdate))
  
  # Set Output Parameters
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
  
  # Write Parameter File
  write.ss.prm(ss_folder, "Parameters")
  
  # ---- 3: Run SatScan ----
  satscan = satscan(ss_folder, 
                    "Parameters", 
                    sslocation = sslocation, 
                    ssbatchfilename = "SaTScanBatch64",
                    verbose = TRUE)
  print(paste0("satscan_sample", sample_size, "_rep", iteration, " has been written!"))
  
  # Save the result with iteration in the filename
  save(satscan, file = paste0(output_dir, "/satscan_sample", sample_size, "_rep", iteration, ".rda"))
}

# Monte Carlo SatScan
for (i in 1:99) {
  sample_ss(casfile, geofile, sample_size=25000, iteration=i)
}



# ---- 3: Save Shapefiles ----

# Set the path to your folder containing the .rda files
rda_folder <- here::here('ss', 'MonteCarloOutputs')

# Set the path to the folder where you want to save the shapefiles
dir.create(here::here('ss', 'mc_shp'))
output_folder <- here::here('ss', 'mc_shp')

# Create the output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# List all .rda files in the folder
rda_files <- list.files(rda_folder, pattern = "\\.rda$", full.names = TRUE)

# Function to process each file
process_file <- function(file_path) {
  # Load the .rda file
  load(file_path)
  
  # Check if 'satscan$shapeclust' exists
  if (exists("satscan") && !is.null(satscan$shapeclust)) {
    # Extract the file name without extension
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Create the path for the new shapefile
    shp_path <- file.path(output_folder, paste0(file_name, ".shp"))
    
    # Write to shapefile
    st_write(satscan$shapeclust, shp_path)
    cat("Written shapefile for:", file_name, "\n")
  } else {
    cat("No 'satscan$shapeclust' found in:", file_name, "\n")
  }
}

# Apply the function to each file
results <- lapply(rda_files, process_file)

# ---- 4: Analyze MC Consistency ----

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
