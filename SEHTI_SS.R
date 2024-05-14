sehti_ss = function(PTS_Study, Study_Dates, 
                    counties, waterbodies, roads, structures, roads_buffered,
                    hs_nlcd, nlcd_study_val, trails,
                    CaseType, ControltoCase_Ratio, year, seasons,
                    loop_iteration){
  # ---- 1: Setup Output
  
  ss_folder_name = paste0("ss_", CaseType, "_", ControltoCase_Ratio, "_", paste(seasons, collapse = ""), "_", paste(year, collapse = ""), "_", loop_iteration)
  # Create storage directory for script outputs
  dir.create(ss_folder_name, showWarnings = FALSE)
  ss_folder = here(ss_folder_name)
  
 # ---- 2: Pre-Process Data cont. ----
  
  # 2.4: Filter Data by Case vs. Control
  # Filter Case & Control
  PTS_Case    = PTS_Study %>% filter(status == "Case")
  PTS_Control = PTS_Study %>% filter(status == "Control")
  
  # Filter by Restricted Case (RC) & General Case (GC)
  # RC = attached tick was reason for participation: screen_q4_case == 1
  # GC = attached tick since start of study period:  screen_q4_case == 0 | screen_q4_case_2 == 1 | feedback_tickexp == 1 | tickexp_yn == 1
  if (CaseType == "RC") {
    PTS_Case = PTS_Case %>% filter(CaseType == "RC")
  } else if (CaseType == "GC") {
    PTS_Case = PTS_Case %>% filter(CaseType == "GC")
  } else {
    print("All case points will be used, no filtering has been done to distinguish the Restricted Case (RC) & General Case (GC).")
  }
  
  # Optionally Create a Combination of Cases & Controls
  # Calculate the number of samples needed from PTS_Control
  Control_Sample_Size <- nrow(PTS_Case) * ControltoCase_Ratio
  
  # Randomly sample from PTS_Control
  PTS_Control_sampled <- PTS_Control %>% 
    sample_n(size = Control_Sample_Size)
  
  # Combine PTS_Case and the sampled PTS_Control
  PTS_Model <- bind_rows(PTS_Case, PTS_Control_sampled)
  
  # Inspect Distribution of Participants Cases vs. Controls
  status_summary <- Study_Dates %>%
    group_by(record_id) %>%
    summarise(
      status = first(status),
      screen_qr_case_2 = first(screen_q4_case_2),
      feedback_tickexp = first(feedback_tickexp),
      tickexp_yn = first(tickexp_yn),
      screen_q4_case = first(screen_q4_case),
      .groups = 'drop'
    )
  n_distinct(PTS_Case$ID)
  n_distinct(PTS_Control$ID)
  n_distinct(PTS_Model$ID)
  
  # 2.5: Filter Study Area & Remove Houses, Roads, & Waterbodies
  # Clip by county study area to reduce computational time
  counties    = counties %>% filter(COUNTY == c('HAMPSHIRE', 'HAMPDEN', 'FRANKLIN'))
  waterbodies = st_intersection(waterbodies, counties)
  
  # Ensure Data Layer Projections Match PTS
  counties    = st_transform(counties,    crs = "EPSG:4326")
  waterbodies = st_transform(waterbodies, crs = "EPSG:4326")
  st_crs(structures) == st_crs(counties)
  
  # Process Buffered Roads Uniquely
  if (exists("roads_buffered")){
    print("roads_buffered already exists in the environment.")
  } else if (file.exists(here('data', 'roads_buffered', 'roads_buffered.shp'))){
    roads_buffered = st_read(here('data', 'roads_buffered', 'roads_buffered.shp'))
  } else {
    roads = st_intersection(roads, counties)
    roads = st_transform(roads, crs = "EPSG:4326")
    # Buffer Road Types 1-4
    #   RT1: Limited Access Highway
    #   RT2: Multi-lane Highway, not limited access
    #   RT3: Other numbered route
    #   RT4: Major road - arterials and collectors
    roads_1234 = roads %>% filter(RDTYPE < 5)
    buffer_by_rdtype <- function(df) {
      df %>%
        mutate(geometry = case_when(
          RDTYPE == 1 ~ st_buffer(geometry, dist = 100),
          RDTYPE == 2 ~ st_buffer(geometry, dist = 100),
          RDTYPE == 3 ~ st_buffer(geometry, dist = 20),
          RDTYPE == 4 ~ st_buffer(geometry, dist = 10),
          TRUE ~ geometry  # Keep original geometry if none of the conditions are met
        ))
    }
    roads_buffered <- buffer_by_rdtype(roads_1234)
    st_write(roads_buffered, here('data', 'roads_buffered', 'roads_buffered.shp'))
  }
  
  # Clip & Erase PTS Layer
  # Combine sf objects into one for a simple erase
  erase_sf <- rbind(
    st_sf(geometry = st_geometry(roads_buffered)),
    st_sf(geometry = st_geometry(waterbodies)),
    st_sf(geometry = st_geometry(structures))
  )
  
  # Ensure Valid Geometries
  invalid_geoms <- !st_is_valid(erase_sf)
  table(invalid_geoms)
  
  # If there are any invalid geometries, attempt to repair them
  if(any(invalid_geoms)) {
    erase_sf$geometry[invalid_geoms] <- st_make_valid(erase_sf$geometry[invalid_geoms])
  }
  if(any(!st_is_valid(PTS_Model))) {
    PTS_Model$geometry[!st_is_valid(PTS_Model)] <- st_make_valid(PTS_Model$geometry[!st_is_valid(PTS_Model)])
  } #HUH
  
  # Clip to Study Area
  PTS_Study = st_intersection(PTS_Model, counties)
  
  # Efficient Erase
  PTS_Erase = st_intersects(PTS_Study, erase_sf)
  non_intersecting_points <- which(lengths(PTS_Erase) == 0)
  PTS_Final <- PTS_Study[non_intersecting_points, ]
  
  # Write Output SHP
  PTS_folder = paste0(ss_folder, "/PTS")
  dir.create(paste0(ss_folder, "/PTS"), showWarnings = FALSE)
  st_write(PTS_Final, paste0(PTS_folder, "/PTS.shp"))
  
  # 2.6: Create Case & Geo Files
  ss_case = PTS_Final %>%
    st_drop_geometry() %>%
    mutate(cases = 1) %>%
    dplyr::select(cases, date)
  
  coords <- st_coordinates(PTS_Final)
  PTS_Final$long <- coords[, "X"]
  PTS_Final$lat  <- coords[, "Y"]
  ss_geo = PTS_Final %>%
    st_drop_geometry() %>%
    dplyr::select(lat, long)
  
  # Write Case & Geo File to CSV
  write.csv(ss_case, paste0(ss_folder, "/ss_case.csv"))
  write.csv(ss_geo, paste0(ss_folder, "/ss_geo.csv"))
  
  # ---- 3: Seasonal SatScan ----
  
  # 3.1: Prepare Input
  # Read casfile and geofiles
  casfile = read.csv(paste0(ss_folder, "/ss_case.csv"))
  geofile = read.csv(paste0(ss_folder, "/ss_geo.csv"))
  
  # Subset to year of choice
  if (year[1] == "All") {
    print("No yearwise subset for this run, all years are being used!")
  } else {
    casfile <- casfile %>%
      filter(year(ymd_hms(date)) == year)
    geofile <- geofile %>%
      filter(X %in% casfile$X) 
  }
  
  # 3.2: Create SatScan Function
  ss_run = function(casfile, geofile, season) {
    # 1: Get Date Ranges
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
    
    # 2: Write Input Files
    write.cas(casfile, ss_folder, "CasFile")  # Case File
    write.geo(geofile, ss_folder, "GeoFile")  # Geo File
    
    # 3: Assign Input Files
    casfile_location = paste0(ss_folder, "/CasFile.cas")
    geofile_location = paste0(ss_folder, "/GeoFile.geo")
    
    # 4: Set Parameters
    # Reset Parameter File
    invisible(ss.options(reset=TRUE))
    
    # 5: Set Input Parameters
    ss.options(list(CaseFile =        casfile_location, 
                    CoordinatesFile = geofile_location,
                    CoordinatesType = 1      # Latitude/Longitude 
    ))
    
    # 6: Set Analysis Parameters
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
      MaxTemporalSize = 14,  # Maximum size of 3 weeks
      MaxSpatialSizeInPopulationAtRisk = 10, # x0% alt (50% default)
      CriteriaForReportingSecondaryClusters=0,
      MonteCarloReps = 999 # 99 rep minimum
    ))
    ss.options(c(mindate, maxdate))
    ss.options(c("UseDistanceFromCenterOption=y",
                 "MaxSpatialSizeInDistanceFromCenter=0.125"))
    
    # 7: Set Output Parameters
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
    
    # 8: Write Parameter File
    write.ss.prm(ss_folder, "Parameters")
    
    # 9: Run SatScan
    satscan = satscan(ss_folder, 
                      "Parameters", 
                      sslocation = sslocation, 
                      ssbatchfilename = "SaTScanBatch64",
                      verbose = TRUE)
    
    # 10: Save the result with iteration in the file name
    save(satscan, file = paste0(ss_folder, "/satscan_", season,"_999mcreps", ".rda"))
  }
  
  # 3.3: Run Models for Each Season
  # 3.3.1: Function to determine the season based on the date
  get_season <- function(date, season) {
    month <- as.integer(format(as.Date(date), "%m")) # Extract month and convert to integer
    # Define season ranges
    if (season == "Fall") {
      # Fall includes months 9, 10, and 11
      season_subset <- month %in% c(9, 10, 11)
    } else if (season == "Winter") {
      # Winter includes months 12, 1, and 2
      season_subset <- month %in% c(12, 1, 2)
    } else if (season == "Spring") {
      # Spring includes months 3, 4, and 5
      season_subset <- month %in% c(3, 4, 5)
    } else if (season == "Summer") {
      # Summer includes months 6, 7, and 8
      season_subset <- month %in% c(6, 7, 8)
    } else if (season == "Fall_Winter") {
      # Fall_Winter includes months 9, 10, 11, 12, 1, and 2
      season_subset <- month %in% c(9, 10, 11, 12, 1, 2)
    } else if (season == "Fall_Winter") {
      # Fall_Winter includes months 9, 10, 11, 12, 1, and 2
      season_subset <- month %in% c(9, 10, 11, 12, 1, 2)
    } else if (season == "Winter_Spring") {
      # Winter_Spring includes months 12, 1, 2, 3, 4, and 5
      season_subset <- month %in% c(12, 1, 2, 3, 4, 5)
    } else if (season == "Spring_Summer") {
      # Spring_Summer includes months 3, 4, 5, 6, 7, and 8
      season_subset <- month %in% c(3, 4, 5, 6, 7, 8)
    } else if (season == "Summer_Fall") {
      # Summer_Fall includes months 6, 7, 8, 9, 10, and 11
      season_subset <- month %in% c(6, 7, 8, 9, 10, 11)
    } else if (season == "All") {
      # Include all months for full run
      season_subset <- month %in% c(1:12)
    } else {
      season_subset <- FALSE
    }
    
    return(season_subset)
  }
  get_season_vectorized <- Vectorize(get_season, vectorize.args = "date")
  
  # 3.3.2: Run Seasonal SatScan Analysis Loop
  casfile$Date <- as.Date(casfile$date, format = "%Y-%m-%d %H:%M:%S")
  for (season in seasons) {
    print(paste0("Processing ", season, " Points"))
    # Filter data for the season
    season_casfile <- subset(casfile, get_season_vectorized(Date, season), select = -Date)
    
    # Use the row indices of season_casfile to subset geofile
    season_geofile <- geofile %>%
      filter(X %in% season_casfile$X)
    
    # Determine your max sample size by inspecting number of points per season
    season_counts <- casfile %>%
      mutate(Season = get_season_vectorized(Date, season)) %>%
      group_by(Season) %>%
      summarise(Count = n())
    print(season_counts)
    
    # Run SatScan for each season
    ss_run(season_casfile, season_geofile, season)
  }
  
  # 3.4: Save Shapefiles
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
      shp_path <- file.path(shp_folder, paste0(file_name, ".shp"))
      
      # Write to shapefile
      tryCatch({
        st_write(satscan$shapeclust, shp_path, overwrite = TRUE)
        cat("Successfully written to", shp_path, "\n")
      }, error = function(e) {
        cat("Error in writing to", shp_path, ":", e$message, "\n")
      })
    } else {
      cat("No 'satscan$shapeclust' found in:", file_name, "\n")
    }
  }
  
  # Set the path to the folder where you want to save the shapefiles
  shp_folder <- here(paste0(ss_folder, '/shp'))
  
  # Create the output folder if it doesn't exist
  if (!dir.exists(shp_folder)) {
    dir.create(shp_folder, recursive = TRUE)
  }
  
  # Process each season within the season folder
  process_folder(ss_folder)
  
  # ---- 4: Filter Developed Areas ----
  
  # 4.1: Setup HotSpot Analysis Dependencies
  # Initialize an empty dataframe to store results
  local_proportions <- data.frame(CLUSTER = integer(),
                                  Class = character(),
                                  Count = integer(),
                                  Proportion = numeric())
  
  # NLCD Legend Mapping
  nlcd_legend <- data.frame(
    Red = c(11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95),
    Class = c("Open Water", "Perennial Ice/Snow", 
              "Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity", 
              "Barren Land",
              "Deciduous Forest", "Evergreen Forest", "Mixed Forest", 
              "Dwarf Scrub", "Shrub/Scrub",
              "Herbaceous", "Sedge/Herbaceous", "Lichens", "Moss", 
              "Hay/Pasture", "Cultivated Crops",
              "Woody Wetlands", "Emergent Herbaceous Wetlands")
  )
  
  # 4.2: Calculate HotSpot NLCD Proportions
  # Calculate Proportions Custom Function
  calculate_local_proportions <- function(hs, hs_nlcd, nlcd_legend) {
    local_proportions <- data.frame(CLUSTER = integer(), 
                                    Class = character(), 
                                    Count = integer(), 
                                    Proportion = numeric())
    for (i in 1:nrow(hs)) {
      single_poly <- hs[i, ]
      clipped_raster <- crop(hs_nlcd, single_poly) %>% mask(single_poly)
      nlcd_val_single <- values(clipped_raster)
      prop_single <- na.omit(data.frame(Red = nlcd_val_single)) %>%
        group_by(Red) %>%
        summarise(Count = n()) %>%
        mutate(Proportion = (Count / sum(Count))) %>%
        left_join(nlcd_legend, by = "Red")
      prop_single$CLUSTER <- i
      local_proportions <- rbind(local_proportions, prop_single)
    }
    return(local_proportions)
  }
  
  # 4.3: Remove Developed Land Class HotSpots
  rm_developed <- function(df) {
    # Filter the dataframe to include only the specified classes
    relevant_classes <- c('Developed, Open Space', 'Developed, Low Intensity', 
                          'Developed, Medium Intensity', 'Developed, High Intensity')
    filtered_df <- df %>%
      filter(Class %in% relevant_classes)
    
    # Calculate the sum of Proportion for each CLUSTER for the relevant classes
    proportion_sums <- filtered_df %>%
      group_by(CLUSTER) %>%
      summarize(TotalProportion = sum(Proportion)) %>%
      ungroup()
    
    # Filter CLUSTERs where the total proportion is 0.5 or more
    valid_polygon_ids <- proportion_sums %>%
      filter(TotalProportion < 0.5) %>%
      pull(CLUSTER)
    
    # Filter the original dataframe to include only the valid CLUSTERs
    final_df <- df %>%
      filter(CLUSTER %in% valid_polygon_ids)
    return(final_df)
  }
  
  # 4.4: Combine 4.2 & 4.3 for folder-wide shapefile filtering!
  # Get a list of all shapefile names in the directory
  shapefiles <- list.files(shp_folder, pattern = "\\.shp$", full.names = TRUE)
  
  # Initialize an empty sf object to store combined data
  hs_nlcd_filt <- NULL
  hs_nlcd_repo <- NULL
  
  # Loop through each shapefile in the shp folder
  for(shapefile in shapefiles) {
    # Read the shapefile
    hs <- st_read(shapefile)
    if (nrow(hs) == 0) {
      next
    }
    
    # Calculate local proportions and filter developed land class hotspots
    hs_local_proportions      <- calculate_local_proportions(hs, hs_nlcd, nlcd_legend)
    hs_local_proportions_filt <- rm_developed(hs_local_proportions)
    hs_filt <- unique(hs_local_proportions_filt$CLUSTER)
    hs_filt <- hs %>% filter(CLUSTER %in% hs_filt)
    
    if (nrow(hs_filt) == 0) {
      next
    }
    
    # Store full results for summary statistics
    if(is.null(hs_nlcd_repo)) {
      hs_nlcd_repo <- hs_local_proportions_filt
    } else {
      hs_nlcd_repo <- rbind(hs_nlcd_repo, hs_local_proportions_filt)
    }
    
    # Extract SEASON from the shapefile name
    season <- unlist(regmatches(shapefile, gregexpr("(?<=satscan_)[^_]+", shapefile, perl=TRUE)))
    
    # Add the SEASON column to the filtered shapefile
    hs_filt$SEASON <- season
    
    # Combine with the previously processed shapefiles
    if(is.null(hs_nlcd_filt)) {
      hs_nlcd_filt <- hs_filt
    } else {
      hs_nlcd_filt <- rbind(hs_nlcd_filt, hs_filt)
    }
  }
  
  # 4.5: Report Summary Statistics
  if (is.null(hs_nlcd_filt)) {
    print("There are no hotspots for this data!")
    return()
  } else {
    # 4.5.1: Calculate Summary Statistics for each Class
    class_summary <- hs_nlcd_repo %>%
      group_by(Class) %>%
      summarise(
        MeanProportion   = mean(Proportion, na.rm = TRUE),
        MedianProportion = median(Proportion, na.rm = TRUE),
        MinProportion    = min(Proportion, na.rm = TRUE),
        MaxProportion    = max(Proportion, na.rm = TRUE),
        SDProportion     = sd(Proportion, na.rm = TRUE),
      )
    
    # 4.5.2: Compare to Study Area Values
    nlcd_study_prop = data.frame(Class = nlcd_study_val$NLCD_Land,
                                 Count = nlcd_study_val$Count) %>%
      mutate(Proportion = (nlcd_study_val$Count / sum(nlcd_study_val$Count)))
    class_summary = left_join(class_summary, nlcd_study_prop, by = "Class")
    cols_to_multiply <- c("MeanProportion", "MedianProportion", "MinProportion", "MaxProportion", "SDProportion", "Proportion")
    class_summary[cols_to_multiply] <- class_summary[cols_to_multiply] * 100
    class_summary$Count <- NULL
    names(class_summary)[names(class_summary) == "Proportion"] <- "Study_Area_Proportion"
    write.csv(class_summary, paste0(ss_folder, "/HotSpot_NLCD_Class_Summary.csv"))
  }
  
  # ---- 5: Determine Nearest Trailheads ----
  
  # Determine nearest trailhead to each hotspot
  nearest_trailhead <- function(hs, trails) {
    # Ensure hs and trails are sf objects
    stopifnot("sf" %in% class(hs), "sf" %in% class(trails))
    
    # Find the nearest trail for each polygon in hs
    nearest_indices <- st_nearest_feature(hs, trails)
    
    # Calculate the minimum distance to the nearest trail
    distances <- st_distance(hs, trails[nearest_indices, ], by_element = TRUE)
    
    # Extract the TRAIL_NAME of the nearest trail using the nearest_indices
    nearest_trail_names <- trails$TRAIL_NAME[nearest_indices]
    
    # Add the nearest trail name and distance to hs
    hs$nearest_trail_name <- nearest_trail_names
    hs$nearest_distance <- distances
    
    return(hs)
  }
  trails = st_transform(trails, st_crs(hs_nlcd_filt))
  hs_trails = nearest_trailhead(hs_nlcd_filt, trails)
  hs_trails_500m = hs_trails %>% filter(nearest_distance < units::set_units(500, m))
  st_write(hs_trails, paste0(shp_folder, "/hs_trails_all.shp"))
  st_write(hs_trails_500m, paste0(shp_folder, "/hs_trails_500m.shp"))
  write.csv(st_drop_geometry(hs_trails), paste0(ss_folder, "/hs_trails.csv"))
  write.csv(st_drop_geometry(hs_trails_500m), paste0(ss_folder, "/hs_trails_500m.csv"))
  
}
