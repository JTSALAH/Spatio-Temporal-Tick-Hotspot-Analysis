# Input Data & Parameters

Data
- Folder of all participants .gpx files 
- MassGIS Data: Counties
- MassGIS Data: Major Ponds and Major Streams
- MassGIS Data: Building Structures (2-D)
- MassGIS Data: MassGIS-MassDOT Roads
- MassGIS Data: Tracks and Trails
- MRLC's National Land Cover Database

SatScan Macro Parameters
- Determine Case Type: Attached vs. Detected
- Specify Control to Case Ratio
- Choose analytical year range: For example, just 2022, 2022 & 2023, or 'All' years
- Select seasons of interest: For example, just Fall, Fall and Winter, Spring Summer and Fall etc.

# Format Data

1. Read all .gpx files from a folder and combine them together into a single master sf object
2. Utilize adehabitatLT package to create a trajectory dataframe
3. Removal of any data with lag time greater than 900 seconds
4. Read in Cases & Contacts reference table
5. Join Cases & Contacts table to points
6. Filter observations based on study dates
7. Categorize participants as cases or controls based on survey data
8. Read counties polygon, water bodies polygon, building footprint polygon, and roads line.
9. Clip points to the extent of W. Mass counties (Hampshire, Hampden, Franklin)
10. Erase all points in water bodies
11. Buffer structures by 10m and erase all points within 10m of structures
12. Erase points on road types 1-4, with each road type being buffered appropriately.
- 1 - Limited Access Highway
  - 100 m
- 2 - Multi-lane Highway, not limited access
  - 100 m
- 3 - Other numbered route
  - 20 m
- 4 - Major road - arterials and collectors
  - 10 m
13. Format a SatScan case and geo file based on the final output.

# SatScan Parameters

## Model Type
- AnalysisType = 3, # 3 = Retrospective Space-Time
- ModelType = 2,    # 2 = Space-Time Permutation
- ScanAreas = 3     # 3 = Both High & Low Rates Detect

## Model Settings
- PrecisionCaseTimes = 3,                  # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
- TimeAggregationUnits = 3,                # 3 = Day
- TimeAggregationLength = 7,               # 7 days = 1wk
- MinimumTemporalClusterSize = 7,          # Minimum size of 7 days
- MaxTemporalSize = 14,                    # Maximum size of 14 days (2wks)
- MaxSpatialSizeInPopulationAtRisk = 10,   # (50% default)
- CriteriaForReportingSecondaryClusters = 0,
- MonteCarloReps = 999
- MaxSpatialSizeInDistanceFromCenter=0.125 # Maximum hotspot radius = 125 meters

```r
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
      PrecisionCaseTimes = 3,                # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
      TimeAggregationUnits = 3,              # 3 = Day
      TimeAggregationLength = 7,             # 7 days = 1wk
      MinimumTemporalClusterSize = 7,        # Minimum size of 7 days
      MaxTemporalSize = 14,                  # Maximum size of 3 weeks
      MaxSpatialSizeInPopulationAtRisk = 10, # (50% default)
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
```

# Hotspot Analytics
## Remove Developed Hotspots

It is difficult to organize tick treatments on private property, so we chose to remove hotspots which were present on more than 50% developed land cover.

```r
 # 1: Setup HotSpot Analysis Dependencies
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
  
  # 2: Calculate HotSpot NLCD Proportions
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
  
  # 3: Remove Developed Land Class HotSpots
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
```

## Determine Nearest Trailhead to Hotspots

In looking for appropriate tick sites, we chose to find local trails to the tick hotspots. Our criteria was that the spatio-temporal hotspot must be within 500 meters of the trail!
 
```r
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
```

# Iterative SatScan Hotspot Runs

The analytical framework custom function sehti_ss() is stored in the SEHTI_SS.R file for organizational purposes, which we can import.
This sehti_ss() function contains a random sampling logic which allows us to run our hotspot analysis with different samples of the controls to account for noise.

```r
  # Create a Combination of Cases & Controls
  # Calculate the number of samples needed from PTS_Control
  Control_Sample_Size <- nrow(PTS_Case) * ControltoCase_Ratio
  
  # Randomly sample from PTS_Control
  PTS_Control_sampled <- PTS_Control %>% 
    sample_n(size = Control_Sample_Size)
  
  # Combine PTS_Case and the sampled PTS_Control
  PTS_Model <- bind_rows(PTS_Case, PTS_Control_sampled)
```

Each run of the sehti_ss() function will result in a different hotspot result, allowing us to iteratively run the SatScan models and rule out noise from controls at various ratios.

```r
  source(here('SEHTI_SS.R'))
  
  num_iters = 5
  for (i in 1:num_iters) {
    sehti_ss(PTS_Study, Study_Dates, 
             counties, waterbodies, roads, structures, roads_buffered,
             hs_nlcd, nlcd_study_val, trails,
             CaseType, ControltoCase_Ratio, year, seasons,
             loop_iteration = i)
  }
```

## Iterations Done

Case Types
- Restricted cases (RC): Attached tick was reason for visit/participation.
- General cases (GC): Attached tick since start of study period.

Iterations
- RC 1:1 Case:Control Results (10 Repetitions)
- RC 1:2 Case:Control Results (5 Repetitions)
- GC 1:1 Case:Control Results (5 Repetitions)
- GC 1:2 Case:Control Results (5 Repetitions) 

# Example Output

I configured this script to create a folders/subfolders for each iteration in order to organize and summarize model outputs.

In the example output below, we store:
- A folder with a shapefile of our model points used 
- A folder of shapefiles for our hotspot
- The SatScan Case & Geo files in both .cas/.geo & .csv
- Land Cover Class summary for hotspots
- A table of all hotspots and their closest trails
- An R friendly .rda file for each SatScan output

![ExampleOutput](https://raw.githubusercontent.com/JTSALAH/Spatio-Temporal-Tick-Hotspot-Analysis/main/IMAGES/ExampleOutput.png)







