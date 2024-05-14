# ---- 0: Load Packages ----

  require(sf)
  require(tidyverse)
  require(terra)

# ---- 1: Prepare Data ----
# ---- 1.1: Single Run Option ----
  
  # Load Data
  hs = st_read(here::here('ss', 'mc_shp', 'satscan_sample25000_rep1.shp'))
  hs = st_read(here::here('seasonal_hs_trails','seasonal_hs_trails.shp'))
  nlcd = rast(here::here('data', 'nlcd_2021', 
                         'nlcd_wgs84.tif'))
  
  # Ensure CRS
  crs = crs(hs)
  nlcd = project(nlcd, crs)
  
# ---- 1.2: Batch Run Option ----
  
  # shapefile_folder <- here::here('ss', 'mc_shp')
  # shapefile_folder <- here::here('ss', 'shp')
  # shapefile_folder <- here::here('ss', 'shp', 'Fall')
  shapefile_folder <- here::here('ss', 'shp', 'Spring')
  # shapefile_folder <- here::here('ss', 'shp', 'Summer')
  shapefiles <- list.files(shapefile_folder, pattern = "\\.shp$", full.names = TRUE)
  
  sf_list <- shapefiles %>% 
    map(~st_read(.x, quiet = TRUE))
  
  hs <- bind_rows(sf_list)
  
# ---- 2.1.1: NLCD Statistics (Local) ----
  
  # Clip NLCD with the multipolygon
  # hs_nlcd <- crop(nlcd, hs)
  # hs_nlcd <- mask(hs_nlcd, hs)
  # writeRaster(hs_nlcd, 'hs_nlcd.tif')
  
  # Optional - Read in Clipped NLCD
  hs_nlcd = rast(here::here('data', 'nlcd_2021', 
                            'nlcd_wgs84.tif'))
  
  # Initialize an empty dataframe to store results
  local_proportions <- data.frame(PolygonID = integer(),
                                  Class = character(),
                                  Count = integer(),
                                  Proportion = numeric())
  
  # NLCD legend mapping
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
  
  # Loop through each polygon in the multipolygon
  for (i in 1:nrow(hs)) {
    # Extract the single polygon
    single_poly <- hs[i, ]
    
    # Clip and mask the raster with the single polygon
    clipped_raster <- crop(hs_nlcd, single_poly) %>% mask(single_poly)
    
    # Calculate statistics for the single polygon
    nlcd_val_single <- values(clipped_raster)
    prop_single <- na.omit(data.frame(Red = nlcd_val_single)) %>%
      group_by(Red) %>%
      summarise(Count = n()) %>%
      mutate(Proportion = (Count / sum(Count))) %>% # Opt: Can multiply by 100
      left_join(nlcd_legend, by = "Red") %>%
      select(Class, Count, Proportion)
    
    # Add the polygon ID to the dataframe
    prop_single$PolygonID <- i
    
    # Bind the results to the main dataframe
    local_proportions <- rbind(local_proportions, prop_single)
  }
  
  # Ensure proportions are calculated correctly (add up to 100)
  total_proportions <- local_proportions %>%
    group_by(PolygonID) %>%
    summarise(TotalProportion = sum(Proportion))
  
  # Print the local proportions dataframe
  print(local_proportions)

# ---- 2.1.2: NLCD Proportions Summary Statistics ----
  
  # Summary statistics for each Class
  class_summary <- local_proportions %>%
    group_by(Class) %>%
    summarise(
      MeanProportion   = mean(Proportion, na.rm = TRUE),
      MedianProportion = median(Proportion, na.rm = TRUE),
      MinProportion    = min(Proportion, na.rm = TRUE),
      MaxProportion    = max(Proportion, na.rm = TRUE),
      SDProportion     = sd(Proportion, na.rm = TRUE),
    )
  
  # Compare to Study Area Values
  nlcd_study_val = read.csv(here::here('data', 'nlcd_2021','NLCD_Study_Values.csv'))
  nlcd_study_prop = data.frame(Class = nlcd_study_val$NLCD_Land,
                               Count = nlcd_study_val$Count) %>%
    mutate(Proportion = (nlcd_study_val$Count / sum(nlcd_study_val$Count))) # Opt: Can multiply by 100
  class_summary = left_join(class_summary, nlcd_study_prop, by = "Class")
  
  # Print the summary statistics
  view(class_summary)
  
  # Assuming your dataframe is named class_summary
  
  # Multiply columns by 100
  cols_to_multiply <- c("MeanProportion", "MedianProportion", "MinProportion", "MaxProportion", "SDProportion", "Proportion")
  class_summary[cols_to_multiply] <- class_summary[cols_to_multiply] * 100
  
  # Remove the Count column
  class_summary$Count <- NULL
  
  # Rename the Proportion column to Study_Proportion
  names(class_summary)[names(class_summary) == "Proportion"] <- "Study_Proportion"
  
  # Now your dataframe class_summary is updated as per your requirements
  
  
# ---- 2.1.3: NLCD Proportions Chi-Sqr ----
  
  # Merge the local proportions with the study area proportions
  # comparison_data <- local_proportions %>%
    # left_join(nlcd_study_prop, by = "Class", suffix = c("_hotspot", "_study"))
  
  # ChiSqr is a count data metric, so we need to standardize the proportion percentages
  base_number <- sum(nlcd_study_prop$Count) # Standardize hotspots to study area proportions
  standardize_proportion <- function(proportion, base) {
    count <- proportion * base
    if (count - floor(count) > 0.5) { # Round up/down to maintain count format
      return(ceiling(count))
    } else {
      return(floor(count))
    }
  }
  
  # Apply to local_proportions & nlcd_study_prop
  local_proportions$STDCount <- mapply(standardize_proportion, 
                                       local_proportions$Proportion, 
                                       MoreArgs = list(base = base_number))
  nlcd_study_prop$STDCount <- mapply(standardize_proportion, 
                                     nlcd_study_prop$Proportion, 
                                     MoreArgs = list(base = base_number)) # This is the calculation proof
  
  # Create Results DF to store output
  chisqr_results <- data.frame(Class = character(), 
                               ChiSquared = numeric(), 
                               pvalue = numeric(),
                               stringsAsFactors = FALSE)
  
  for (class in unique(local_proportions$Class)) {
    # Extract counts for the current class
    hotspot_count <- mean(local_proportions$STDCount[local_proportions$Class == class])
    study_count <- nlcd_study_prop$STDCount[nlcd_study_prop$Class == class]
    
    # Perform chi-square test
    chisq_test <- chisq.test(c(hotspot_count, study_count))
    
    # Store the results
    chisqr_results <- rbind(chisqr_results, 
                            data.frame(Class = class, 
                                       ChiSquared = chisq_test$statistic,
                                       pvalue = chisq_test$p.value))
  }
  
  # View results
  view(chisqr_results)
  
# ---- 2.2: NLCD Statistics (Global) ----
  
  # Clip NLCD with the multipolygon
  # hs_nlcd <- crop(nlcd, hs)
  # hs_nlcd <- mask(hs_nlcd, hs)
  # writeRaster(hs_nlcd, 'hs_nlcd.tif')
  
  # Optional - Read in Clipped NLCD
  hs_nlcd = rast(here::here('data', 'nlcd_2021', 
                            'hs_nlcd.tif'))
  
  # Plot Clipped NLCD
  plot(hs_nlcd)
  
  # Calculate statistics
  nlcd_val <- values(hs_nlcd)
  proportions <- na.omit(data.frame(nlcd_val)) %>%
    group_by(Red) %>% # May be "Red"
    summarise(count = n()) %>%
    mutate(proportion = (count / sum(count)) * 100)
  
  # NLCD legend mapping
  nlcd_legend <- data.frame(
    Red = c(11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95),
    Class = c("Open Water", "Perennial Ice/Snow", "Developed, Open Space", "Developed, Low Intensity", 
              "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
              "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Dwarf Scrub", "Shrub/Scrub",
              "Herbaceous", "Sedge/Herbaceous", "Lichens", "Moss", "Hay/Pasture", "Cultivated Crops",
              "Woody Wetlands", "Emergent Herbaceous Wetlands")
  )
  
  proportions <- left_join(proportions, nlcd_legend, 
                           by = "Red") %>%
    select(Class, count, proportion)
  
  print(proportions)
  
  