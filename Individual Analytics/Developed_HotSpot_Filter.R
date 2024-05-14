# ---- 0: Load Packages ----

require(sf)
require(terra)
require(tidyverse)

# ---- 1: Prepare Data ----

  hs_full = st_read(here::here('SS_Runs', 'ss_2024_Pull', 'ss_1234rm_seasonal23_r125m', 
                               'shp', 'satscan_Summer_999mcreps.shp'))
  # Optional - Read in Clipped NLCD
  hs_nlcd = rast(here::here('data', 'nlcd_2021', 
                            'nlcd_wgs84.tif'))
  
  
  # ---- 2: Run HotSpot NLCD Analysis ----
  # Setup HotSpot Analysis Dependencies
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

# ---- 3: Calculate HotSpot NLCD Proportions ----
  
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
        left_join(nlcd_legend, by = "Red") %>%
        select(Class, Count, Proportion)
      prop_single$CLUSTER <- i
      local_proportions <- rbind(local_proportions, prop_single)
    }
    return(local_proportions)
  }
  
  # Run Proportion Calculations
  hs_local_proportions  = calculate_local_proportions(hs_full, hs_nlcd, nlcd_legend)

# ---- 4: Remove Developed Land Class HotSpots ----

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
  
  hs_local_proportions_filt = rm_developed(hs_local_proportions)
  hs_filt = unique(hs_local_proportions_filt$CLUSTER)
  hs_full_filt = hs_full %>% filter(CLUSTER %in% hs_filt)
  st_write(hs_full_filt,
           here::here('SS_Runs', 'ss_2024_Pull', 'ss_1234rm_seasonal23_r125m', 'hs_summer23_devrm.shp'),
    delete_layer = TRUE)
