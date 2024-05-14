# ---- 0: Load Packages ----

require(sf)
require(tidyverse)

# ---- 1: Prepare Data ----

hs_fall = st_read(here::here('SS_Runs', 'ss_2024_Pull', 'ss_1234rm_seasonal23_r125m', 
                             'devrm', 'hs_fall23_devrm.shp'))
hs_spring = st_read(here::here('SS_Runs', 'ss_2024_Pull', 'ss_1234rm_seasonal23_r125m', 
                               'devrm', 'hs_spring23_devrm.shp'))
hs_summer = st_read(here::here('SS_Runs', 'ss_2024_Pull', 'ss_1234rm_seasonal23_r125m', 
                               'devrm', 'hs_summer23_devrm.shp'))

trails = st_read(here::here('data', 'MAD_Trails', 'MAD_TRAILS.shp'))
trails = trails[!is.na(trails$TRAIL_NAME), ]

trails = st_transform(trails, st_crs(hs))


# ---- 2: Determine Nearest Trailhead ----

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

# Apply the function
spring_trails <- nearest_trailhead(hs_spring, trails)
summer_trails <- nearest_trailhead(hs_summer, trails)
fall_trails   <- nearest_trailhead(hs_fall, trails)

# Filter for trails within 500m of hs
spring_trails_500m = spring_trails %>% filter(nearest_distance < units::set_units(500, m))
summer_trails_500m = summer_trails %>% filter(nearest_distance < units::set_units(500, m))
fall_trails_500m   = fall_trails %>% filter(nearest_distance < units::set_units(500, m))

# Merge output for simple viewing
spring_trails_500m$Season = "SPRING"
summer_trails_500m$Season = "SUMMER"
fall_trails_500m$Season   = "FALL"
seasonal_hs_trails = rbind(spring_trails_500m, 
                           summer_trails_500m,
                           fall_trails_500m)

# Write the output
write.csv(seasonal_hs_trails, 'seasonal_hs_trails.csv')
st_write(seasonal_hs_trails, 'seasonal_hs_trails.shp')

display = select(seasonal_hs_trails, P_VALUE, nearest_trail_name, nearest_distance, Season)
