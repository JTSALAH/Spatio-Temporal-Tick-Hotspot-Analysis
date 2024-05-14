# SEHTI SatScan Master Analysis
# INSTRUCTIONS: Specify all input data and preferred settings in Part 1, and run all subsequent code!

# ---- 0: Load in Packages ----

require(rsatscan)
require(sf)
require(terra)
require(tidyverse)
require(adehabitatLT)
require(here)

# ---- 1: Input Data ----

# 2: Pre-Process Data
gpx_folder_path = here('data', '2024_Data_Pull', "All_gpx")
Study_Dates     = read.csv(here('data', 'SpatialEpiOfTBD-CasesAndContacts_DATA_2024-02-09_1452.csv'))
counties        = st_read(here('data', 'counties', 'COUNTIES_POLY.shp'))
waterbodies     = st_read(here('data', 'majorhydro', 'MAJPOND_POLY.shp'))
roads           = st_read(here('data', 'MassDOT_Roads_SHP', 'EOTROADS_ARC.shp'))
structures      = st_read(here('data', 'SP_Buf10m_Study', 'SP_Buf10m_Study_wgs84.shp'))
# Optional - Read in Buffered Roads (TIME SAVE)
roads_buffered  = st_read(here('data', 'roads_buffered', 'roads_buffered.shp'))

# 3: Seasonal SatScan
sslocation = "C:/Program Files/SaTScan" # Ensure this is consistent with your machine
# Choose Case Type to Filter by Restricted Case (RC) & General Case (GC), or don't!
CaseType = "RC" # Options: RC, GC, All
# Desired ratio of Cases to Controls
ControltoCase_Ratio <- 2 # For example, 1 is a 1:1 ratio, 2 is a 2:1 ratio, 0 means no controls used
year = c(2022, 2023)     # You can specify each year you want, or "All" to use all years
seasons <- c("Fall", "Winter")

# 4: Filter Developed Areas
hs_nlcd = rast(here('data', 'nlcd_2021', 'nlcd_wgs84.tif'))
nlcd_study_val = read.csv(here::here('data', 'nlcd_2021','NLCD_Study_Values.csv'))

# 5: Determine Nearest Trailheads
trails = st_read(here::here('data', 'MAD_Trails', 'MAD_TRAILS.shp'))
trails = trails[!is.na(trails$TRAIL_NAME), ]

# ---- 2: Pre-Process Data ----

# 2.1: Load GPX Data From Folder
read_gpx_folder = function(folder_path) {
  # List all GPX files in the folder
  gpx_files = list.files(path = folder_path, pattern = "\\.gpx$", full.names = TRUE)
  
  # Extract just the file names without the .gpx extension
  file_names = sapply(gpx_files, function(file) {
    file_name = basename(file)
    sub("\\.gpx$", "", file_name)
  })
  
  # Read each file into an sf object and store them in a list
  sf_list = setNames(lapply(gpx_files, function(file) {
    sf_object = tryCatch({
      st_read(file, layer = "track_points", quiet = TRUE)
    }, error = function(e) {
      warning(paste("Error reading file:", file, "\n", e))
      return(NULL)
    })
    
    if (is.null(sf_object)) {
      return(NULL)
    }
    
    coords = st_coordinates(sf_object)
    sf_object$Longitude = coords[, "X"]
    sf_object$Latitude = coords[, "Y"]
    return(sf_object)
  }), file_names)
  
  # Remove NULL elements
  sf_list = sf_list[!sapply(sf_list, is.null)]
  
  if (length(sf_list) == 0) {
    stop("No valid GPX files found in the directory.")
  }
  
  # Standardize columns in all sf objects
  all_cols = unique(unlist(lapply(sf_list, names)))
  sf_list = lapply(sf_list, function(sf) {
    missing_cols = setdiff(all_cols, names(sf))
    for (col in missing_cols) {
      sf[[col]] = NA
    }
    sf[all_cols]
  })
  
  # Combine all sf objects into a single sf object
  combined_sf = do.call(rbind, sf_list)
  
  # Convert the combined sf object to a dataframe
  combined_df = as.data.frame(st_drop_geometry(combined_sf))
  
  # Return both the combined dataframe and the list of individual sf objects
  return(list(combined_df = combined_df, sf_list = sf_list))
}

# Extract Data from Folder
length(list.files(gpx_folder_path))
gpx_folder = read_gpx_folder(gpx_folder_path)
gpx_list = gpx_folder$sf_list
track_all = gpx_folder$combined_df

# 2.2: Isolate Trajectory Movement Points
# 2.2.1: Batch Process Individual Trajectories
batch_ltraj = function(gpx_list) {
  # Apply as.ltraj to each individual
  ltraj_list = lapply(names(gpx_list), function(name) {
    # Load Dataframe
    df = gpx_list[[name]]
    df = st_drop_geometry(df)
    
    # Quality Control Dataframe
    duplicates = duplicated(df$time)
    df = df[!duplicated(df$time), ]
    df = df[complete.cases(df$Longitude, df$Latitude, df$time), ]
    
    # Run Trajectory Function
    as.ltraj(xy = df[, c("Longitude", "Latitude")], 
             date = df$time, 
             id = name, 
             typeII = TRUE)
  })
  
  # Set the names of the ltraj_list to be the same as those of gpx_list
  names(ltraj_list) = names(gpx_list)
  
  return(ltraj_list)
}
ltraj_list = batch_ltraj(gpx_list)

# 2.2.2: Convert Trajectories to Dataframes & Combine
ltraj_df_list = lapply(names(ltraj_list), function(name) {
  # Access the ltraj object in the nested structure
  ltraj_obj = ltraj_list[[name]][[1]]
  # Convert the ltraj object to a dataframe
  df = as.data.frame(ltraj_obj)
  # Add the ID column
  df$ID = name
  return(df)
})
PTS = do.call(rbind, ltraj_df_list)

# Remove Points w/ Long Idle Times
PTS_Filt = PTS %>% 
  filter(dt < 900)

PTS_Filt_SF = PTS %>% 
  filter(dt < 900) %>%    
  st_as_sf(coords = c("x", "y"), 
           crs = 4326)

# 2.3: Filter Data by Study Dates
# Wrangle Study Dates Dataframe
Study_Dates = Study_Dates %>%
  mutate(
    datacollection_startdate = as.POSIXct(datacollection_startdate, format = "%Y-%m-%d"),
    datacollection_enddate   = as.POSIXct(datacollection_enddate,   format = "%Y-%m-%d")
  )
Study_Dates$record_id <- sprintf("%03d", Study_Dates$record_id)
PTS_Filt_SF$record_id <- substr(PTS_Filt_SF$ID, 1, 3)

# Determine if Participant is a Case or Control
Study_Dates <- Study_Dates %>%
  mutate(status = case_when(
    (screen_q4_case == 1 | screen_q4_case_2 == 1 | feedback_tickexp == 1 | tickexp_yn == 1) ~ "Case",
    TRUE ~ "Control"),
    CaseType = case_when(
      screen_q4_case == 1 ~ "RC",
      screen_q4_case == 0 | screen_q4_case_2 == 1 | feedback_tickexp == 1 | tickexp_yn == 1 ~ "GC",
      TRUE ~ NA_character_))

# Join the Dataframes
PTS_Filt_SF <- left_join(PTS_Filt_SF, Study_Dates, by = "record_id")

# Inspect the dataset
table(PTS_Filt_SF$CaseType)
PTS_Filt_SF %>%
  filter(CaseType %in% c("RC", "GC")) %>%  # Filter for RC and GC CaseTypes
  group_by(CaseType) %>%                   # Group by CaseType
  summarise(Unique_IDs = n_distinct(ID))  

# Filter based on the time criteria
PTS_Study <- PTS_Filt_SF %>%
  filter(date >= datacollection_startdate & date <= datacollection_enddate)

# ---- 3: Run multiple iterations of the analysis ----

source(here('SEHTI_SS.R'))

num_iters = 5
for (i in 1:num_iters) {
  sehti_ss(PTS_Study, Study_Dates, 
           counties, waterbodies, roads, structures, roads_buffered,
           hs_nlcd, nlcd_study_val, trails,
           CaseType, ControltoCase_Ratio, year, seasons,
           loop_iteration = i)
}
