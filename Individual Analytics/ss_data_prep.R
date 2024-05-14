# ---- 0: Load Packages ----

  require(sf)
  require(tidyverse)
  require(adehabitatLT)

# ---- 1: Load GPX Data From Folder ----

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
  gpx_folder_path = here::here('data', '2024_Data_Pull', "All_gpx")
  length(list.files(gpx_folder_path))
  gpx_folder = read_gpx_folder(gpx_folder_path)
  gpx_list = gpx_folder$sf_list
  track_all = gpx_folder$combined_df

# ---- 2: Isolate Trajectory Movement Points ----
  
  # 2.1: Batch Process Individual Trajectories
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
  
  # 2.2: Convert Trajectories to Dataframes & Combine
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
  
  # 2.2: Remove Points w/ Long Idle Times
  PTS_Filt = PTS %>% 
    filter(dt < 900)
  
  PTS_Filt_SF = PTS %>% 
    filter(dt < 900) %>%    
    st_as_sf(coords = c("x", "y"), 
             crs = 4326)
  
  # OPTIONAL: Write PTS Shapefile
  # st_write(PTS_Filt_SF, "PTS.shp")

# ---- 3: Filter Data by Study Dates ----
  
  # Read Study Dates Dataframe
  Study_Dates = read.csv(here::here('data', 'SpatialEpiOfTBD-CasesAndContacts_DATA_2024-02-09_1452.csv'))
  Study_Dates = Study_Dates %>%
    mutate(
      datacollection_startdate = as.POSIXct(datacollection_startdate, format = "%Y-%m-%d"),
      datacollection_enddate   = as.POSIXct(datacollection_enddate, format = "%Y-%m-%d")
    )
  Study_Dates$record_id <- sprintf("%03d", Study_Dates$record_id)
  PTS_Filt_SF$record_id <- substr(PTS_Filt_SF$ID, 1, 3)
  
  # Determine if Participant is a Case or Control
  Study_Dates <- Study_Dates %>%
    mutate(status = case_when(
      (screen_q4_case == 1 | screen_q4_case_2 == 1 | feedback_tickexp == 1 | tickexp_yn == 1) ~ "Case",
      TRUE ~ "Control"
    ))
  
  # Join the Dataframes
  PTS_Filt_SF <- left_join(PTS_Filt_SF, Study_Dates, by = "record_id")
  
  # Filter based on the time criteria
  PTS_Study <- PTS_Filt_SF %>%
    filter(date >= datacollection_startdate & date <= datacollection_enddate)
  
# ---- 4: Filter Data by Case vs. Control ----
  
  # Filter Case & Control
  PTS_Case    = PTS_Study %>% filter(status == "Case")
  PTS_Control = PTS_Study %>% filter(status == "Control")
  
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
  
  # Note: 
  
# ---- 5: Filter Study Area & Remove Houses, Roads, & Waterbodies ----
  
  # Read in Filter Data Layers
  counties    = st_read(here::here('data', 'counties', 'COUNTIES_POLY.shp'))
  waterbodies = st_read(here::here('data', 'majorhydro', 'MAJPOND_POLY.shp'))
  roads       = st_read(here::here('data', 'MassDOT_Roads_SHP', 'EOTROADS_ARC.shp'))
  structures  = st_read(here::here('data', 'SP_Buf10m_Study', 'SP_Buf10m_Study_wgs84.shp'))
  
  # Clip by county study area to reduce computational time
  counties    = counties %>% filter(COUNTY == c('HAMPSHIRE', 'HAMPDEN', 'FRANKLIN'))
  waterbodies = st_intersection(waterbodies, counties)
  roads       = st_intersection(roads, counties)
  
  # Ensure Data Layer Projections Match PTS
  counties    = st_transform(counties, crs = "EPSG:4326")
  waterbodies = st_transform(waterbodies, crs = "EPSG:4326")
  roads       = st_transform(roads, crs = "EPSG:4326")
  st_crs(structures) == st_crs(roads)
  
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
  
  # Apply the function to the sf object
  if (file.exists(here::here('data', 'roads_buffered', 'roads_buffered.shp'))){
    roads_buffered = st_read(here::here('data', 'roads_buffered', 'roads_buffered.shp'))
  } else {
    roads_buffered <- buffer_by_rdtype(roads_1234)
    st_write(roads_buffered, 'roads_buffered.shp')
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
  if(any(!st_is_valid(PTS_Filt_SF))) {
    PTS_Filt_SF$geometry[!st_is_valid(PTS_Filt_SF)] <- st_make_valid(PTS_Filt_SF$geometry[!st_is_valid(PTS_Filt_SF)])
  }
  
  # Clip to Study Area
  PTS_Study = st_intersection(PTS_Case, counties)
  
  # Efficient Erase
  PTS_Erase = st_intersects(PTS_Study, erase_sf)
  non_intersecting_points <- which(lengths(PTS_Erase) == 0)
  PTS_Final <- PTS_Study[non_intersecting_points, ]
  st_write(PTS_Final, 'PTS_Final.shp')
  
# ---- 6: Create Case & Geo Files ----
  
  # Optional - Read in Filtered CSV
  # PTS_Final = read.csv(file.choose())
  
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
  write.csv(ss_case, "ss_case.csv")
  write.csv(ss_geo, "ss_geo.csv")
