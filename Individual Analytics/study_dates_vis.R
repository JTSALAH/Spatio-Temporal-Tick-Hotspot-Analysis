# ---- 0: Load Packages ----

  require(tidyverse)

# ---- 1: Prepare Data ----

  pts = read.csv(file.choose())
  pts$ID <- substr(pts$ID, 1, 3)
  class(pts$time)
  
  # dates = read.csv(file.choose())
  dates = read.csv(here::here('data', 'SpatialEpiOfTBD-CasesAndContacts_DATA_2023-12-12_1134.csv'))
  dates$record_id <- sprintf("%03d", dates$record_id)

# ---- 2: Filter ----

  pts$time <- as.Date(mdy_hms(pts$time))
  
  # Ensure the start and end dates in dates are in Date format
  dates$datacollection_startdate <- as.Date(dates$datacollection_startdate)
  dates$datacollection_enddate <- as.Date(dates$datacollection_enddate)
  
  # Join the dataframes
  combined_df <- left_join(pts, dates, by = c("ID" = "record_id"))
  
  # Filter based on the time criteria
  pts_study <- combined_df %>%
    filter(time >= datacollection_startdate & time <= datacollection_enddate)
  
  write.csv(pts_study, "TP_WMA_RD1234RM_Study.csv")

# ---- 3: Number of PTS per Week ----

  pts_study = read.csv("TP_WMA_RD1234RM_Study.csv")
  
  # Convert 'character' class dates to 'POSIXct' date class
  pts_study$date <- mdy_hms(pts_study$date)
  class(pts_study$date)
  
  # Create a formatted Year-Week label
  pts_study$year_week <- paste0("'", year(pts_study$date) %% 100, " wk", week(pts_study$date))
  class(pts_study$year_week)
  
  # Plot a histogram
  ggplot(pts_study, aes(x = year_week)) +
    geom_bar(stat = "count") +  # Using geom_bar for categorical data
    labs(x = "Year-Week", y = "Count of Data Points", title = "Data Points by Week in SEHTI Study") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 0.5))
