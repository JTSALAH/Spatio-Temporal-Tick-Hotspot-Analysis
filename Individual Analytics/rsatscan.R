# ---- 0: Load Packages & Data ----

  require(rsatscan)
  require(sf)
  require(ggplot2)
  require(rnaturalearth)
  require(rnaturalearthdata)

# ---- 1: Set Parameters ----

  # Reset Parameter File
  invisible(ss.options(reset=TRUE))
  
  # Set Parameters
  ss.options(list(CaseFile="NYCfever.cas", PrecisionCaseTimes=3))
  ss.options(c("StartDate=2001/11/1","EndDate=2001/11/24"))
  ss.options(list(CoordinatesFile="NYCfever.geo", 
                  AnalysisType=4, 
                  ModelType=2, 
                  TimeAggregationUnits=3))
  ss.options(list(UseDistanceFromCenterOption="y", 
                  MaxSpatialSizeInDistanceFromCenter=3, 
                  NonCompactnessPenalty=0))
  ss.options(list(MaxTemporalSizeInterpretation=1, 
                  MaxTemporalSize=7))
  ss.options(list(ProspectiveStartDate="2001/11/24", 
                  ReportGiniClusters="n", 
                  LogRunToHistoryFile="n"))
  ss.options(list(SaveSimLLRsDBase="y"))
  
  # Inspect Parameter File
  head(ss.options(),3)
  
  # Write Parameter File
  dir.create("ss")
  ss_folder = here::here("ss")
  write.ss.prm(ss_folder, "NYCfever")           # SatScan Parameter File
  write.cas(NYCfevercas, ss_folder, "NYCfever") # Case File
  write.geo(NYCfevergeo, ss_folder, "NYCfever") # Geo File

# ---- 2: Run SatScan ----
  
  # Input SatScan Program File Location
  sslocation = "C:/Program Files/SaTScan"
  
  # Run SatScan
  NYCfever = satscan(ss_folder, 
                     "NYCfever", 
                     sslocation = sslocation, 
                     ssbatchfilename = "SaTScanBatch64",
                     verbose = TRUE)
  
  # View Summary
  summary(NYCfever)
  
# ---- 3: Plot Output ----
  
  # Load basemap
  world <- ne_countries(scale = "medium", returnclass = "sf")

  # Plot Hot-Spot Clusters
  ggplot() +
    geom_sf(data = world) +
    geom_sf(data = NYCfever$shapeclust, color = "red", size = 3) +
    coord_sf(xlim = c(min(st_bbox(NYCfever$shapeclust)[c("xmin", "xmax")]), 
                      max(st_bbox(NYCfever$shapeclust)[c("xmin", "xmax")])), 
             ylim = c(min(st_bbox(NYCfever$shapeclust)[c("ymin", "ymax")]), 
                      max(st_bbox(NYCfever$shapeclust)[c("ymin", "ymax")]))) +
    theme_dark()
  