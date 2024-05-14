# ---- 0: Load Packages & Data ----

# WARNING: Clearing your R environment will delete the 'ssenv' object 
#          created when loading rsatscan and cause an error
require(rsatscan)
require(sf)
require(tidyverse)
require(rnaturalearth)
require(rnaturalearthdata)

# ---- 1: Prepare Input Files ----

# CasFile Requires: locationid, numcases
casfile = read.csv(here::here("ss_case.csv"))
geofile = read.csv(here::here("ss_geo.csv"))

# Subset to year of choice
casfile <- casfile %>%
  filter(year(mdy_hms(date)) == 2022) # 2023 Data Only
geofile <- geofile %>%
  filter(X %in% casfile$X)

# Get Date Ranges
maxdate <- casfile$date %>%
  as.Date(format = "%m/%d/%Y %H:%M:%S") %>%
  max() %>%
  format("%Y/%m/%d") %>%
  gsub("/0", "/", .) 
maxdate = paste("EndDate=", maxdate, sep="")

mindate <- casfile$date %>%
  as.Date(format = "%m/%d/%Y %H:%M:%S") %>%
  min() %>%
  format("%Y/%m/%d") %>%
  gsub("/0", "/", .)
mindate = paste("StartDate=", mindate, sep="")

# Write Input Files
dir.create("ss")
ss_folder = here::here("ss")
write.cas(casfile, ss_folder, "CasFile")  # Case File
write.geo(geofile, ss_folder, "GeoFile")  # Geo File

# Assign Input Files
casfile = paste0(ss_folder, "/CasFile.cas")
geofile = paste0(ss_folder, "/GeoFile.geo")


# ---- 2: Set Parameters ----

# Reset Parameter File
invisible(ss.options(reset=TRUE))

# Set Input Parameters
ss.options(list(CaseFile =        casfile, 
                # ControlFile =     "CtlFile.ctl",
                CoordinatesFile = geofile,
                CoordinatesType = 1      # Latitude/Longitude 
))

# Set Analysis Parameters
ss.options(list(# For the Space-Time Permutation model, the analysis type must be either Retrospective or Prospective Space-Time.
  AnalysisType = 3, # 3 = Retrospective Space-Time, 7 = Seasonal Temporal
  ModelType = 2,    # 1 = Bernoulli, 2 = Space-Time Permutation
  ScanAreas = 3     # Both High & Low Rates Detect
))
ss.options(list(# Try none, generic, and month
  PrecisionCaseTimes = 3,    # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
  TimeAggregationUnits = 3,  # 3 = Day
  TimeAggregationLength = 7, # 7 days = 1wk
  MinimumTemporalClusterSize = 7,  # Minimum size of 7 days
  MaxTemporalSize = 30,  # Maximum size of 30 days
  MaxSpatialSizeInPopulationAtRisk = 10, # 10% alt (50% default)
  CriteriaForReportingSecondaryClusters=0,
  MonteCarloReps = 99 # 99 rep minimum
))
ss.options(c(mindate, maxdate))
ss.options(c("UseDistanceFromCenterOption=y",
             "MaxSpatialSizeInDistanceFromCenter=0.125"))


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

# Inspect Parameter File
head(ss.options(),3)

# Write Parameter File
write.ss.prm(ss_folder, "Parameters")

# ---- 3: Run SatScan ----

# Input SatScan Program File Location
sslocation = "C:/Program Files/SaTScan"

# Run SatScan
satscan = satscan(ss_folder, 
                  "Parameters", 
                  sslocation = sslocation, 
                  ssbatchfilename = "SaTScanBatch64",
                  verbose = TRUE)

# Save SatScan Outputs
ss_layers = paste0(ss_folder, "/ss_layers")
dir.create(ss_layers)
save(satscan, file = paste0(ss_layers, "/satscan.rda"))
st_write(satscan$shapeclust, paste0(ss_layers, "/shapeclust.shp"))

# View Summary
summary(satscan)

# ---- 4: Plot Output ----

# Load basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot Hot-Spot Clusters
ggplot() +
  geom_sf(data = world) +
  geom_sf(data = satscan$shapeclust, color = "red", size = 3) +
  coord_sf(xlim = c(min(st_bbox(satscan$shapeclust)[c("xmin", "xmax")]), 
                    max(st_bbox(satscan$shapeclust)[c("xmin", "xmax")])), 
           ylim = c(min(st_bbox(satscan$shapeclust)[c("ymin", "ymax")]), 
                    max(st_bbox(satscan$shapeclust)[c("ymin", "ymax")]))) +
  theme_dark()
