### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of biomass for EBS snow crab for males >95mm
# and mature females. Time range is 1980-present.

# Author: Emily Ryznar

# TO DOs:
# 1) MAKE SURE PACKAGE VERSIONS (sdmTMB, glmmTMB, Matrix, TMB...) ARE THE SAME BETWEEN DESKTOP AND VM!! 

### LOAD LIBRARIES/PARAMS --------------------------------------------------------

library(INLA)
library(tidyverse)
library(sdmTMB)
library(glmmTMB)
library(broom)
library(sf)
library(gstat)
library(rnaturalearth)
library(raster)
library(concaveman)
library(png)
library(crabpack)
library(data.table)
library(future)
library(future.apply)

# Read in spatial layers
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

# Set coordinate reference system
ncrs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# Set directory on the Y:: drive to save models, outputs, etc.
dir <- "Y:/KOD_Research/Ryznar/Model-based indices/SNOW/"

# Load prediction grid
ebsnbs_grid <- read.csv(here::here(paste0(dir, "Data/bering_coarse_grid.csv")))
ebs_grid <- read.csv(here::here(paste0(dir, "Data/ebs_coarse_grid.csv")))
nbs_grid <- read.csv(here::here(paste0(dir, "Data/nbs_coarse_grid.csv")))


# Plot EBS pred grid
ggplot(ebs_grid, aes(X, Y))+
  geom_point(size = 1, color = "darkgrey")+
  theme_bw()+
  ggtitle("EBS prediction grid")+
  ylab("Latitude")+
  xlab("Longitude")

ggsave("./SNOW/Figures/EBS_predgrid.png", height = 4, width =6, units = "in")


### PROCESS DATA -----------------------------------------------------------------

# Load and process response data
snow.matfem.cpue <- readRDS(paste0(dir, "Data/snow_survey_cpue_matfem_EBSNBS.rda")) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(cpue_km = CPUE/3.429904, # num crab per nmi2 to km2
         cpue_kg_km = (CPUE_MT * 1000)/3.429904) %>% # crab kg per km2
  dplyr::rename(year = YEAR, lat = Y, lon = X,
         cpue = CPUE, cpue_mt = CPUE_MT, category = CATEGORY, region = REGION) %>%
  mutate(lat = lat/1000, # scale to km so values don't get too large
         lon = lon/1000) %>%
  dplyr::select(year, category, region, cpue, cpue_mt, cpue_km, cpue_kg_km, lon, lat) %>%
  mutate(category = "Mature female")

snow.male95.cpue <- readRDS(paste0(dir, "Data/snow_survey_cpue_male_EBSNBS.rda")) %>%
  filter(SIZE_1MM >95) %>%
  group_by(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE, SEX_TEXT, DISTRICT, STRATUM,
           TOTAL_AREA) %>%
  reframe(COUNT = sum(COUNT),
          CPUE = sum(CPUE),
          CPUE_MT = sum(CPUE_MT),
          CPUE_LBS = sum(CPUE_LBS)) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(cpue_km = CPUE/3.429904, # num crab per nmi2 to km2
         cpue_kg_km = (CPUE_MT * 1000)/3.429904) %>% # crab kg per km2
  dplyr::rename(year = YEAR, lat = Y, lon = X,
                cpue = CPUE, cpue_mt = CPUE_MT, region = REGION) %>%
  mutate(lat = lat/1000, # scale to km so values don't get too large
         lon = lon/1000) %>%
  dplyr::select(year, region, cpue, cpue_mt, cpue_km, cpue_kg_km, lon, lat) %>%
  mutate(category = "Male95")


# Load and process survey data
m.surv <- readRDS(paste0(dir, "Data/snow_survey_biomass_male_EBSNBS.rda")) %>%
            filter(SIZE_1MM >95) %>%
        group_by(YEAR, REGION) %>%
        reframe(BIOMASS_MT = sum(BIOMASS_MT),
                BIOMASS_MT_CV = sum(BIOMASS_MT_CV),
                BIOMASS_MT_CI = sum(BIOMASS_MT_CI)) %>%
        mutate(category = "Male95",
               BIOMASS_KG = BIOMASS_MT,
               BIOMASS_KG_CV = BIOMASS_MT_CV,
               BIOMASS_KG_CI = BIOMASS_MT_CI)

mf.surv <- readRDS(paste0(dir, "Data/snow_survey_biomass_matfem_EBSNBS.rda")) %>%
        mutate(category = "Mature female",
               BIOMASS_KG = BIOMASS_MT,
               BIOMASS_KG_CV = BIOMASS_MT_CV,
               BIOMASS_KG_CI = BIOMASS_MT_CI) %>%
        dplyr::select(YEAR, REGION, BIOMASS_MT, BIOMASS_MT_CV, BIOMASS_MT_CI, category, BIOMASS_KG, BIOMASS_KG_CV,
                      BIOMASS_KG_CI)
