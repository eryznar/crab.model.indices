### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

#PLOTS: 1) model 100 between knots 2) QQ plots and spatial residuals 3) index
# other different distribution, ar1 vs. iid, joining timeseries together, depth

# TO DOs:
# 1) Look at residuals
# 2) Add in scripts to load new survey data (CPUE, BIO/ABUND) and process each year (CPUE script is in TECHMEMONEW)

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
ggplot(ebsnbs_grid, aes(X, Y))+
  geom_point(size = 1, color = "darkgrey")+
  theme_bw()+
  ggtitle("EBS/NBS prediction grid")+
  ylab("Latitude")+
  xlab("Longitude")

#ggsave(plot = EBS.pg, "./BAIRDI/Figures/EBS_predgrid.png", height = 4, width =6, units = "in")


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
  dplyr::select(year, category, region, cpue, cpue_mt, cpue_km, cpue_kg_km, lon, lat)

# Create dummy 2020 data
# tan.cpue %>% 
#   filter(year == 2024) %>%
#   mutate(year = 2020) -> dummy
# 
# rbind(tan.cpue, dummy) -> tan.cpue

# group male categories
tan.cpue2 <- tan.cpue %>%
  mutate(matsex = case_when((matsex %in% c("Immature Male", "Mature Male")) ~ "Male",
                            TRUE ~ matsex)) %>%
  group_by(year, lon, lat, matsex, stock) %>%
  reframe(cpue = sum(cpue),
          cpue_kg = sum(cpue_kg),
          cpue_km = sum(cpue_km),
          cpue_kg_km = sum(cpue_kg_km)) %>%
  rename(mat.sex = matsex, stck = stock)


# Load observed abundance/biomass
tan.obs <- right_join(rbind(read.csv(paste0(dir, "Data/E166_CB_OBSERVEDabundbio.csv")),
                            read.csv(paste0(dir, "Data/W166_CB_OBSERVEDabundbio.csv"))) %>%
                        rename(Year = AKFIN_SURVEY_YEAR, matsex = MAT_SEX, stock = STOCK, abundance= ABUNDANCE, biomass = BIOMASS) %>%
                        dplyr::select(Year, matsex, abundance, biomass, stock) %>%
                        mutate(abundance = abundance/1e6, biomass = biomass/1000) %>%
                        pivot_longer(., c("abundance", "biomass"), names_to = "type", values_to = "value"),
                      rbind(read.csv(paste0(dir, "Data/E166_CB_OBSERVEDabundbio.csv")),
                            read.csv(paste0(dir, "Data/W166_CB_OBSERVEDabundbio.csv"))) %>%
                        dplyr::select(AKFIN_SURVEY_YEAR,MAT_SEX, ABUNDANCE_CI, BIOMASS_CI, STOCK) %>%
                        rename(Year = AKFIN_SURVEY_YEAR, matsex = MAT_SEX, stock = STOCK, abundance = ABUNDANCE_CI, biomass = BIOMASS_CI) %>%
                        mutate(abundance = abundance/1e6, biomass = biomass/1000) %>%
                        pivot_longer(., c("abundance", "biomass"), names_to = "type", values_to = "CI")) %>%
  filter(matsex %in% c("Immature Female", "Mature Female", "Immature Male", "Mature Male")) %>%
  mutate(matsex = case_when(matsex %in% c("Immature Male", "Mature Male") ~ "Male",
                            TRUE ~ matsex)) %>%
  group_by(Year, type, matsex, stock) %>%
  reframe(value = sum(value),
          CI = sum(CI))


# ggplot(tan.obs %>% filter(stock == "TannerW", type == "abundance"), aes(Year, value))+
#   geom_line()+
#   facet_wrap(~matsex, scales = "free_y")
