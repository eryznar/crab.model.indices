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

# Read in spatial layers
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

# Set coordinate reference system
ncrs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# Set directory on the Y:: drive to save models, outputs, etc.
dir <- "Y:/KOD_Research/Ryznar/Model-based indices/BAIRDI/"

# Load prediction grid
pred_grid <- readRDS(here::here(paste0(dir, "Data/EBS_bairdi_grid_5km_No_Land.rds")))

pred_grid %>% dplyr::select(Lon, Lat) %>% distinct() -> pg

# Plot EBS pred grid
ggplot(pg, aes(Lon, Lat))+
  geom_point(size = 0.10, color = "darkgrey")+
  theme_bw()+
  ggtitle("EBS prediction grid")+
  ylab("Latitude")+
  xlab("Longitude") -> EBS.pg

ggsave(plot = EBS.pg, "./BAIRDI/Figures/EBS_predgrid.png", height = 4, width =6, units = "in")

# Filter Tanner W pred grid, save and plot
pred_grid %>% filter(Lon < -166) -> pg.W

#saveRDS(pg.W, paste0(dir, "Data/West_predgrid.rda"))

ggplot(pg.W, aes(Lon, Lat))+
  geom_point(size = 0.10, color = "darkgrey")+
  theme_bw()+
  ggtitle("Tanner crab West prediction grid")+
  ylab("Latitude")+
  xlab("Longitude") -> W.pg

ggsave(plot = W.pg, "./BAIRDI/Figures/west_predgrid.png", height = 4, width =6, units = "in")

# Filter Tanner E pred grid, save and plot
pred_grid %>% filter(Lon > -166) -> pg.E

#saveRDS(pg.E, paste0(dir, "Data/East_predgrid.rda"))

ggplot(pg.E, aes(Lon, Lat))+
  geom_point(size = 0.10, color = "darkgrey")+
  theme_bw()+
  ggtitle("Tanner crab East prediction grid")+
  ylab("Latitude")+
  xlab("Longitude") -> E.pg

ggsave(plot = E.pg, "./BAIRDI/Figures/east_predgrid.png", height = 4, width =6, units = "in")

### PROCESS DATA -----------------------------------------------------------------

# Load and process response data
tan.cpue <- read.csv(paste0(dir, "Data/bairdi_cpue2.csv")) %>%
  dplyr::select(!X) %>%
  st_as_sf(., coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  rename(year = AKFIN_SURVEY_YEAR, lat = Y, lon = X,
         cpue_km = CPUE_KM, cpue_kg_km = CPUE_KG_KM, stock = STOCK, cpue = CPUE, cpue_kg = CPUE_KG, matsex = MAT_SEX) %>%
  filter(matsex %in% c("Immature Female", "Immature Male", "Mature Male", "Mature Female")) %>%
  mutate(lat = lat/1000, # scale to km so values don't get too large
         lon = lon/1000) %>%
  dplyr::select(year, matsex, cpue, cpue_kg, cpue_km, cpue_kg_km, stock, lon, lat)

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
