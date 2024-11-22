### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Extend year range back to 1975
# 2) Limit CW to >=25mm

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

# Read in spatial layers
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

# Set coordinate reference system
ncrs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# polyon for cropping objects
ebs <- st_read(here::here("BAIRDI/data/ebs_friendly.kml")) %>% 
  st_transform(., crs = ncrs)

### PROCESS DATA -----------------------------------------------------------------
# shoreline
shoreline <- 
  ne_states(country = "United States of America") %>% 
  filter(name == "Alaska") %>% 
  st_union() %>% 
  st_transform(., crs = ncrs) %>% 
  st_intersection(., ebs) 


sample <- read_csv(file = here::here("BAIRDI/data/bairdi_cpue.csv")) %>% 
  # filter(AKFIN_SURVEY_YEAR == 2024) %>% 
  st_as_sf(., coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = in.crs) %>%
  st_transform(., crs = ncrs) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  rename(year = AKFIN_SURVEY_YEAR, lat = Y, lon = X,
       cpue = CPUE, cpue_kg = CPUE_KG) %>%
  mutate(lon_group = factor(
    cut_number(lon, n = 10))) %>%
  filter(MAT_SEX %in% c("Immature Female", "Immature Male", "Mature Male", "Mature Female"))

males <- sample %>%
          filter(MAT_SEX %in% c("Immature Male", "Mature Male")) %>%
          group_by(GIS_STATION, REGION, AREA_SWEPT, HAUL_TYPE, lon, lat, geometry, lon_group) %>%
          reframe(cpue = sum(cpue),
                  cpue_kg = sum(cpue_kg)) %>%
          mutate(MAT_SEX = "All males")


imfem <- sample %>%
          filter(MAT_SEX == "Immature Female")


matfem <- sample %>%
            filter(MAT_SEX == "Mature Female")

### CREATE MESH ------------------------------------------------------------------
  # create buffer region around survey locs to use in mesh
  survey_buffed <-
    sample %>% 
    st_as_sf(., coords = c('lon','lat'), crs = ncrs) %>% 
    concaveman::concaveman(.) %>%
    st_union() %>%
    st_buffer(dist = 15) %>%
    {. ->> survey_buffed_sf} %>% 
    as_Spatial()
  
  plot(survey_buffed)
  
  # rasterize the polygon to create the prediction grid
  r <- raster(extent(survey_buffed),
              res = 10,
              crs = crs(ncrs))
  
  # prediction grid 
  pred_grid <- 
    fasterize::fasterize(survey_buffed_sf, r) %>% 
    as("SpatialPixelsDataFrame") %>%
    as_tibble() %>% 
    dplyr::select(lon = x, 
                  lat = y)
  
  ## INLA mesh----
  mesh2 <- INLA::inla.mesh.2d(boundary = survey_buffed,
                             max.edge = c(20, 50), # max allowed triangle edge length inside and outside boundary
                             offset = c(10, 40), # extension of inner and outer boundary
                             cutoff = 25) #minimum distance between points (was 20)
  plot(mesh2)
  
  ## convert for sdmTMB----
  ebs_spde <- sdmTMB::make_mesh(data = sample,
                                xy_cols = c("lon", "lat"),
                                mesh = mesh2)
  plot(ebs_spde)

### FIT MODELS -------------------------------------------------------------------
# 1) All Males

# Tweedie GLMM w/ spatial RE-----


# fit a spatial GLMM
m3 <- sdmTMB(cpue_kg ~ 1, #the 0 is there so there is a factor predictor for each time slice
             spatial = "on",
             mesh = ebs_spde,
             family = tweedie(),
             data = sample)

# m3 <- sdmTMB(cpue_kg ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
#              spatial = "on",
#              mesh = ebs_spde,
#              family = tweedie(),
#              time = "year",
#              data = sample)

# evaluate residuals
# Look at residual patterning
sample$s_glmm_resids <- residuals(m3, type = "mle-mvn")

# visualize residuals across the EBS
ggplot(sample) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, x = lon, color = s_glmm_resids), size = 5) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16)) 


# visualize using a variogram
v3_out <- tibble()
for (i in 1:500){
  sample$s_glmm_resids <- residuals(m3, type = "mle-mvn")
  ex3_sp <- sample
  coordinates(ex3_sp) <- c("lon", "lat")

  v3 <- variogram(object = s_glmm_resids ~ 1,
                  data = ex3_sp,
                  cressie = T,
                  cutoff = 500) %>%
    mutate(draw = i)
  assign("v3_out", rbind(v3, v3_out))
}

ggplot(v3_out) +
  geom_boxplot(aes(y = gamma,
                   x = dist,
                   group = dist)) +
  labs(y = "Sample Variogram",
       x = "Distance (km)") +
  theme_gray() +
  theme(axis.title = element_text(size = 15))

# predict from the model to extract spatial random field (LOG SCALE)
pred_df <- predict(m3, newdata = pred_grid)

# model predictions (global intercept + spatial random field)
ggplot(pred_df) +
  geom_raster(aes(y = lat, x = lon, fill = est)) + 
  geom_sf(data = shoreline) +
  scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       fill = "Log bairdi per sq. km") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16),
        legend.position = "bottom")
