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

# shoreline
shoreline <- 
  ne_states(country = "United States of America") %>% 
  filter(name == "Alaska") %>% 
  st_union() %>% 
  st_transform(., crs = ncrs) %>% 
  st_intersection(., ebs) 

sample <- read_csv(file = here::here("BAIRDI/data/bairdi_cpue.csv")) %>% 
  filter(AKFIN_SURVEY_YEAR == 2010) %>% 
  st_as_sf(., coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = in.crs) %>%
  st_transform(., crs = ncrs) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  rename(year = AKFIN_SURVEY_YEAR, lat = Y, lon = X,
       cpue = CPUE) %>%
  mutate(lon_group = factor(
    cut_number(lon, n = 10))
  )

# Question 1: What is the mean density of Tanner crab in the EBS?

# Example 1: Tweedie GLM. We estimate crab densities using an intercept-only
# GLM

# Tweedie GLM----
m1 <- sdmTMB(cpue ~ 1, data = sample, family = tweedie(), spatial = "off")

# Extract predictions and SEs. In the intercept only model, the intercept is the 
# estimate of the mean density
glm_pred_df <- predict(m1, se_fit = T) %>% 
  mutate(glm_est = exp(est),
         glm_upr = exp(est + 1.96 * est_se),
         glm_lwr = exp(est - 1.96 * est_se)) %>% 
  dplyr::select(-1)

ggplot(glm_pred_df) + 
  geom_histogram(aes(cpue)) + 
  geom_vline(aes(xintercept = glm_est), color = "purple", linewidth = 1) + 
  geom_vline(aes(xintercept = glm_lwr), color = "purple", linetype = 2, linewidth = 1) +
  geom_vline(aes(xintercept = glm_upr), color = "purple", linetype = 2, linewidth = 1) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_minimal() + 
  labs(y = "Frequency",
       x = "Bairdi per sq. km") + 
  theme(axis.title = element_text(size = 16))

# Look at residual patterning
sample$glm_resids <- residuals(m1, type = "mle-mvn")

# visualize residuals across the EBS
ggplot(sample) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, x = lon, color = glm_resids), size = 5) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16))
  # coord_sf(ylim = c(550, 1000),
  #          xlim = c(-850, -250))
