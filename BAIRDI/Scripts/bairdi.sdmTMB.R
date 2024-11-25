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

### LOAD FUNCTIONS -------------------------------------------------------------
eval_resid <- function(data, model, type, matsex){
  data$s_glmm_resids <- residuals(model, type = "mle-mvn")
  
  # visualize residuals across the EBS
  ggplot(data) + 
    geom_sf(data = shoreline) +
    geom_point(aes(y = lat, x = lon, color = s_glmm_resids), size = 1) +
    scale_color_gradient2(midpoint = 0) + 
    labs(y = "Latitude",
         x = "Longitude") +
    theme_gray() + 
    ggtitle(paste("Bairdi", matsex, type, "residuals"))+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom") -> res_plot
  
  ggsave(plot = res_plot, paste0("./BAIRDI/Figures/bairdi_", matsex, "_", type, "_resid.png"), height = 9, width = 8.5, units = "in")
  
  return(res_plot)
}

predict_model <- function(newdat, model, type, matsex, years){
  newdat %>%
    filter(year %in% years) -> newdat2
 
  out <- predict(model, newdata= newdat2, return_tmb_object = T)
    
  lab <- ifelse(type == "abundance", "Log bairdi per sq.km", "Log bairdi biomass (kg) per sq.skm")
  
  # model predictions (global interept + spatial random field)
  ggplot(out$data) +
    geom_raster(aes(y = lat, x = lon, fill = est)) + 
    geom_sf(data = shoreline) +
    scale_fill_gradient2() + 
    labs(y = "Latitude",
         x = "Longitude",
         fill = lab) +
    theme_gray() + 
    facet_wrap(~year)+
    ggtitle(paste("Bairdi", matsex, "predicted", type))+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom") -> pred_plot
  
  ggsave(plot = pred_plot, paste0("./BAIRDI/Figures/bairdi_", matsex, "_", type, "_predicted.png"), height = 9, width = 8.5, units = "in")
  
  
  return(list(pred = out, pred_plot = pred_plot))
}

### PROCESS DATA -----------------------------------------------------------------
# load observed total abundance/biomass
tanEbio<- read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/bio_EBS_TannerE.csv")
tanWbio <- read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/bio_EBS_TannerW.csv")

tanEabund<- read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/abund_EBS_TannerE.csv")
tanWabund <- read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/abund_EBS_TannerW.csv")

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
          group_by(year, GIS_STATION, REGION, AREA_SWEPT, HAUL_TYPE, lon, lat, geometry, lon_group) %>%
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
  
  ## INLA mesh
  mesh2 <- INLA::inla.mesh.2d(boundary = survey_buffed,
                             max.edge = c(20, 50), # max allowed triangle edge length inside and outside boundary
                             offset = c(10, 40), # extension of inner and outer boundary
                             cutoff = 25) #minimum distance between points (was 20)
  plot(mesh2)
  
  ## convert for sdmTMB
  male.mesh <- sdmTMB::make_mesh(data = males,
                                xy_cols = c("lon", "lat"),
                                mesh = mesh2)
 
  imfem.mesh <- sdmTMB::make_mesh(data = imfem,
                                 xy_cols = c("lon", "lat"),
                                 mesh = mesh2)
  
  matfem.mesh <- sdmTMB::make_mesh(data = matfem,
                                 xy_cols = c("lon", "lat"),
                                 mesh = mesh2)
  
### FIT MODELS -------------------------------------------------------------------

# 1) All Males ----

# Fit Tweedie GLMM w/ spatial RE
male.abund <- sdmTMB(cpue ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                     spatial = "on",
                     mesh = male.mesh,
                     family = tweedie(),
                     time = "year",
                     anisotropy = TRUE,
                     data = males)
  
saveRDS(male.abund, "./BAIRDI/bairdi_male_abundTMB.rda")

male.bio <- sdmTMB(cpue_kg ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
             spatial = "on",
             mesh = male.mesh,
             family = tweedie(),
             time = "year",
             anisotropy = TRUE,
             data = males)

saveRDS(male.abund, "./BAIRDI/bairdi_male_bioTMB.rda")



# Look at residuals
eval_resid(males, male.abund, "abundance", "male")
eval_resid(males, male.bio, "biomass", "male")

# set plotting years
years <- c(1988:2019, 2021:2024)

# Replicate pred_grid (newdata) by years
pred_grid2 <- replicate_df(pred_grid, "year", years)

# Predict from model
predict_model(pred_grid2, male.abund, "abundance", "male", years) -> pred.male.abund
predict_model(pred_grid2, male.bio, "biomass", "male", years) -> pred.male.bio

saveRDS(pred.male.abund, "./BAIRDI/Output/predicted_male_abundance(spatial).rda")
saveRDS(pred.male.bio, "./BAIRDI/Output/predicted_male_biomass(spatial).rda")

# Extract the total abundance/biomass calculations and standard errors
get_index(pred.male.abund$pred, area = 1, bias_correct = TRUE) -> male.abund.ind
get_index(pred.male.bio$pred, area = 1, bias_correct = TRUE) -> male.bio.ind

# Process observed abundance/biomass
rbind(tanEbio %>%
      dplyr::select(Year, Small.male...113.mm., Large.male...113.mm.) %>%
      rename(smmale = Small.male...113.mm., lgmale = Large.male...113.mm.) %>%
      pivot_longer(., cols = c(2:3), names_to = "matsex", values_to = "biomass") %>%
      mutate(CI = as.numeric(gsub(",", "", sapply(str_extract_all(biomass, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =","))),
             biomass = as.numeric(gsub(",", "", gsub("\\([^()]*\\)", "", biomass)))),
      tanWbio %>%
        dplyr::select(Year, Small.male...103.mm., Large.male...103.mm.) %>%
        rename(smmale = Small.male...103.mm., lgmale = Large.male...103.mm.) %>%
        pivot_longer(., cols = c(2:3), names_to = "matsex", values_to = "biomass") %>%
        mutate(CI = as.numeric(gsub(",", "", sapply(str_extract_all(biomass, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =","))),
               biomass = as.numeric(gsub(",", "", gsub("\\([^()]*\\)", "", biomass))))) %>%
  group_by(Year) %>%
  reframe(biomass = sum(biomass),
          CI = sum(CI)) -> male.bio.obs

rbind(tanEabund %>%
        dplyr::select(Year, Small.male...113.mm., Large.male...113.mm.) %>%
        rename(smmale = Small.male...113.mm., lgmale = Large.male...113.mm.) %>%
        pivot_longer(., cols = c(2:3), names_to = "matsex", values_to = "abundance") %>%
        mutate(CI = as.numeric(gsub(",", "", sapply(str_extract_all(abundance, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =","))),
               abundance = as.numeric(gsub(",", "", gsub("\\([^()]*\\)", "", abundance)))),
      tanWabund %>%
        dplyr::select(Year, Small.male...103.mm., Large.male...103.mm.) %>%
        rename(smmale = Small.male...103.mm., lgmale = Large.male...103.mm.) %>%
        pivot_longer(., cols = c(2:3), names_to = "matsex", values_to = "abundance") %>%
        mutate(CI = as.numeric(gsub(",", "", sapply(str_extract_all(abundance, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =","))),
               abundance = as.numeric(gsub(",", "", gsub("\\([^()]*\\)", "", abundance))))) %>%
  group_by(Year) %>%
  reframe(abundance = sum(abundance),
          CI = sum(CI)) -> male.abund.obs



# Plot predicted vs. observed
male.abund.ind %>% 
  mutate(est = est/1e6, lwr = lwr/1e6, upr = upr/1e6) %>% 
  rename(Year = year, abundance = est) -> male.abund.ind
  
ggplot()+
  geom_line(male.abund.ind, mapping = aes(Year, abundance))+
  #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
  geom_line(male.abund.obs, mapping = aes(Year))
                                                                                                            

