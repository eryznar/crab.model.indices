### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Look at residuals
# 2) Troubleshoot immature female abundance
# 3) Get abund/bio obs back to 1975

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

pred_grid <- readRDS(here::here("BAIRDI/data/EBS_bairdi_grid_5km_No_Land.rds")) 

### LOAD FUNCTIONS -------------------------------------------------------------
fit_models <- function(data, matsex, stock, years, period, knots){
  
  # Filter data by params
  if(stock != "All"){
    data %>%
      filter(mat.sex == matsex, stck == stock, year %in% years) -> data2
  }else{
    data %>%
      filter(mat.sex == matsex, year %in% years) -> data2
  }
  
  # Make mesh
  mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = knots, type = "kmeans")
  
  
  # Fit models
  print("Fitting abundance model")
  abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                  spatial = "on",
                  spatiotemporal = "iid",
                  mesh = mesh2,
                  family = tweedie(link = "log"),
                  time = "year",
                  anisotropy = TRUE,
                  data = data2)
  
  saveRDS(abund, paste0("./BAIRDI/Models/bairdi_", matsex, "_", stock, "_", period, "_", knots, "_abundTMB.rda"))
  
  print("Fitting biomass model")
  bio <- sdmTMB(cpue_kg_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                spatial = "on",
                spatiotemporal = "iid",
                mesh = mesh2,
                family = tweedie(link = "log"),
                time = "year",
                anisotropy = TRUE,
                data = data2)
  
  saveRDS(bio, paste0("./BAIRDI/Models/bairdi_", matsex, "_", stock, "_", period, "_", knots, "_bioTMB.rda"))
  
  return(list(abundTMB = abund, bioTMB = bio, mesh = mesh2))
}

eval_resid <- function(data, model, type, matsex){
  data$s_glmm_resids <- residuals(model, type = "mle-mvn")
  
  # visualize residuals across the EBS
  ggplot(data) + 
    #geom_sf(data = shoreline) +
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

predict_and_getindex <- function(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots){
  newdat %>%
    filter(year %in% years) -> newdat2
  
  print("predicting abundance")
  pred.abund <- predict(abund.mod, newdata= newdat2, return_tmb_object = T)
  print("predicting biomass")
  pred.bio <- predict(bio.mod, newdata= newdat2, return_tmb_object = T)
  
  write.csv(pred.abund$data, paste0("./BAIRDI/Output/", matsex, "_abundance_", stock, "_", period, "_", knots, "_spatialpreds.csv"))
  write.csv(pred.bio$data, paste0("./BAIRDI/Output/", matsex, "_biomass_", stock, "_", period, "_", knots, "_spatialpreds.csv"))
  
  # Get index
  print("getting abundance index")
  get_index(pred.abund, area = unique(newdat$Area_km2), bias_correct = TRUE) -> ind.abund
  print("getting biomass index")
  get_index(pred.bio, area = unique(newdat$Area_km2), bias_correct = TRUE) -> ind.bio
  
  
  write.csv(ind.abund, paste0("./BAIRDI/Output/", matsex, "_abundance_", stock, "_", period, "_", knots, "_index.csv"))
  write.csv(ind.bio, paste0("./BAIRDI/Output/", matsex, "_biomass_", stock, "_", period, "_", knots, "_index.csv"))
  
  
  return(list(pred.abund = pred.abund, pred.bio = pred.bio))
}

### PROCESS DATA -----------------------------------------------------------------

# Load and process response data
tan.cpue <- read.csv("./BAIRDI/Data/bairdi_cpue2.csv") %>%
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
tanE.obs <- rbind(read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/abund_EBS_TannerE.csv") %>% mutate(type = "abundance"),
                  read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/bio_EBS_TannerE.csv") %>% mutate(type = "biomass"))

tanW.obs <- rbind(read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/abund_EBS_TannerW.csv") %>% mutate(type = "abundance"),
                  read.csv("Y:/KOD_Survey/EBS Shelf/2024/Tech Memo/Outputs/bio_EBS_TannerW.csv") %>% mutate(type = "biomass"))

# Process observed abundance/biomass
tanE.obs %>%
  dplyr::select(Year, Small.male...113.mm., Large.male...113.mm., Immature.female, Mature.female, type) %>%
  rename(smmale = Small.male...113.mm., lgmale = Large.male...113.mm., 'Immature Female' = Immature.female, 'Mature Female' = Mature.female) %>%
  pivot_longer(., cols = c(2:5), names_to = "matsex", values_to = "value") %>%
  mutate(CI = as.numeric(gsub(",", "", sapply(str_extract_all(value, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =","))),
         value = as.numeric(gsub(",", "", gsub("\\([^()]*\\)", "", value))),
         stock = "TannerE",
         matsex = case_when((matsex %in% c("smmale", "lgmale")) ~ "Male",
                            TRUE ~ matsex)) %>%
  group_by(Year, type, matsex, stock) %>%
  reframe(value = sum(value),
          CI = sum(CI)) -> tanE.obs2


tanW.obs %>%
  dplyr::select(Year, Small.male...103.mm., Large.male...103.mm., Immature.female, Mature.female, type) %>%
  rename(smmale = Small.male...103.mm., lgmale = Large.male...103.mm., 'Immature Female' = Immature.female, 'Mature Female' = Mature.female) %>%
  pivot_longer(., cols = c(2:5), names_to = "matsex", values_to = "value") %>%
  mutate(CI = as.numeric(gsub(",", "", sapply(str_extract_all(value, "(?<=\\()[^)(]+(?=\\))"), paste0, collapse =","))),
         value = as.numeric(gsub(",", "", gsub("\\([^()]*\\)", "", value))),
         stock = "TannerW",
         matsex = case_when((matsex %in% c("smmale", "lgmale")) ~ "Male",
                            TRUE ~ matsex)) %>%
  group_by(Year, type, matsex, stock) %>%
  reframe(value = sum(value),
          CI = sum(CI)) -> tanW.obs2


tan.obs <- rbind(tanE.obs2, tanW.obs2)

