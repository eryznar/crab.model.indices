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
      filter(matsex == matsex, stock == stock, year %in% years) -> data2
  }else{
    data %>%
      filter(matsex == matsex, year %in% years) -> data2
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
  get_index(pred.abund, area = newdat$Area_km2, bias_correct = TRUE) -> ind.abund
  print("getting biomass index")
  get_index(pred.bio, area = newdat$Area_km2, bias_correct = TRUE) -> ind.bio
  
  
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
          cpue_kg_km = sum(cpue_kg_km))


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
               stock = "tanE",
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
             stock = "tanW",
             matsex = case_when((matsex %in% c("smmale", "lgmale")) ~ "Male",
                                TRUE ~ matsex)) %>%
        group_by(Year, type, matsex, stock) %>%
        reframe(value = sum(value),
                CI = sum(CI)) -> tanW.obs2
    
  
    tan.obs <- rbind(tanE.obs2, tanW.obs2)


### TANNER W ------------------------------------------------------------------
years <- c(1975:2019, 2021:2024)

# Filter pred_grid by stock, transform to UTM, replicate by number of years
W.grid <- pred_grid %>% filter(Lon < -166)

W.grid2 <- W.grid %>%
  st_as_sf(., coords = c("Lon", "Lat"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  dplyr::select(Area_km2, X, Y) %>%
  replicate_df(., "year", years) %>%
  mutate(X = X/1000, Y = Y/1000) %>%
  rename(lon = X, lat = Y)

## Males -----
data <- tan.cpue2
matsex <- "Male"
stock <- "TannerW"

  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- W.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out.Wmale1
    
    # Predict and get index
    abund.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_abundTMB.rda")
    bio.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots = 120) ->  ind.Wmale1

  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- W.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out.Wmale2
    
    # Predict models
    abund.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_abundTMB.rda")
    bio.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots = 120) -> ind.Wmale2

  ### Join timeseries
  
    # Spatial predictions
    Wmale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_spatialpreds.csv"))
    
    Wmale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_spatialpreds.csv"),
                            read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_spatialpreds.csv"))
    
    Wmale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wmale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wmale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_index.csv"),
                               read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_index.csv")) %>%
                         rename(abundance = est, Year = year) %>%
                         mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_index.csv"),
                             read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_index.csv")) %>%
                        rename(biomass = est, Year = year) %>%
                        mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wmale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_index.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)

  ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Wmale.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wmale.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")

    # Indices 
    ggplot()+
      geom_point(male.abund.obs %>% filter(stock == "tanW"), mapping = aes(Year, abundance), color = "grey20", size = 1.5)+
      geom_ribbon(Wmale.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmale.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 1)+
      geom_ribbon(Wmale.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmale.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 1)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(male.abund.obs %>% filter(stock == "tanW"), mapping = aes(x = Year, ymin = abundance - CI, ymax = abundance+CI), color = "grey20", width = 0)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("Tanner W male estimate abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.male.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_point(male.bio.obs %>% filter(stock == "tanW"), mapping = aes(Year, biomass), color = "grey20", size = 1.5)+
      geom_ribbon(Wmale.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmale.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 1)+
      geom_ribbon(Wmale.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmale.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 1)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(male.bio.obs %>% filter(stock == "tanW"), mapping = aes(x = Year, ymin = biomass - CI, ymax = biomass+CI), color = "grey20", width = 0)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW male estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.male.bio.index.png", width = 6, height = 4, units = "in")
    ### Join timeseries
  
    # Spatial predictions
    Wmale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_spatialpreds.csv"))
    
    Wmale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_spatialpreds.csv"),
                            read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_spatialpreds.csv"))
    
    Wmale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wmale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wmale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_index.csv"),
                               read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_index.csv")) %>%
                         rename(abundance = est, Year = year) %>%
                         mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_index.csv"),
                             read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_index.csv")) %>%
                        rename(biomass = est, Year = year) %>%
                        mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wmale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_index.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)

  ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Wmale.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wmale.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")

    # Indices 
    ggplot()+
      geom_point(male.abund.obs %>% filter(stock == "tanW"), mapping = aes(Year, abundance), color = "grey20", size = 1.5)+
      geom_ribbon(Wmale.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmale.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 1)+
      geom_ribbon(Wmale.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmale.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 1)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(male.abund.obs %>% filter(stock == "tanW"), mapping = aes(x = Year, ymin = abundance - CI, ymax = abundance+CI), color = "grey20", width = 0)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("Tanner W male estimate abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.male.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_point(male.bio.obs %>% filter(stock == "tanW"), mapping = aes(Year, biomass), color = "grey20", size = 1.5)+
      geom_ribbon(Wmale.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmale.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 1)+
      geom_ribbon(Wmale.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmale.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 1)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(male.bio.obs %>% filter(stock == "tanW"), mapping = aes(x = Year, ymin = biomass - CI, ymax = biomass+CI), color = "grey20", width = 0)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW male estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.male.bio.index.png", width = 6, height = 4, units = "in")
    

## Immature females -----
data <- tan.cpue2
matsex <- "Immature Female"
stock <- "TannerW"
    
    ### Pre-1988
    years <- c(1975:1987)
    period <- "pre-1988"
    newdat <- W.grid2 %>% filter(year %in% years)
    
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
    
    ### Post-1988
    years <- c(1988:2019, 2021:2024)
    period <- "post-1988"
    newdat <- W.grid2 %>% filter(year %in% years)
    
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
    

    ### Join timeseries
    
    # Spatial predictions
    Wimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_index.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_index.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Wimfem.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wimfem.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    
    # Indices 
    ggplot()+
      geom_ribbon(Wimfem.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wimfem.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wimfem.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wimfem.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "abundance"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "abundance"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("TannerW immature female estimated abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.imfem.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_ribbon(Wimfem.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wimfem.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wimfem.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wimfem.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "biomass"),
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "biomass"),
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW immature female estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.imfem.bio.index.png", width = 6, height = 4, units = "in")
    
    
## Mature females -----
data <- tan.cpue2
matsex <- "Mature Female"
stock <- "TannerW"
    
    ### Pre-1988
    years <- c(1975:1987)
    period <- "pre-1988"
    newdat <- W.grid2 %>% filter(year %in% years)
    
      # Fit models
      fit_models(data, matsex, stock, years, period, knots = 120) -> out
      fit_models(data, matsex, stock, years, period, knots = 50) -> out
      
      # Predict and get index
      abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_120_abundTMB.rda")
      bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_120_bioTMB.rda")
      
      abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_50_abundTMB.rda")
      bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_50_bioTMB.rda")
      
      predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
      predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
      
    ### Post-1988
    years <- c(1988:2019, 2021:2024)
    period <- "post-1988"
    newdat <- W.grid2 %>% filter(year %in% years)
    
      # Fit models
      fit_models(data, matsex, stock, years, period, knots = 120) -> out
      fit_models(data, matsex, stock, years, period, knots = 50) -> out
      
      # Predict and get index
      abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_120_abundTMB.rda")
      bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_120_bioTMB.rda")
      
      abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_50_abundTMB.rda")
      bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_50_bioTMB.rda")
      
      predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
      predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
      
    
    ### Join timeseries
    
    # Spatial predictions
    Wimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_index.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Wimfem.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wimfem.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    
    # Indices 
    ggplot()+
      geom_ribbon(Wimfem.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wimfem.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wimfem.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wimfem.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "abundance"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "abundance"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("TannerW immature female estimated abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.imfem.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_ribbon(Wimfem.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wimfem.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wimfem.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wimfem.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "biomass"),
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == "tanW", matsex == "Immature Female", type == "biomass"),
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW immature female estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.imfem.bio.index.png", width = 6, height = 4, units = "in")
    
    