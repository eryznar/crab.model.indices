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
pred_grid <- readRDS(here::here(paste0(dir, "Data/EBS_bairdi_grid_5km_No_Land.rds")) )

### LOAD FUNCTIONS -------------------------------------------------------------
fit_models <- function(data, matsex, stock, years, period, knots){
  
  # Filter data by params
  if(stock != "All"){
    data %>%
      filter(mat.sex == matsex, stck == stock, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
  }else{
    data %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
  }
  
  # Make mesh
  mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = knots, type = "kmeans")
  
  
  # Fit models
  print("Fitting abundance model")
  abund <- sdmTMB(cpue_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                  spatial = "on",
                  spatiotemporal = "iid",
                  mesh = mesh2,
                  family = tweedie(link = "log"),
                  #family = delta_lognormal(),
                  time = "year",
                  anisotropy = TRUE,
                  data = data2)
  
  saveRDS(abund, paste0(dir, "Models/bairdi_", matsex, "_", stock, "_", period, "_", knots, "_abundTMB.rda"))
  
  print("Fitting biomass model")
  bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                spatial = "on",
                spatiotemporal = "iid",
                mesh = mesh2,
                family = tweedie(link = "log"),
                #family = delta_lognormal(),
                time = "year",
                anisotropy = TRUE,
                data = data2)
  
  saveRDS(bio, paste0(dir, "Models/bairdi_", matsex, "_", stock, "_", period, "_", knots, "_bioTMB.rda"))
  
  return(list(abundTMB = abund, bioTMB = bio, mesh = mesh2))
}

evaluate_diagnostics <- function(data, pre.model, post.model, stock2, type, knots, family, method, matsex2){
  
  mod <- paste0(type, "-", knots, "-", family)
  
  # Run sanity check
  print("Pre-1982 model")
  sanity_check_pre <- sanity(pre.model) 
  print("Post-1982 model")
  sanity_check_post <- sanity(post.model) 
  
  # Calculate Dharma residuals
  resid1 <- simulate(pre.model, nsim = 300, type= "mle-mvn")|>
                dharma_residuals(pre.model, plot = FALSE)
  
  ggplot()+
    theme_bw()+
    geom_point(resid1, mapping = aes(expected, observed), size = 2, fill = "black")+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle("<1982")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) -> r1.plot
    # 
  # png(filename = paste0("./BAIRDI/Figures/DHARMa_pre1988_", type, "_", stock2, "_", matsex2, "_", knots, ".png"), width=7, height=5, units="in", res=600)
  # 
  # plot(resid1, title= paste0("DHARMa residuals (Pre-1988, ", type, ",", stock2, " ", matsex2, ", knots=", knots, ")"))
  # 
  # dev.off()
  
  resid2 <- simulate(post.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(post.model, plot= FALSE)
  
  ggplot()+
    theme_bw()+
    geom_point(resid2, mapping = aes(expected, observed), size = 2, fill = "black")+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle("≥1982")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) -> r2.plot
  
  rbind(resid1 %>% mutate(period = "<1982"), resid2 %>% mutate(period = "≥1982")) %>%
    mutate(matsex = matsex2) -> all.resids
  
 # Test residuals
  resid1 <- simulate(pre.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(pre.model, return_DHARMa = TRUE)
  resid2 <- simulate(post.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(post.model, return_DHARMa = TRUE)
  
  round(DHARMa::testQuantiles(resid1, plot = FALSE)$p.value, 2) -> qq1
  round(DHARMa::testQuantiles(resid2, plot = FALSE)$p.value, 2) -> qq2
  
  round(DHARMa::testDispersion(resid1, plot = FALSE)$p.value, 2) -> dd1
  round(DHARMa::testDispersion(resid2, plot = FALSE)$p.value,2) -> dd2
  
  round(DHARMa::testOutliers(resid1, plot = FALSE)$p.value, 2) -> oo1
  round(DHARMa::testOutliers(resid2, plot = FALSE)$p.value, 2) -> oo2
  
  round(DHARMa::testZeroInflation(resid1, plot = FALSE)$p.value, 2) -> zz1
  round(DHARMa::testZeroInflation(resid2, plot = FALSE)$p.value, 2) -> zz2
  
  
  if(stock2 != "All"){
    data %>% 
      filter(mat.sex == matsex2, stck == stock2) %>%
      mutate(resids = c(resid1$scaledResiduals, resid2$scaledResiduals)) -> data2
  } else{
    data %>% 
      filter(mat.sex == matsex2) %>%
      mutate(resids = c(resid1$scaledResiduals, resid2$scaledResiduals)) -> data2
  }
  
  
  # visualize residuals across the EBS
  ggplot(data2) + 
    #geom_sf(data = shoreline) +
    geom_point(aes(y = lat, x = lon, color = resids), size = 1) +
    scale_color_gradient2(midpoint = 0) + 
    labs(y = "Latitude",
         x = "Longitude") +
    theme_gray() + 
    ggtitle(paste0(stock2," ", matsex2," ", type, " residuals (knots=", knots, ")"))+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom") -> res_plot
  
  ggsave(plot = res_plot, paste0("./BAIRDI/Figures/DHARMa", stock2, "_", matsex2, "_", mod, "_SPATIAL.png"), height = 9, width = 8.5, units = "in")
  
  # Calculate log-likelihood
  clust <- sample(seq_len(10), size = nrow(pre.model$data), replace = TRUE)
  
   pre.ll <- sdmTMB_cv(
      data = pre.model$data, 
      formula = cpue_km ~ 0 + year_fac, 
      spatial = "on",
      time = "year", 
      mesh = pre.model$spde, 
      spatiotemporal = "iid",
      silent = FALSE,
      anisotropy = TRUE,
      family = tweedie(link = "log"),
      fold_ids = clust
    )
   
   clust <- sample(seq_len(10), size = nrow(post.model$data), replace = TRUE)
   
   post.ll <- sdmTMB_cv(
     data = post.model$data, 
     formula = cpue_km ~ 0 + year_fac, 
     spatial = "on",
     time = "year", 
     mesh = post.model$spde, 
     spatiotemporal = "iid",
     silent = FALSE,
     anisotropy = TRUE,
     family = tweedie(link = "log"),
     fold_ids = clust
   )
   
   eval.df <- data.frame(matsex = rep(matsex2,2), 
                         type = rep(type,2),
                         knots = rep(knots,2), 
                         family = rep(family, 2), 
                         method = rep(method,2), 
                         period = c("<1982", "≥1982"),
                         loglik = c(pre.ll$sum_loglik, post.ll$sum_loglik),
                         quantiles = c(qq1, qq2),
                         dispersion = c(dd1, dd1),
                         outliers = c(oo1, oo2),
                         zeroinf = c(zz1, zz2))
              
  
  return(list(sanity_check_pre, sanity_check_post, all.resids, eval.df))
}

predict_and_getindex <- function(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots){
  newdat %>%
    filter(year %in% years) %>%
    mutate(year_fac = as.factor(year)) -> newdat2
  
  print("predicting abundance")
  pred.abund <- predict(abund.mod, newdata= newdat2, return_tmb_object = T)
  print("predicting biomass")
  pred.bio <- predict(bio.mod, newdata= newdat2, return_tmb_object = T)
  
  write.csv(pred.abund$data, paste0("./BAIRDI/Output/", matsex, "_abundance_", stock, "_", period, "_", knots, "_spatialpreds.csv"))
  write.csv(pred.bio$data, paste0("./BAIRDI/Output/", matsex, "_biomass_", stock, "_", period, "_", knots, "_spatialpreds.csv"))
  
  gc()
  # Get index
  print("getting abundance index")
  get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
  gc()
  print("getting biomass index")
  get_index(pred.bio, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.bio
  
  
  write.csv(ind.abund, paste0("./BAIRDI/Output/", matsex, "_abundance_", stock, "_", period, "_", knots, "_index.csv"))
  write.csv(ind.bio, paste0("./BAIRDI/Output/", matsex, "_biomass_", stock, "_", period, "_", knots, "_index.csv"))
  
  
  return(list(pred.abund = pred.abund, pred.bio = pred.bio))
}

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
