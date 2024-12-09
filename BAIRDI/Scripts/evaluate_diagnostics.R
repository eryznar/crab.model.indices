### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Look at residuals
# 2) Add in scripts to load new survey data (CPUE, BIO/ABUND) and process each year (CPUE script is in TECHMEMONEW)

### LOAD LIBRARIES/FUNCTIONS/DATA --------------------------------------------------------
source("./BAIRDI/Scripts/load_libs_functions.R")

### TANNER W ------------------------------------------------------------------
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "TannerW"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
 
  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerW"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  ## Mature females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerW"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
 
### TANNER E ------------------------------------------------------------------
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "TannerE"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerE"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  ## Mature females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerE"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
### EBS-wide ------------
years <- c(1975:2019, 2021:2024)

pred_grid2 <- pred_grid %>%
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
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_All_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_All_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out

  ## Immature Females -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  