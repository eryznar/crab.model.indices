### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Look at residuals
# 2) Troubleshoot immature female abundance
# 3) Get abund/bio obs back to 1975

### LOAD LIBRARIES/FUNCTIONS/DATA --------------------------------------------------------
source("./BAIRDI/Scripts/load_libs_functions.R")

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
    fit_models(data, matsex, stock, years, period, knots = 50) -> out.Wmale2
    fit_models(data, matsex, stock, years, period, knots = 120) -> out.Wmale1
    
    # Predict and get index
    abund.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_50_abundTMB.rda")
    bio.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_50_bioTMB.rda")
    
    abund.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_abundTMB.rda")
    bio.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots = 120) ->  ind.Wmale1
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- W.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 50) -> out.Wmale2
    fit_models(data, matsex, stock, years, period, knots = 120) -> out.Wmale2
    
    # Predict models
    abund.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_50_abundTMB.rda")
    bio.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_50_bioTMB.rda")
    
    abund.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_120_abundTMB.rda")
    bio.mod <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_120_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots = 120) -> ind.Wmale2

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
  
  
  
### TANNER E ------------------------------------------------------------------
years <- c(1975:2019, 2021:2024)
  
# Filter pred_grid by stock, transform to UTM, replicate by number of years
E.grid <- pred_grid %>% filter(Lon > -166)

E.grid2 <- E.grid %>%
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
  stock <- "TannerE"
  
  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- E.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- E.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out

  ## Immature females -----
  data <- tan.cpue2
  matsex <- "Immature Female"
  stock <- "TannerE"
  
  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- E.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- E.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ## Mature females -----  
  data <- tan.cpue2
  matsex <- "Mature Female"
  stock <- "TannerE"
  
  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- E.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- E.grid2 %>% filter(year %in% years)
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
### ALL ---------------------------------------------------------------------------------------------------
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
  matsex <- "Male"
  stock <- "All"
  
  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_All_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_All_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_All_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_All_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ## Immature Females -----  
  data <- tan.cpue2
  matsex <- "Immature Female"
  stock <- "All"
  
  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex <- "Mature Female"
  stock <- "All"
  
  ### Pre-1988
  years <- c(1975:1987)
  period <- "pre-1988"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out
  
  ### Post-1988
  years <- c(1988:2019, 2021:2024)
  period <- "post-1988"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, knots = 120) -> out
    fit_models(data, matsex, stock, years, period, knots = 50) -> out
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1988_120_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1988_120_bioTMB.rda")
    
    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1988_50_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1988_50_bioTMB.rda")
    
    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 50) ->  out  
  
