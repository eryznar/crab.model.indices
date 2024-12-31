### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Look at residuals
# 2) Troubleshoot immature female abundance
# 3) Get abund/bio obs back to 1975

### LOAD LIBRARIES/DATA -----------------------------------------------------------
source("./BAIRDI/Scripts/load_libs_functions.R")

### LOAD FUNCTION -----------------------------------------------------------------
fit_models <- function(data, matsex, stock, years, period, dist, knots){
  
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
  
  if(dist == "DG"){
    # Fit models
    print("Fitting abundance model")
    abund <- sdmTMB(cpue_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "iid",
                    mesh = mesh2,
                    family = delta_gamma(type = "poisson-link"),
                    time = "year",
                    anisotropy = TRUE,
                    data = data2)
    
    print("Fitting biomass model")
    bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                  spatial = "on",
                  spatiotemporal = "iid",
                  mesh = mesh2,
                  family = delta_gamma(type = "poisson-link"),
                  time = "year",
                  anisotropy = TRUE,
                  data = data2)
    
  } else if(dist == "TW"){
    # Fit models
    print("Fitting abundance model")
    abund <- sdmTMB(cpue_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "iid",
                    mesh = mesh2,
                    family = tweedie(link = "log"),
                    time = "year",
                    anisotropy = TRUE,
                    data = data2)
    
    print("Fitting biomass model")
    bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                  spatial = "on",
                  spatiotemporal = "iid",
                  mesh = mesh2,
                  family = tweedie(link = "log"),
                  time = "year",
                  anisotropy = TRUE,
                  data = data2)
  } else{
    # Fit models
    print("Fitting abundance model")
    abund <- sdmTMB(cpue_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "iid",
                    mesh = mesh2,
                    family = delta_lognormal(),
                    time = "year",
                    anisotropy = TRUE,
                    data = data2)
    
    print("Fitting biomass model")
    bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                  spatial = "on",
                  spatiotemporal = "iid",
                  mesh = mesh2,
                  family = delta_lognormal(),
                  time = "year",
                  anisotropy = TRUE,
                  data = data2)
    
  }
  
  
  
  
  saveRDS(abund, paste0(dir, "Models/bairdi_", matsex, "_", stock, "_", period, "_", knots, "_", dist, "_abundTMB.rda"))
  
  
  saveRDS(bio, paste0(dir, "Models/bairdi_", matsex, "_", stock, "_", period, "_", knots, "_", dist, "_bioTMB.rda"))
  
  return(list(abundTMB = abund, bioTMB = bio, mesh = mesh2))
  
}

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
    fit_models(data, matsex, stock, years, period, knots = 40) -> out.Wmale2
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
  
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, dist = "DG", knots = 120) -> out
    fit_models(data, matsex, stock, years, period, dist = "DG", knots = 90) -> out
    fit_models(data, matsex, stock, years, period, dist = "DG", knots = 50) -> out

  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
    # Fit models
    fit_models(data, matsex, stock, years, period, dist = "DG", knots = 120) -> out
    fit_models(data, matsex, stock, years, period, dist = "DG", knots = 90) -> out
    fit_models(data, matsex, stock, years, period, dist = "DG", knots = 50) -> out
    
   
# 1) Fit male biomass DG models across all knots, 2) Fit female bio/abund DG models across all knots
# 3) Eval diagnostics for female bio/abund @50 knots 4) Predict/get index for all DG models at 50 knots