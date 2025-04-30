### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of biomass for EBS snow crab for males >95mm
# and mature females. Time range is 1980-present.

# Author: Emily Ryznar

# TO DOs:
# 1) MAKE SURE PACKAGE VERSIONS (sdmTMB, glmmTMB, Matrix, TMB...) ARE THE SAME BETWEEN DESKTOP AND VM!! 

### LOAD LIBRARIES/DATA -----------------------------------------------------------
source("./SNOW/Scripts/snow_load_libs_functions.R")

### LOAD FUNCTION -----------------------------------------------------------------
fit_models <- function(data, category, years, dist, knots, reg){
  
  # Filter data by params
  if(region != "All"){
    data %>%
      filter(category == category, year %in% years, region == reg) %>%
      mutate(year_fac = as.factor(year)) -> data2
  
  }else{
    data %>%
      filter(category == category, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
  }
  
  # # Set spatiotemporal estimator
  # if(region != "EBS"){
  #   sptmp <- "ar1"
  # }else{
  #   sptmp <- "iid"
  # }
  #  
  # Make mesh
  mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = knots, type = "kmeans")
  
  if(dist == "DG"){
    if(reg != "EBS"){
      # Fit models
      print("Fitting biomass model")
      bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "ar1",
                    mesh = mesh2,
                    family = delta_gamma(type = "poisson-link"),
                    time = "year",
                    extra_time = c(2020),
                    anisotropy = TRUE,
                    data = data2)
    } else{
      print("Fitting biomass model")
      bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "iid",
                    mesh = mesh2,
                    family = delta_gamma(type = "poisson-link"),
                    time = "year",
                    anisotropy = TRUE,
                    data = data2)
    }
   
    
  } else if(dist == "TW"){
    if(reg != "EBS"){
      # Fit models
      print("Fitting biomass model")
      bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "ar1",
                    mesh = mesh2,
                    family = tweedie(link = "log"),
                    time = "year",
                    extra_time = c(2020),
                    anisotropy = TRUE,
                    data = data2)
    } else{
      print("Fitting biomass model")
      bio <- sdmTMB(cpue_kg_km ~ 0 + year_fac, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    spatiotemporal = "iid",
                    mesh = mesh2,
                    family = tweedie(link = "log"),
                    time = "year",
                    anisotropy = TRUE,
                    data = data2)
    }
    
  } 
  
  saveRDS(bio, paste0(dir, "Models/snow_", reg, "_", category, "_", knots, "_", dist, "_bioTMB.rda"))
  
  return(list(bioTMB = bio, mesh = mesh2))
  
}

### SNOW EBS ------------------------------------------------------------------
years <- c(1980:2019, 2021:2024)

# Filter pred_grid by stock, transform to UTM, replicate by number of years
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", years) %>%
  rename(lon = X, lat = Y)

## Mature female EBS DG -----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "EBS"
dist <- "DG"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> ebs.mf.50
fit_models(data, category, years, dist, knots = 90, reg) -> ebs.mf.90
fit_models(data, category, years, dist, knots = 120, reg) -> ebs.mf.120

## Mature female EBS TW -----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "EBS"
dist <- "TW"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> ebs.mf.50
fit_models(data, category, years, dist, knots = 90, reg) -> ebs.mf.90 # doesn't fit!!!
fit_models(data, category, years, dist, knots = 120, reg) -> ebs.mf.120 # doesn't fit!!!

## Male EBS DG -----
data <- snow.male95.cpue 
category <- "Male95"
reg <- "EBS"
dist <- "DG"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> ebs.m.50
fit_models(data, category, years, dist, knots = 90, reg) -> ebs.m.90
fit_models(data, category, years, dist, knots = 120, reg) -> ebs.m.120

## Male EBS TW -----
data <- snow.male95.cpue 
category <- "Male95"
reg <- "EBS"
dist <- "TW"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> ebs.m.50
fit_models(data, category, years, dist, knots = 90, reg) -> ebs.m.90
fit_models(data, category, years, dist, knots = 120, reg) -> ebs.m.120


## Mature female EBS-NBS DG -----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "All"
dist <- "DG"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> all.mf.50
fit_models(data, category, years, dist, knots = 90, reg) -> all.mf.90
fit_models(data, category, years, dist, knots = 120, reg) -> all.mf.120

## Mature female EBS-NBS TW -----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "All"
dist <- "TW"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> all.mf.50
fit_models(data, category, years, dist, knots = 90, reg) -> all.mf.90
fit_models(data, category, years, dist, knots = 120, reg) -> all.mf.120

## Male EBS-NBS DG -----
data <- snow.male95.cpue 
category <- "Male95"
reg <- "All"
dist <- "DG"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> all.m.50
fit_models(data, category, years, dist, knots = 90, reg) -> all.m.90
fit_models(data, category, years, dist, knots = 120, reg) -> all.m.120

## Male EBS-NBS TW -----
data <- snow.male95.cpue 
category <- "Male95"
reg <- "All"
dist <- "TW"

# Fit models
fit_models(data, category, years, dist, knots = 50, reg) -> all.m.50 # doesn't fit!!!
fit_models(data, category, years, dist, knots = 90, reg) -> all.m.90
fit_models(data, category, years, dist, knots = 120, reg) -> all.m.120

