### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of biomass for EBS snow crab for males >95mm
# and mature females. Time range is 1980-present.

# Author: Emily Ryznar

# TO DOs:
# 1) MAKE SURE PACKAGE VERSIONS (sdmTMB, glmmTMB, Matrix, TMB...) ARE THE SAME BETWEEN DESKTOP AND VM!! 

### LOAD LIBS/PARAMBS
source("./SNOW/Scripts/snow_load_libs_functions.R")

### LOAD FUNCTION -----------------------------------------------------------------------------------------

predict_and_getindex <- function(newdat, model, category, reg, years, knots, dist){
  
  mod <- paste0(category, "-", reg, "-", knots, "-", dist)
  
  newdat %>%
    filter(year %in% years) %>%
    mutate(year_fac = as.factor(year)) -> newdat2
  
  # Predict from model, get index
    print("predicting biomass")
    pred.bio <- predict(model, newdata= newdat2, return_tmb_object = T)
    
    gc()
    print("getting biomass index")
    get_index_split(model, newdata = newdat2, area = newdat2$area_km2, bias_correct = TRUE, nsplit = 3) -> ind.bio
   
  write.csv(ind.bio, paste0(dir, "Output/Indices/snow_", reg, "_", category, "_", knots, "_", dist, "_biomassindex.csv"))
  
  
  return(list(pred.bio = pred.bio$data))
}

### EBS ----
years <- c(1980:2019, 2021:2024)

# Filter pred_grid by stock, transform to UTM, replicate by number of years
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", years) %>%
  rename(lon = X, lat = Y)

## Mature female EBS 50 knots Delta gamma ----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_120_DG_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "DG") -> ebs.mf.DG50
predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "DG") -> ebs.mf.DG90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "DG") -> ebs.mf.DG120

## Mature female EBS tweedie (90 and 120 knots doesn't fit)----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_50_TW_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "TW") -> ebs.mf.TW50

## Male95 EBS Delta gamma ----
data <- snow.male95.cpue
category <- "Male95"
reg <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_DG_bioTMB.rda"))


predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "DG") -> ebs.m.DG50
predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "DG") -> ebs.m.DG90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "DG") -> ebs.m.DG120

## Male95 EBS Tweedie ----
data <- snow.male95.cpue
category <- "Male95"
reg <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_TW_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_90_TW_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_TW_bioTMB.rda"))


predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "TW") -> ebs.m.TW50
predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "TW") -> ebs.m.TW90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "TW") -> ebs.m.TW120

## Mature female EBS-NBS 50 knots Delta gamma ----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "All"
model50 <- readRDS(paste0(dir, "Models/snow_All_Mature female_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_All_Mature female_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_All_Mature female_120_DG_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "DG") -> all.mf.DG50
predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "DG") -> all.mf.DG90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "DG") -> all.mf.DG120

## Mature female EBS-NBS tweedie ----
data <- snow.matfem.cpue
category <- "Mature female"
reg <- "All"
model50 <- readRDS(paste0(dir, "Models/snow_All_Mature female_50_TW_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_All_Mature female_90_TW_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_All_Mature female_120_TW_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "TW") -> all.mf.TW50
predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "TW") -> all.mf.TW90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "TW") -> all.mf.TW120


## Male95 EBS-NBS Delta gamma ----
data <- snow.male95.cpue
category <- "Male95"
reg <- "All"
model50 <- readRDS(paste0(dir, "Models/snow_All_Male95_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_All_Male95_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_All_Male95_120_DG_bioTMB.rda"))


predict_and_getindex(ebs_grid2, model50, category, reg, years, knots = 50, dist = "DG") -> all.m.DG50
predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "DG") -> all.m.DG90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "DG") -> all.m.DG120

## Male95 EBS-NBS Tweedie (50 knots doesn't fit) ----
data <- snow.male95.cpue
category <- "Male95"
reg <- "All"
model90 <- readRDS(paste0(dir, "Models/snow_All_Male95_90_TW_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_All_Male95_120_TW_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model90, category, reg, years, knots = 90, dist = "TW") -> all.m.TW90
predict_and_getindex(ebs_grid2, model120, category, reg, years, knots = 120, dist = "TW") -> all.m.TW120
