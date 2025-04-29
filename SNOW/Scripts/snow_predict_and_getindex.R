source("./SNOW/Scripts/load_libs_functions.R")

### LOAD FUNCTION -----------------------------------------------------------------------------------------

predict_and_getindex <- function(newdat, model, category, region, years, knots, dist){
  
  mod <- paste0(category, "-", region, "-", knots, "-", dist)
  
  newdat %>%
    filter(year %in% years) %>%
    mutate(year_fac = as.factor(year)) -> newdat2
  
  # Predict from model, get index
    print("predicting biomass")
    pred.bio <- predict(model, newdata= newdat2, return_tmb_object = T)
    
    gc()
    print("getting biomass index")
    get_index_split(model, newdata = newdat2, area = newdat2$area_km2, bias_correct = TRUE, nsplit = 3) -> ind.bio
   
  write.csv(ind.bio, paste0(dir, "Output/snow_", region, "_", category, "_", knots, "_", dist, "_biomassindex.csv"))
  
  
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
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_50_DG_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 50, dist = "DG") -> ebs.mf.DG50

## Mature female EBS 120 knots Delta gamma ----
data <- snow.matfem.cpue
category <- "Mature female"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_120_DG_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 120, dist = "DG") -> ebs.mf.DG120

## Mature female EBS 50 knots tweedie (120 knots doesn't fit)----
data <- snow.matfem.cpue
category <- "Mature female"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_50_TW_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 50, dist = "TW") -> ebs.mf.TW50

## Male95 EBS 50 knots Delta gamma ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_DG_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 50, dist = "DG") -> ebs.m.DG50

## Male95 EBS 120 knots Delta gamma ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_DG_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 120, dist = "DG") -> ebs.m.DG120

## Male95 EBS 50 knots tweedie ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_TW_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 50, dist = "TW") -> ebs.m.TW50

## Male95 EBS 120 knots tweedie ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_TW_bioTMB.rda"))

predict_and_getindex(ebs_grid2, model, category, region, years, knots = 120, dist = "TW") -> ebs.m.TW120