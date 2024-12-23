source("./BAIRDI/load_libs_functions.R")

### LOAD FUNCTION -----------------------------------------------------------------------------------------

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
    
    # Predict and get index
    abund.mod1 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_DG_abundTMB.rda"))
    bio.mod1 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_DG_bioTMB.rda"))

    abund.mod2 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_DG_abundTMB.rda"))
    bio.mod2 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_DG_abundTMB.rda"))
    bio.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_DG_bioTMB.rda"))

    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90) ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50) ->  out

    
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  

    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1982_120_DG_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1982_120_DG_bioTMB.rda")

    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1982_90_DG_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1982_90_DG_bioTMB.rda")

    abund.mod3 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1982_50_DG_abundTMB.rda")
    bio.mod3 <- readRDS("./BAIRDI/Models/bairdi_Male_All_post-1982_50_DG_bioTMB.rda")

    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90) ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50) ->  out

    
  ## Immature Females -----  
  data <- tan.cpue2
  matsex <- "Immature Female"
  stock <- "All"
  
    ### Pre-1988
    years <- c(1975:1981)
    period <- "pre-1982"
    newdat <- pred_grid2
    

    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1982_120_DG_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1982_120_DG_bioTMB.rda")

    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1982_90_DG_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1982_90_DG_bioTMB.rda")

    abund.mod3 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1982_50_DG_abundTMB.rda")
    bio.mod3 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_pre-1982_50_DG_bioTMB.rda")

    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90) ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50) ->  out

    
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
    
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1982_120_DG_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1982_120_DG_bioTMB.rda")

    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1982_90_DG_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1982_90_DG_bioTMB.rda")

    abund.mod3 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1982_50_DG_abundTMB.rda")
    bio.mod3 <- readRDS("./BAIRDI/Models/bairdi_Immature Female_All_post-1982_50_DG_bioTMB.rda")

    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90) ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50) ->  out
    
  ## Mature Females -----  
  data <- tan.cpue2
  matsex <- "Mature Female"
  stock <- "All"
  
    ### Pre-1982
    years <- c(1975:1981)
    period <- "pre-1982"
    newdat <- pred_grid2
     

    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1982_120_DG_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1982_120_DG_bioTMB.rda")

    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1982_90_DG_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1982_90_DG_bioTMB.rda")

    abund.mod3 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1982_50_DG_abundTMB.rda")
    bio.mod3 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_pre-1982_50_DG_bioTMB.rda")

    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90) ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50) ->  out

    
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
      
    # Predict and get index
    abund.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1982_120_DG_abundTMB.rda")
    bio.mod1 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1982_120_DG_bioTMB.rda")

    abund.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1982_90_DG_abundTMB.rda")
    bio.mod2 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1982_90_DG_bioTMB.rda")

    abund.mod3 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1982_50_DG_abundTMB.rda")
    bio.mod3 <- readRDS("./BAIRDI/Models/bairdi_Mature Female_All_post-1982_50_DG_bioTMB.rda")

    predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120) ->  out
    predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90) ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50) ->  out

    
