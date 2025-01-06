source("./BAIRDI/Scripts/load_libs_functions.R")

### LOAD FUNCTION -----------------------------------------------------------------------------------------

predict_and_getindex <- function(newdat, abund.mod, bio.mod, matsex, stock, years, period, knots, dist){
  
  mod <- paste0(knots, "-", dist)
  
  newdat %>%
    filter(year %in% years) %>%
    mutate(year_fac = as.factor(year)) -> newdat2
 
  if(period == "pre-1982"){
    print("predicting abundance")
    pred.abund <- predict(abund.mod, newdata= newdat2, return_tmb_object = T)
    print("predicting biomass")
    pred.bio <- predict(bio.mod, newdata= newdat2, return_tmb_object = T)

    gc()
    # Get index
    print("getting abundance index")
    get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
    gc()
    print("getting biomass index")
    get_index(pred.bio, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.bio
  } else{
    print("predicting abundance")
    pred.abund <- predict(abund.mod, newdata= newdat2, return_tmb_object = T)
    print("predicting biomass")
    pred.bio <- predict(bio.mod, newdata= newdat2, return_tmb_object = T)

    gc()
    print("getting abundance index")
    get_index_split(abund.mod, newdata = newdat2, area = unique(newdat2$Area_km2), bias_correct = TRUE, nsplit = 3) -> ind.abund

    gc()
    print("getting biomass index")
    get_index_split(bio.mod, newdata = newdat2, area = unique(newdat2$Area_km2), bias_correct = TRUE, nsplit = 3) -> ind.bio

  }



  # write.csv(pred.abund$data, paste0(dir, "Output/", matsex, "_abundance_", stock, "_", period, "_", mod, "_spatialpreds.csv"))
  # write.csv(pred.bio$data, paste0(dir, "Output/", matsex, "_biomass_", stock, "_", period, "_",mod, "_spatialpreds.csv"))


  write.csv(ind.abund, paste0(dir, "Output/", matsex, "_abundance_", stock, "_", period, "_", mod, "_index.csv"))
  write.csv(ind.bio, paste0(dir, "Output/", matsex, "_biomass_", stock, "_", period, "_",mod, "_index.csv"))


return(list(pred.abund = pred.abund$data, pred.bio = pred.bio$data))
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
    # abund.mod1 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_abund_DG_IID.rda"))
    # bio.mod1 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_DG_bioTMB.rda"))
    # 
    # abund.mod2 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abund_DG_IID.rda"))
    # bio.mod2 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abund_DG_IID.rda"))
    bio.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_DG_bioTMB.rda"))

    # predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120, "Delta_gamma") ->  out
    # predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90, "Delta_gamma") ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.male

    
    ### Post-1982
    years <- c(1982:2019, 2021:2024)
    period <- "post-1982"
    newdat <- pred_grid2
  

    # # Predict and get index
    # abund.mod1 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_120_abund_DG_IID.rda"))
    # bio.mod1 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_120_DG_bioTMB.rda"))
    # 
    # abund.mod2 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_90_abund_DG_IID.rda"))
    # bio.mod2 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_50_abund_DG_IID.rda"))
    bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_50_DG_bioTMB.rda"))

    # predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120, "Delta_gamma") ->  out
    # predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90, "Delta_gamma") ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.male

    ### Plot spatial predictions 
    abund <- rbind(pre.male$pred.abund %>% filter(lat<6841), post.male$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
             mutate(value = plogis(est1) * exp(est2))
    bio <- rbind(pre.male$pred.bio %>% filter(lat<6841), post.male$pred.bio) %>%
             mutate(value =plogis(est1) * exp(est2))
    
    
    ggplot(abund) +
      geom_tile(aes(y = lat, x = lon, fill = log(value))) +
      scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
      labs(y = "Latitude",
           x = "Longitude") +
      theme_bw() +
      scale_x_continuous(breaks = c(250, 750, 1250))+
      ggtitle("EBS predicted male abundance")+
      facet_wrap(~year)+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    ggsave("./BAIRDI/Figures/EBS_male_spatabund.png", width = 8.5, height = 9.5)
    
    ggplot(bio) +
      geom_tile(aes(y = lat, x = lon, fill = log(value))) +
      scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
      labs(y = "Latitude",
           x = "Longitude") +
      theme_bw() +
      scale_x_continuous(breaks = c(250, 750, 1250))+
      ggtitle("EBS predicted male biomass")+
      facet_wrap(~year)+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    ggsave("./BAIRDI/Figures/EBS_male_spatbio.png", width = 8.5, height = 9.5)
    
  ## Immature Females -----  
  data <- tan.cpue2
  matsex <- "Immature Female"
  stock <- "All"
  
    ### Pre-1988
    years <- c(1975:1981)
    period <- "pre-1982"
    newdat <- pred_grid2
    

    # Predict and get index
    # abund.mod1 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_120_DG_abundTMB.rda"))
    # bio.mod1 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_120_DG_bioTMB.rda"))
    # 
    # abund.mod2 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_90_DG_abundTMB.rda"))
    # bio.mod2 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_50_DG_abundTMB.rda"))
    bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_50_DG_bioTMB.rda"))

    # predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120, "Delta_gamma") ->  out
    # predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90,"Delta_gamma") ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.imfem

    
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
    
    # # Predict and get index
    # abund.mod1 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_120_DG_abundTMB.rda"))
    # bio.mod1 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_120_DG_bioTMB.rda"))
    # 
    # abund.mod2 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_90_DG_abundTMB.rda"))
    # bio.mod2 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_50_DG_abundTMB.rda"))
    bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_50_DG_bioTMB.rda"))

    # predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120, "Delta_gamma") ->  out
    # predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90, "Delta_gamma") ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.imfem
   
  ### Plot spatial predictions 
    abund <- rbind(pre.imfem$pred.abund %>% filter(lat<6841), post.imfem$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
      mutate(value = plogis(est1) * exp(est2))
    bio <- rbind(pre.imfem$pred.bio %>% filter(lat<6841), post.imfem$pred.bio) %>%
      mutate(value =plogis(est1) * exp(est2))
    
    
    ggplot(abund) +
      geom_tile(aes(y = lat, x = lon, fill = log(value))) +
      scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
      labs(y = "Latitude",
           x = "Longitude") +
      theme_bw() +
      scale_x_continuous(breaks = c(250, 750, 1250))+
      ggtitle("EBS predicted immature female abundance")+
      facet_wrap(~year)+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    ggsave("./BAIRDI/Figures/EBS_imfem_spatabund.png", width = 8.5, height = 9.5)
    
    ggplot(bio) +
      geom_tile(aes(y = lat, x = lon, fill = log(value))) +
      scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
      labs(y = "Latitude",
           x = "Longitude") +
      theme_bw() +
      scale_x_continuous(breaks = c(250, 750, 1250))+
      ggtitle("EBS predicted immature female biomass")+
      facet_wrap(~year)+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    ggsave("./BAIRDI/Figures/EBS_imfem_spatbio.png", width = 8.5, height = 9.5)  
  
     
  ## Mature Females -----  
  data <- tan.cpue2
  matsex <- "Mature Female"
  stock <- "All"
  
    ### Pre-1982
    years <- c(1975:1981)
    period <- "pre-1982"
    newdat <- pred_grid2
     

    # Predict and get index
    # abund.mod1 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_120_DG_abundTMB.rda"))
    # bio.mod1 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_120_DG_bioTMB.rda"))
    # 
    # abund.mod2 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_90_DG_abundTMB.rda"))
    # bio.mod2 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_50_DG_abundTMB.rda"))
    bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_50_DG_bioTMB.rda"))

    # predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120, "Delta_gamma") ->  out
    # predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90, "Delta_gamma") ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.matfem

    
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
      
    # # Predict and get index
    # abund.mod1 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_120_DG_abundTMB.rda"))
    # bio.mod1 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_120_DG_bioTMB.rda"))
    # 
    # abund.mod2 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_90_DG_abundTMB.rda"))
    # bio.mod2 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_90_DG_bioTMB.rda"))

    abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_50_DG_abundTMB.rda"))
    bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_50_DG_bioTMB.rda"))

    # predict_and_getindex(newdat, abund.mod1, bio.mod1, matsex, stock, years, period, knots = 120, "Delta_gamma") ->  out
    # predict_and_getindex(newdat, abund.mod2, bio.mod2, matsex, stock, years, period, knots = 90, "Delta_gamma") ->  out
    predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.matfem

  ### Plot spatial predictions 
    abund <- rbind(pre.matfem$pred.abund %>% filter(lat<6841), post.matfem$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
      mutate(value = plogis(est1) * exp(est2))
    bio <- rbind(pre.matfem$pred.bio %>% filter(lat<6841), post.matfem$pred.bio) %>%
      mutate(value =plogis(est1) * exp(est2))
    
    
    ggplot(abund) +
      #geom_sf(data = shoreline) +
      geom_tile(aes(y = lat, x = lon, fill = log(value))) +
      scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
      labs(y = "Latitude",
           x = "Longitude") +
      theme_bw() +
      scale_x_continuous(breaks = c(250, 750, 1250))+
      ggtitle("EBS predicted mature female abundance")+
      facet_wrap(~year)+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    ggsave("./BAIRDI/Figures/EBS_matfem_spatabund.png", width = 8.5, height = 9.5)
    
    ggplot(bio) +
      geom_tile(aes(y = lat, x = lon, fill = log(value))) +
      scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
      labs(y = "Latitude",
           x = "Longitude") +
      theme_bw() +
      scale_x_continuous(breaks = c(250, 750, 1250))+
      ggtitle("EBS predicted mature female biomass")+
      facet_wrap(~year)+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    ggsave("./BAIRDI/Figures/EBS_matfem_spatbio.png", width = 8.5, height = 9.5)  

### Tanner W --------------------------------------------------------------------
years <- c(1975:2019, 2021:2024)

pred_grid2 <- pg.W %>%
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
  stock <- "West"
  
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
  # Predict and get index
  abund.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abund_DG_IID.rda"))
  bio.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.male
  
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
  
  # # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_50_abund_DG_IID.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_50_DG_bioTMB.rda"))

  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.male
  
  ### Plot spatial predictions 
  abund <- rbind(pre.male$pred.abund %>% filter(lat<6841), post.male$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
    mutate(value = plogis(est1) * exp(est2))
  bio <- rbind(pre.male$pred.bio %>% filter(lat<6841), post.male$pred.bio) %>%
    mutate(value =plogis(est1) * exp(est2))
  
  
  ggplot(abund) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner West predicted male abundance")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerW_male_spatabund.png", width = 8.5, height = 9.5)
  
  ggplot(bio) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner West predicted male biomass")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerW_male_spatbio.png", width = 8.5, height = 9.5)
  
  ## Immature Females -----  
  data <- tan.cpue2
  matsex <- "Immature Female"
  stock <- "West"
  
  ### Pre-1988
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
  
  # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.imfem
  
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
  
  # # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.imfem
  
  ### Plot spatial predictions 
  abund <- rbind(pre.imfem$pred.abund %>% filter(lat<6841), post.imfem$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
    mutate(value = plogis(est1) * exp(est2))
  bio <- rbind(pre.imfem$pred.bio %>% filter(lat<6841), post.imfem$pred.bio) %>%
    mutate(value =plogis(est1) * exp(est2))
  
  
  ggplot(abund) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner West predicted immature female abundance")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerW_imfem_spatabund.png", width = 8.5, height = 9.5)
  
  ggplot(bio) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner West predicted immature female biomass")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerW_imfem_spatbio.png", width = 8.5, height = 9.5)  
  
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex <- "Mature Female"
  stock <- "West"
  
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
  
  # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.matfem
  
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
  
  # # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.matfem
  
  ### Plot spatial predictions 
  abund <- rbind(pre.matfem$pred.abund %>% filter(lat<6841), post.matfem$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
    mutate(value = plogis(est1) * exp(est2))
  bio <- rbind(pre.matfem$pred.bio %>% filter(lat<6841), post.matfem$pred.bio) %>%
    mutate(value =plogis(est1) * exp(est2))
  
  
  ggplot(abund) +
    #geom_sf(data = shoreline) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner West predicted mature female abundance")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerW_matfem_spatabund.png", width = 8.5, height = 9.5)
  
  ggplot(bio) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner West predicted mature female biomass")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerW_matfem_spatbio.png", width = 8.5, height = 9.5)  

### Tanner E --------------------------------------------------------------------
  years <- c(1975:2019, 2021:2024)
  
  pred_grid2 <- pg.E %>%
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
  stock <- "East"
  
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
  # Predict and get index
  abund.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abund_DG_IID.rda"))
  bio.mod3 <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.male
  
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
  
  # # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_50_abund_DG_IID.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Male_All_post-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.male
  
  ### Plot spatial predictions 
  abund <- rbind(pre.male$pred.abund %>% filter(lat<6841), post.male$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
    mutate(value = plogis(est1) * exp(est2))
  bio <- rbind(pre.male$pred.bio %>% filter(lat<6841), post.male$pred.bio) %>%
    mutate(value =plogis(est1) * exp(est2))
  
  
  ggplot(abund) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner East predicted male abundance")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerE_male_spatabund.png", width = 8.5, height = 9.5)
  
  ggplot(bio) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner East predicted male biomass")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerE_male_spatbio.png", width = 8.5, height = 9.5)
  
  ## Immature Females -----  
  data <- tan.cpue2
  matsex <- "Immature Female"
  stock <- "East"
  
  ### Pre-1988
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
  
  # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_pre-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.imfem
  
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
  
  # # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Immature Female_All_post-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.imfem
  
  ### Plot spatial predictions 
  abund <- rbind(pre.imfem$pred.abund %>% filter(lat<6841), post.imfem$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
    mutate(value = plogis(est1) * exp(est2))
  bio <- rbind(pre.imfem$pred.bio %>% filter(lat<6841), post.imfem$pred.bio) %>%
    mutate(value =plogis(est1) * exp(est2))
  
  
  ggplot(abund) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner East predicted immature female abundance")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerE_imfem_spatabund.png", width = 8.5, height = 9.5)
  
  ggplot(bio) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner East predicted immature female biomass")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerE_imfem_spatbio.png", width = 8.5, height = 9.5)  
  
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex <- "Mature Female"
  stock <- "East"
  
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
  
  # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_pre-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  pre.matfem
  
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
  
  # # Predict and get index
  abund.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_50_DG_abundTMB.rda"))
  bio.mod3 <- readRDS(paste0(dir,"Models/bairdi_Mature Female_All_post-1982_50_DG_bioTMB.rda"))
  
  predict_and_getindex(newdat, abund.mod3, bio.mod3, matsex, stock, years, period, knots = 50, "Delta_gamma") ->  post.matfem
  
  ### Plot spatial predictions 
  abund <- rbind(pre.matfem$pred.abund %>% filter(lat<6841), post.matfem$pred.abund) %>% # filtering lat so predictions do not extend beyond mesh extend for <1982 models
    mutate(value = plogis(est1) * exp(est2))
  bio <- rbind(pre.matfem$pred.bio %>% filter(lat<6841), post.matfem$pred.bio) %>%
    mutate(value =plogis(est1) * exp(est2))
  
  
  ggplot(abund) +
    #geom_sf(data = shoreline) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(num ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner East predicted mature female abundance")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerE_matfem_spatabund.png", width = 8.5, height = 9.5)
  
  ggplot(bio) +
    geom_tile(aes(y = lat, x = lon, fill = log(value))) +
    scale_fill_viridis_c(name = expression(paste("log(kg ", km^-2, ")")))+
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    scale_x_continuous(breaks = c(250, 750, 1250))+
    ggtitle("Tanner East predicted mature female biomass")+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ggsave("./BAIRDI/Figures/TannerE_matfem_spatbio.png", width = 8.5, height = 9.5)  
  