source("./BAIRDI/Scripts/load_libs_functions.R")

dir <- "Y:/KOD_Research/Ryznar/Model-based indices/BAIRDI/"

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

# Delta gamma() X IID ----
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
    tan.cpue2 %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
    
    # Make mesh
    mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = 90, type = "kmeans")
    
    # Fit models
      abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                      spatial = "on",
                      spatiotemporal = "iid",
                      mesh = mesh2,
                      family = delta_gamma(type = "poisson-link"),
                      time = "year",
                      anisotropy = TRUE,
                      data = data2)
      
      saveRDS(abund, paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abund_DG_IID.rda"))
    
    # Predict and get index  
      newdat %>%
        filter(year %in% years) %>%
          mutate(year_fac = as.factor(year)) -> newdat2
        
      pred.abund <- predict(abund, newdata= newdat2, return_tmb_object = T)
      
      gc()
      get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
      
      write.csv(ind.abund, paste0(dir, "Output/Male_abundance_All_pre-1982_90_DG_IID_index.csv"))
      
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
    tan.cpue2 %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
    
    # Make mesh
    mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = 90, type = "kmeans")
  
      # Fit models
      abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                      spatial = "on",
                      spatiotemporal = "iid",
                      mesh = mesh2,
                      family = delta_gamma(type = "poisson-link"),
                      time = "year",
                      anisotropy = TRUE,
                      data = data2)
      
      saveRDS(abund, paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_DG_IID.rda"))
      
      abund <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_DG_IID.rda"))
      
    # Predict and get index 
      newdat %>%
        filter(year %in% years) %>%
        mutate(year_fac = as.factor(year)) -> newdat2
      
      pred.abund <- predict(abund, newdata= newdat2, return_tmb_object = T)
      
      gc()
      get_index_split(abund, area = unique(newdat2$Area_km2), bias_correct = TRUE, nsplit = 6) -> ind.abund
      
      write.csv(ind.abund, paste0(dir,"Output/Male_abundance_All_post-1982_90_DG_IID_index.csv"))

# Delta lognormal() x IID ----
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
    tan.cpue2 %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
      
      # Make mesh
      mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = 90, type = "kmeans")
    
      # Fit models
      abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                      spatial = "on",
                      spatiotemporal = "iid",
                      mesh = mesh2,
                      family = delta_lognormal(),
                      time = "year",
                      anisotropy = TRUE,
                      data = data2)
      
      saveRDS(abund, paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abund_DLN_IID.rda"))
      
      # Predict and get index  
      newdat %>%
        filter(year %in% years) %>%
        mutate(year_fac = as.factor(year)) -> newdat2
      
      pred.abund <- predict(abund, newdata= newdat2, return_tmb_object = T)
      
      gc()
      get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
      
      write.csv(ind.abund, paste0(dir, "Output/Male_abundance_All_pre-1982_90_DLN_IID_index.csv"))
    
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
    tan.cpue2 %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
    
      # Make mesh
      mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = 90, type = "kmeans")
      
      # Fit models
      abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                      spatial = "on",
                      spatiotemporal = "iid",
                      mesh = mesh2,
                      family = delta_lognormal(),
                      time = "year",
                      anisotropy = TRUE,
                      data = data2)
      
      saveRDS(abund, paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_DG_IID.rda"))
      
      abund <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_DG_IID.rda"))
      
      # Predict and get index  
      newdat %>%
        filter(year %in% years) %>%
        mutate(year_fac = as.factor(year)) -> newdat2
      
      pred.abund <- predict(abund, newdata= newdat2, return_tmb_object = T)
      
      gc()
      get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
      
      write.csv(ind.abund, paste0(dir, "Output/Male_abundance_All_post-1982_90_DLN_IID_index.csv"))
    
# Tweedie() x AR1 ----
  ### Pre-1982
  years <- c(1975:1981)
  period <- "pre-1982"
  newdat <- pred_grid2
  
    tan.cpue2 %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
    
      # Make mesh
      mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = 90, type = "kmeans")
      
      # Fit models
      abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                      spatial = "on",
                      spatiotemporal = "ar1",
                      mesh = mesh2,
                      family = tweedie(link = "log"),
                      time = "year",
                      anisotropy = TRUE,
                      data = data2)
      
      saveRDS(abund, "./BAIRDI/Models/bairdi_Male_All_pre-1982_90_abund_T_AR1.rda")
      
      # Predict and get index  
      newdat %>%
        filter(year %in% years) %>%
        mutate(year_fac = as.factor(year)) -> newdat2
      
      pred.abund <- predict(abund, newdata= newdat2, return_tmb_object = T)
      
      gc()
      get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
      
      write.csv(ind.abund, "./BAIRDI/Output/Male_abundance_All_pre-1982_90_T_AR1_index.csv")
  
  ### Post-1982
  years <- c(1982:2019, 2021:2024)
  period <- "post-1982"
  newdat <- pred_grid2
  
    tan.cpue2 %>%
      filter(mat.sex == matsex, year %in% years) %>%
      mutate(year_fac = as.factor(year)) -> data2
    
      # Make mesh
      mesh2 <- make_mesh(data2, c("lon","lat"), n_knots = 90, type = "kmeans") 
      
      # Fit models
      abund <- sdmTMB(cpue_km ~ 0 + as.factor(year), #the 0 is there so there is a factor predictor for each time slice
                      spatial = "on",
                      spatiotemporal = "ar1",
                      mesh = mesh2,
                      family = tweedie(link = "log"),
                      time = "year",
                      anisotropy = TRUE,
                      data = data2)
      
      saveRDS(abund, "./BAIRDI/Models/bairdi_Male_All_post-1982_90_abund_T_AR1.rda")

      # Predict and get index  
      newdat %>%
        filter(year %in% years) %>%
        mutate(year_fac = as.factor(year)) -> newdat2
      
      pred.abund <- predict(abund, newdata= newdat2, return_tmb_object = T)
      
      gc()
      get_index(pred.abund, area = unique(newdat2$Area_km2), bias_correct = TRUE) -> ind.abund
      
      write.csv(ind.abund, "./BAIRDI/Output/Male_abundance_All_post-1982_90_T_AR1_index.csv")

## PLOT ----
  # Load data
  DG_IID <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_90_DG_IID_index.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_90_DG_IID_index.csv")) %>%
                              rename(abundance = est, Year = year) %>%
                              mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                              mutate(model = "Delta-gamma, IID, 90kn")
  
  DLN_IID <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_90_DLN_IID_index.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_90_DLN_IID_index.csv")) %>%
                              rename(abundance = est, Year = year) %>%
                              mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                              mutate(model = "Delta-lognormal, IID, 90kn")
  
  T_AR1 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_90_T_AR1_index.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_90_T_AR1_index.csv")) %>%
                              rename(abundance = est, Year = year) %>%
                              mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                              mutate(model = "Tweedie, AR1, 90kn")
  
  T_IID <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_90_index.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_90_index.csv")) %>%
                              rename(abundance = est, Year = year) %>%
                              mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                              mutate(model = "Tweedie, IID, 90kn")
  
  
  # Join data
  rbind(DG_IID, DLN_IID, T_AR1, T_IID) -> All.abund.index
  
  All.abund.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, abundance = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.abund.index, dummy) -> All.abund.index2
  
  # Load obs data
  tan.obs %>%
    group_by(Year, type, matsex) %>%
    reframe(value = sum(value),
            CI = sum(CI)) -> tan.obs3
  
  # Plot indices
  ggplot()+
    geom_ribbon(All.abund.index2, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index2, mapping = aes(Year, abundance, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise", "violet", "lightgreen"), labels = c("50", "90", "120"), name = "Model")+
    scale_fill_manual(values = c("salmon", "turquoise", "violet", "lightgreen"), labels = c("50", "90", "120"), name = "Model")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> abund.ind.plot.EBS
