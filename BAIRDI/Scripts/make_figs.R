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
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "TannerW"
 
  ### Join timeseries
    # Spatial predictions
    Wmale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_spatialpreds.csv"))
    
    Wmale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_spatialpreds.csv"),
                            read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_spatialpreds.csv"))
    
    Wmale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wmale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wmale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_index.csv"),
                               read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_index.csv")) %>%
                           rename(abundance = est, Year = year) %>%
                           mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_index.csv"),
                             read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_index.csv")) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wmale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_index.csv")) %>%
                          rename(abundance = est, Year = year) %>%
                          mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_index.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_index.csv")) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)

  ### Plot timeseries
    # Spatial model predictions
    ggplot(Wmale.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wmale.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")

    # Indices 
    ggplot()+
      geom_ribbon(Wmale.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmale.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wmale.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmale.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "abundance"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "abundance"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("TannerW male estimated abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.male.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_ribbon(Wmale.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmale.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wmale.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmale.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "biomass"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "biomass"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW male estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.male.bio.index.png", width = 6, height = 4, units = "in")
  
  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerW"
  
    ### Join timeseries
    
    # Spatial predictions
    Wimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_index.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_index.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Wimfem.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wimfem.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    
    # Indices 
    ggplot()+
      geom_ribbon(Wimfem.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wimfem.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wimfem.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wimfem.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
      #geom_errorbar(abund, mapping = aes(x = Year, ymin = abundance-lwr, ymax = abundance+upr))+
      geom_errorbar(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "abundance"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "abundance"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("TannerW immature female estimated abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.imfem.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_ribbon(Wimfem.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wimfem.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wimfem.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wimfem.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "biomass"),
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "biomass"),
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW immature female estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.imfem.bio.index.png", width = 6, height = 4, units = "in")
    
  ## Mature females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerW"
      
    ### Join timeseries
    
    # Spatial predictions
    Wmatfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wmatfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_50_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_50_spatialpreds.csv"))
    
    Wmatfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_120_spatialpreds.csv"))
    
    Wmatfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Wmatfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmatfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_50_index.csv"),
                                read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Wmatfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Wmatfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Wmatfem.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Wmatfem.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    
    # Indices 
    ggplot()+
      geom_ribbon(Wmatfem.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmatfem.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wmatfem.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmatfem.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "abundance"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "abundance"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("TannerW mature female estimated abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanW.matfem.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_ribbon(Wmatfem.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Wmatfem.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Wmatfem.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Wmatfem.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "biomass"),
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == stock2, matsex == matsex2, type == "biomass"),
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerW mature female estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanW.matfem.bio.index.png", width = 6, height = 4, units = "in")
    
  
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

    
    ### Join timeseries
    
    # Spatial predictions
    Emale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_50_spatialpreds.csv"))
    
    Emale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_50_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_50_spatialpreds.csv"))
    
    Emale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_120_spatialpreds.csv"))
    
    Emale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_120_spatialpreds.csv"))
    
    # Indices
    Emale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_50_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Emale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_50_index.csv"),
                                read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_50_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    Emale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_120_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
    
    Emale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_120_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
    
    ### Plot timeseries
    
    # Spatial model predictions
    ggplot(Emale.spat.abund) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted abundance")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> abund_pred_plot
    
    ggplot(Emale.spat.bio) +
      geom_tile(aes(y = lat, x = lon, fill = est)) + 
      scale_fill_gradient2() + 
      labs(y = "Latitude",
           x = "Longitude",
           fill = "Log bairdi biomass (kg) per sq.km") +
      theme_gray() + 
      facet_wrap(~year)+
      ggtitle("Male TannerW predicted biomass")+
      theme(axis.title = element_text(size = 10),
            legend.position = "bottom") -> bio_pred_plot
    
    
    ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
    ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
    
    # Indices 
    ggplot()+
      geom_ribbon(Emale.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Emale.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Emale.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Emale.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == "tanE", matsex == "Male", type == "abundance"), 
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == "tanE", matsex == "Male", type == "abundance"), 
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Abundance (millions)")+
      ggtitle("TannerE male estimated abundance")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
    
    ggsave(plot = abund.index, "./BAIRDI/Figures/tanE.male.abund.index.png", width = 6, height = 4, units = "in")
    
    ggplot()+
      geom_ribbon(Emale.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
      geom_line(Emale.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
      geom_ribbon(Emale.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
      geom_line(Emale.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
      geom_errorbar(tan.obs %>% filter(stock == "tanE", matsex == "Male", type == "biomass"),
                    mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
      geom_point(tan.obs %>% filter(stock == "tanE", matsex == "Male", type == "biomass"),
                 mapping = aes(Year, value), color = "grey20", size = 1.5)+
      theme_bw()+
      ylab("Biomass (tons)")+
      ggtitle("TannerE male estimated biomass")+
      scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
      scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
    
    ggsave(plot = bio.index, "./BAIRDI/Figures/tanE.male.bio.index.png", width = 6, height = 4, units = "in")
    
    
  
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
      
      
      ### Join timeseries
      
      # Spatial predictions
      Eimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_50_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_50_spatialpreds.csv"))
      
      Eimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_50_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_50_spatialpreds.csv"))
      
      Eimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_120_spatialpreds.csv"),
                                   read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_120_spatialpreds.csv"))
      
      Eimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_120_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_120_spatialpreds.csv"))
      
      # Indices
      Eimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_50_index.csv"),
                                   read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_50_index.csv")) %>%
                              rename(abundance = est, Year = year) %>%
                              mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
      
      Eimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_50_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
      
      Eimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_120_index.csv"),
                                    read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_120_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
      
      Eimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_120_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
      
      ### Plot timeseries
      
      # Spatial model predictions
      ggplot(Eimfem.spat.abund) +
        geom_tile(aes(y = lat, x = lon, fill = est)) + 
        scale_fill_gradient2() + 
        labs(y = "Latitude",
             x = "Longitude",
             fill = "Log bairdi per sq.km") +
        theme_gray() + 
        facet_wrap(~year)+
        ggtitle("Immature female TannerE predicted abundance")+
        theme(axis.title = element_text(size = 10),
              legend.position = "bottom") -> abund_pred_plot
      
      ggplot(Eimfem.spat.bio) +
        geom_tile(aes(y = lat, x = lon, fill = est)) + 
        scale_fill_gradient2() + 
        labs(y = "Latitude",
             x = "Longitude",
             fill = "Log bairdi biomass (kg) per sq.km") +
        theme_gray() + 
        facet_wrap(~year)+
        ggtitle("Immature female TannerE predicted biomass")+
        theme(axis.title = element_text(size = 10),
              legend.position = "bottom") -> bio_pred_plot
      
      
      ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/ImmatureFemale_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
      ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/ImmatureFemale_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
      
      # Indices 
      ggplot()+
        geom_ribbon(Eimfem.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
        geom_line(Eimfem.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
        geom_ribbon(Eimfem.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
        geom_line(Eimfem.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
        geom_errorbar(tan.obs %>% filter(stock == "tanE", matsex == "Immature Female", type == "abundance"), 
                      mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
        geom_point(tan.obs %>% filter(stock == "tanE", matsex == "Immature Female", type == "abundance"), 
                   mapping = aes(Year, value), color = "grey20", size = 1.5)+
        theme_bw()+
        ylab("Abundance (millions)")+
        ggtitle("TannerE immature female estimated abundance")+
        scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
        scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
      
      ggsave(plot = abund.index, "./BAIRDI/Figures/tanE.imfem.abund.index.png", width = 6, height = 4, units = "in")
      
      ggplot()+
        geom_ribbon(Eimfem.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
        geom_line(Eimfem.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
        geom_ribbon(Eimfem.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
        geom_line(Eimfem.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
        geom_errorbar(tan.obs %>% filter(stock == "tanE", matsex == "Immature Female", type == "biomass"),
                      mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
        geom_point(tan.obs %>% filter(stock == "tanE", matsex == "Immature Female", type == "biomass"),
                   mapping = aes(Year, value), color = "grey20", size = 1.5)+
        theme_bw()+
        ylab("Biomass (tons)")+
        ggtitle("TannerE immature female estimated biomass")+
        scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
        scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
      
      ggsave(plot = bio.index, "./BAIRDI/Figures/tanE.imfem.bio.index.png", width = 6, height = 4, units = "in")
    
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
      
      ### Join timeseries
      
      # Spatial predictions
      Ematfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_50_spatialpreds.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_50_spatialpreds.csv"))
      
      Ematfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_50_spatialpreds.csv"))
      
      Ematfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_120_spatialpreds.csv"),
                                    read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_120_spatialpreds.csv"))
      
      Ematfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_120_spatialpreds.csv"))
      
      # Indices
      Ematfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_50_index.csv"),
                                    read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_50_index.csv")) %>%
        rename(abundance = est, Year = year) %>%
        mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
      
      Ematfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_50_index.csv")) %>%
        rename(biomass = est, Year = year) %>%
        mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
      
      Ematfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_120_index.csv"),
                                     read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_120_index.csv")) %>%
        rename(abundance = est, Year = year) %>%
        mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)
      
      Ematfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_120_index.csv")) %>%
        rename(biomass = est, Year = year) %>%
        mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)
      
      ### Plot timeseries
      
      # Spatial model predictions
      ggplot(Ematfem.spat.abund) +
        geom_tile(aes(y = lat, x = lon, fill = est)) + 
        scale_fill_gradient2() + 
        labs(y = "Latitude",
             x = "Longitude",
             fill = "Log bairdi per sq.km") +
        theme_gray() + 
        facet_wrap(~year)+
        ggtitle("Mature female TannerE predicted abundance")+
        theme(axis.title = element_text(size = 10),
              legend.position = "bottom") -> abund_pred_plot
      
      ggplot(Ematfem.spat.bio) +
        geom_tile(aes(y = lat, x = lon, fill = est)) + 
        scale_fill_gradient2() + 
        labs(y = "Latitude",
             x = "Longitude",
             fill = "Log bairdi biomass (kg) per sq.km") +
        theme_gray() + 
        facet_wrap(~year)+
        ggtitle("Mature female TannerE predicted biomass")+
        theme(axis.title = element_text(size = 10),
              legend.position = "bottom") -> bio_pred_plot
      
      
      ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Mature Female_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
      ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Mature Female_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
      
      # Indices 
      ggplot()+
        geom_ribbon(Ematfem.index.abund50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
        geom_line(Ematfem.index.abund50, mapping = aes(Year, abundance,  color = "salmon"), linewidth = 0.5)+
        geom_ribbon(Ematfem.index.abund120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
        geom_line(Ematfem.index.abund120, mapping = aes(Year, abundance, color = "turquoise"), linewidth = 0.5)+
        geom_errorbar(tan.obs %>% filter(stock == stock, matsex == matsex, type == "abundance"), 
                      mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
        geom_point(tan.obs %>% filter(stock == matsex, matsex == matsex, type == "abundance"), 
                   mapping = aes(Year, value), color = "grey20", size = 1.5)+
        theme_bw()+
        ylab("Abundance (millions)")+
        ggtitle("TannerE mature female estimated abundance")+
        scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
        scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> abund.index
      
      ggsave(plot = abund.index, "./BAIRDI/Figures/tanE.matfem.abund.index.png", width = 6, height = 4, units = "in")
      
      ggplot()+
        geom_ribbon(Ematfem.index.bio50, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "salmon"), alpha = 0.4) +
        geom_line(Ematfem.index.bio50, mapping = aes(Year, biomass,  color = "salmon"), linewidth = 0.5)+
        geom_ribbon(Ematfem.index.bio120, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = "turquoise"), alpha = 0.4) +
        geom_line(Ematfem.index.bio120, mapping = aes(Year, biomass, color = "turquoise"), linewidth = 0.5)+
        geom_errorbar(tan.obs %>% filter(stock == stock, matsex == matsex, type == "biomass"),
                      mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
        geom_point(tan.obs %>% filter(stock == stock, matsex == matsex, type == "biomass"),
                   mapping = aes(Year, value), color = "grey20", size = 1.5)+
        theme_bw()+
        ylab("Biomass (tons)")+
        ggtitle("TannerE mature female estimated biomass")+
        scale_color_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "")+
        scale_fill_manual(values = c("salmon", "turquoise"), labels = c("Index - 50", "Index - 120"), name = "") -> bio.index
      
      ggsave(plot = bio.index, "./BAIRDI/Figures/tanE.matfem.bio.index.png", width = 6, height = 4, units = "in")
      
    
### All ------------
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
  
  