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
 
  ### Join timeseries
    # Spatial predictions
    Wmale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_spatialpreds.csv")) %>% 
                          mutate(knots = 50, matsex = matsex2)
    
    Wmale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_spatialpreds.csv"),
                            read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_spatialpreds.csv")) %>% 
                         mutate(knots = 50, matsex = matsex2)
    
    Wmale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_spatialpreds.csv")) %>% 
                            mutate(knots = 120, matsex = matsex2)
    
    Wmale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_spatialpreds.csv")) %>% 
                            mutate(knots = 120, matsex = matsex2)
    
    # Indices
    Wmale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_index.csv"),
                               read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_index.csv")) %>%
                           rename(abundance = est, Year = year) %>%
                           mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6) %>% 
                          mutate(knots = 50, matsex = matsex2)
    
    Wmale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_index.csv"),
                             read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_index.csv")) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000) %>% 
                          mutate(knots = 50, matsex = matsex2)
    
    Wmale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_TannerW_post-1988_120_index.csv")) %>%
                          rename(abundance = est, Year = year) %>%
                          mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6) %>% 
                          mutate(knots = 120, matsex = matsex2)
    
    Wmale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerW_pre-1988_120_index.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_TannerW_post-1988_120_index.csv")) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000) %>% 
                          mutate(knots = 120, matsex = matsex2)

  ### Plot timeseries
    # Spatial model predictions
    # ggplot(Wmale.spat.abund) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted abundance")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> abund_pred_plot
    # 
    # ggplot(Wmale.spat.bio) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi biomass (kg) per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted biomass")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> bio_pred_plot
    # 
    # 
    # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")

  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerW"
  
    ### Join timeseries
    # Spatial predictions
    Wimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_spatialpreds.csv"))%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_spatialpreds.csv"),
                              read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_spatialpreds.csv"))%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_spatialpreds.csv"))%>% 
      mutate(knots = 120, matsex = matsex2)
    
    Wimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_spatialpreds.csv"))%>% 
      mutate(knots = 120, matsex = matsex2)
    
    # Indices
    Wimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_50_index.csv"),
                               read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_pre-1988_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
      mutate(knots = 120, matsex = matsex2)
    
    Wimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_pre-1988_120_index.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
      mutate(knots = 120, matsex = matsex2)
    
    ### Plot timeseries
    
    # # Spatial model predictions
    # ggplot(Wimfem.spat.abund) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted abundance")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> abund_pred_plot
    # 
    # ggplot(Wimfem.spat.bio) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi biomass (kg) per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted biomass")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> bio_pred_plot
    # 
    # 
    # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    # 
 
  ## Mature females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerW"
      
    ### Join timeseries
    
    # Spatial predictions
    Wmatfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_50_spatialpreds.csv")) %>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wmatfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_50_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_50_spatialpreds.csv"))%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wmatfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_120_spatialpreds.csv"))%>% 
      mutate(knots = 120, matsex = matsex2)
    
    Wmatfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_120_spatialpreds.csv"))%>% 
      mutate(knots = 120, matsex = matsex2)
    
    # Indices
    Wmatfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wmatfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_50_index.csv"),
                                read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Wmatfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_abundance_TannerW_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
      mutate(knots = 120, matsex = matsex2)
    
    Wmatfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_biomass_TannerW_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
      mutate(knots = 120, matsex = matsex2)
    
    ### Plot timeseries
    
    # # Spatial model predictions
    # ggplot(Wmatfem.spat.abund) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted abundance")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> abund_pred_plot
    # 
    # ggplot(Wmatfem.spat.bio) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi biomass (kg) per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted biomass")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> bio_pred_plot
    # 
    # 
    # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerW_predicted.png", height = 9, width = 8.5, units = "in")
    # 
    
  ## Join all data and plot -----
  rbind(Wmale.index.abund50, Wmale.index.abund120, Wimfem.index.abund50, 
        Wimfem.index.abund120, Wmatfem.index.abund50, Wmatfem.index.abund120) -> W.abund.index
    
  rbind(Wmale.index.bio50, Wmale.index.bio120, Wimfem.index.bio50, 
        Wimfem.index.bio120, Wmatfem.index.bio50, Wmatfem.index.bio120) -> W.bio.index
  
  rbind(Wmale.spat.abund50, Wmale.spat.abund120, Wimfem.spat.abund50, 
        Wimfem.spat.abund120, Wmatfem.spat.abund50, Wmatfem.spat.abund120) -> W.abund.spatial
  
  rbind(Wmale.spat.bio50, Wmale.spat.bio120, Wimfem.spat.bio50, 
        Wimfem.spat.bio120, Wmatfem.spat.bio50, Wmatfem.spat.bio120) -> W.bio.spatial
  
  
  # Plot indices
  ggplot()+
    geom_line(W.abund.index, mapping = aes(Year, abundance, color = as.factor(knots)))+
    geom_ribbon(W.abund.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4) +
    geom_point(tan.obs %>% filter(stock == stock2, type == "abundance"),
                          mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(tan.obs %>% filter(stock == stock2, type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner W estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal") -> abund.ind.plot.W
  
  ggsave(plot = abund.ind.plot.W, "./BAIRDI/Figures/TannerW.abundance.index.png", height= 6, width = 5, units = "in")
  
  ggplot()+
    geom_line(W.bio.index, mapping = aes(Year, biomass, color = as.factor(knots)))+
    geom_ribbon(W.bio.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4) +
    geom_point(tan.obs %>% filter(stock == stock2, type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(tan.obs %>% filter(stock == stock2, type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Biomass (tons)")+
    ggtitle("Tanner W estimated biomass") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "knots")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal") -> bio.ind.plot.W
  
  ggsave(plot = bio.ind.plot.W, "./BAIRDI/Figures/TannerW.biomass.index.png", height= 6, width = 5, units = "in")
  
    
    
### TANNER E ------------------------------------------------------------------
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "TannerE"
    
    ### Join timeseries
    
    # Spatial predictions
    Emale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_50_spatialpreds.csv")) %>% 
                            mutate(knots = 50, matsex = matsex2)
    
    Emale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_50_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_50_spatialpreds.csv")) %>% 
                            mutate(knots = 50, matsex = matsex2)
    
    Emale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_120_spatialpreds.csv"))%>% 
                          mutate(knots = 120, matsex = matsex2)
    
    Emale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_120_spatialpreds.csv"))%>% 
                            mutate(knots = 120, matsex = matsex2)
    
    # Indices
    Emale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_50_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                              mutate(knots = 50, matsex = matsex2)
    
    Emale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_50_index.csv"),
                                read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_50_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                            mutate(knots = 50, matsex = matsex2)
    
    Emale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_TannerE_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Male_abundance_TannerE_post-1988_120_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                            mutate(knots = 120, matsex = matsex2)
    
    Emale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_TannerE_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_biomass_TannerE_post-1988_120_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                          mutate(knots = 120, matsex = matsex2)
    
    ### Plot timeseries
    
    # # Spatial model predictions
    # ggplot(Emale.spat.abund) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted abundance")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> abund_pred_plot
    # 
    # ggplot(Emale.spat.bio) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi biomass (kg) per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Male TannerW predicted biomass")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> bio_pred_plot
    # 
    # 
    # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
    # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
  
  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerE"
  
      ### Join timeseries
      
      # Spatial predictions
      Eimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_50_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_50_spatialpreds.csv")) %>% 
        mutate(knots = 50, matsex = matsex2)
      
      Eimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_50_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_50_spatialpreds.csv")) %>% 
        mutate(knots = 50, matsex = matsex2)
      
      Eimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_120_spatialpreds.csv"),
                                   read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_120_spatialpreds.csv"))%>% 
        mutate(knots = 120, matsex = matsex2)
      
      Eimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_120_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_120_spatialpreds.csv"))%>% 
        mutate(knots = 120, matsex = matsex2)
      
      # Indices
      Eimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_50_index.csv"),
                                   read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_50_index.csv")) %>%
                              rename(abundance = est, Year = year) %>%
                              mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
        mutate(knots = 50, matsex = matsex2)
      
      Eimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_50_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
        mutate(knots = 50, matsex = matsex2)
      
      Eimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_pre-1988_120_index.csv"),
                                    read.csv("./BAIRDI/Output/Immature Female_abundance_TannerE_post-1988_120_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
        mutate(knots = 120, matsex = matsex2)
      
      Eimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_pre-1988_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_biomass_TannerE_post-1988_120_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
        mutate(knots = 120, matsex = matsex2)
      
      ### Plot timeseries
      
      # # Spatial model predictions
      # ggplot(Eimfem.spat.abund) +
      #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
      #   scale_fill_gradient2() + 
      #   labs(y = "Latitude",
      #        x = "Longitude",
      #        fill = "Log bairdi per sq.km") +
      #   theme_gray() + 
      #   facet_wrap(~year)+
      #   ggtitle("Immature female TannerE predicted abundance")+
      #   theme(axis.title = element_text(size = 10),
      #         legend.position = "bottom") -> abund_pred_plot
      # 
      # ggplot(Eimfem.spat.bio) +
      #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
      #   scale_fill_gradient2() + 
      #   labs(y = "Latitude",
      #        x = "Longitude",
      #        fill = "Log bairdi biomass (kg) per sq.km") +
      #   theme_gray() + 
      #   facet_wrap(~year)+
      #   ggtitle("Immature female TannerE predicted biomass")+
      #   theme(axis.title = element_text(size = 10),
      #         legend.position = "bottom") -> bio_pred_plot
      # 
      # 
      # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/ImmatureFemale_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
      # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/ImmatureFemale_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
      # 
  
  ## Mature females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerE"

    ### Join timeseries
    # Spatial predictions
    Ematfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_50_spatialpreds.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_50_spatialpreds.csv")) %>% 
      mutate(knots = 50, matsex = matsex2)
    
    Ematfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_50_spatialpreds.csv"),
                               read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_50_spatialpreds.csv"))%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Ematfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_120_spatialpreds.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_120_spatialpreds.csv"))%>% 
      mutate(knots = 120, matsex = matsex2)
    
    Ematfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_120_spatialpreds.csv"),
                                read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_120_spatialpreds.csv"))%>% 
      mutate(knots = 120, matsex = matsex2)
    
    # Indices
    Ematfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_50_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Ematfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_50_index.csv"),
                                read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_50_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
      mutate(knots = 50, matsex = matsex2)
    
    Ematfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_pre-1988_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_abundance_TannerE_post-1988_120_index.csv")) %>%
      rename(abundance = est, Year = year) %>%
      mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
      mutate(knots = 120, matsex = matsex2)
    
    Ematfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_pre-1988_120_index.csv"),
                                 read.csv("./BAIRDI/Output/Mature Female_biomass_TannerE_post-1988_120_index.csv")) %>%
      rename(biomass = est, Year = year) %>%
      mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
      mutate(knots = 120, matsex = matsex2)
    
    ### Plot timeseries
    
    # # Spatial model predictions
    # ggplot(Ematfem.spat.abund) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Mature female TannerE predicted abundance")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> abund_pred_plot
    # 
    # ggplot(Ematfem.spat.bio) +
    #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
    #   scale_fill_gradient2() + 
    #   labs(y = "Latitude",
    #        x = "Longitude",
    #        fill = "Log bairdi biomass (kg) per sq.km") +
    #   theme_gray() + 
    #   facet_wrap(~year)+
    #   ggtitle("Mature female TannerE predicted biomass")+
    #   theme(axis.title = element_text(size = 10),
    #         legend.position = "bottom") -> bio_pred_plot
    # 
    # 
    # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Mature Female_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
    # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Mature Female_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
    
  ## Join all data and plot -----
  rbind(Emale.index.abund50, Emale.index.abund120, Eimfem.index.abund50, 
        Eimfem.index.abund120, Ematfem.index.abund50, Ematfem.index.abund120) -> E.abund.index
  
  rbind(Emale.index.bio50, Emale.index.bio120, Eimfem.index.bio50, 
        Eimfem.index.bio120, Ematfem.index.bio50, Ematfem.index.bio120) -> E.bio.index
  
  rbind(Emale.spat.abund50, Emale.spat.abund120, Eimfem.spat.abund50, 
        Eimfem.spat.abund120, Ematfem.spat.abund50, Ematfem.spat.abund120) -> E.abund.spatial
  
  rbind(Emale.spat.bio50, Emale.spat.bio120, Eimfem.spat.bio50, 
        Eimfem.spat.bio120, Ematfem.spat.bio50, Ematfem.spat.bio120) -> E.bio.spatial
  
  
  # Plot indices
  ggplot()+
    geom_line(E.abund.index, mapping = aes(Year, abundance, color = as.factor(knots)))+
    geom_ribbon(E.abund.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4) +
    geom_point(tan.obs %>% filter(stock == stock2, type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(tan.obs %>% filter(stock == stock2, type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner E estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal") -> abund.ind.plot.E
  
  ggsave(plot = abund.ind.plot.E, "./BAIRDI/Figures/TannerE.abundance.index.png", height= 6, width = 5, units = "in")
  
  ggplot()+
    geom_line(E.bio.index, mapping = aes(Year, biomass, color = as.factor(knots)))+
    geom_ribbon(E.bio.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4) +
    geom_point(tan.obs %>% filter(stock == stock2, type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(tan.obs %>% filter(stock == stock2, type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Biomass (tons)")+
    ggtitle("Tanner E estimated biomass") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "knots")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal") -> bio.ind.plot.E
  
  ggsave(plot = bio.ind.plot.e, "./BAIRDI/Figures/TannerE.biomass.index.png", height= 6, width = 5, units = "in")
  
      
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
 
  # Spatial predictions
    # # 50
    # Allmale.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_50_spatialpreds.csv"),
    #                               read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_50_spatialpreds.csv")) %>% 
    #                           mutate(knots = 50, matsex = matsex2)
    # 
    # Allmale.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_All_pre-1982_50_spatialpreds.csv"),
    #                             read.csv("./BAIRDI/Output/Male_biomass_All_post-1982_50_spatialpreds.csv")) %>% 
    #                         mutate(knots = 50, matsex = matsex2)
    # 
    # # 90
    # Allmale.spat.abund90 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_90_spatialpreds.csv"),
    #                             read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_90_spatialpreds.csv")) %>% 
    #                         mutate(knots = 90, matsex = matsex2)
    # 
    # Allmale.spat.bio90 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_All_pre-1982_90_spatialpreds.csv"),
    #                           read.csv("./BAIRDI/Output/Male_biomass_All_post-1982_90_spatialpreds.csv")) %>% 
    #                        mutate(knots = 90, matsex = matsex2)
    # 
    # # 120
    # Allmale.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_120_spatialpreds.csv"),
    #                              read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_120_spatialpreds.csv"))%>% 
    #                       mutate(knots = 120, matsex = matsex2)
    # 
    # Allmale.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_All_pre-1982_120_spatialpreds.csv"),
    #                            read.csv("./BAIRDI/Output/Male_biomass_All_post-1982_120_spatialpreds.csv"))%>% 
    #                       mutate(knots = 120, matsex = matsex2)
  
  # Indices
    # 50
    Allmale.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_50_index.csv"),
                                   read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_50_index.csv")) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                            mutate(knots = 50, matsex = matsex2)
    
    Allmale.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_All_pre-1982_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_biomass_All_post-1982_50_index.csv")) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                            mutate(knots = 50, matsex = matsex2)
    
    # 80
    Allmale.index.abund90 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_90_index.csv"),
                                 read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_90_index.csv")) %>%
                          rename(abundance = est, Year = year) %>%
                          mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                          mutate(knots = 90, matsex = matsex2)
    
    Allmale.index.bio90 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_All_pre-1982_90_index.csv"),
                               read.csv("./BAIRDI/Output/Male_biomass_All_post-1982_90_index.csv")) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                          mutate(knots = 90, matsex = matsex2)
    
    # 120
    Allmale.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Male_abundance_All_pre-1982_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Male_abundance_All_post-1982_120_index.csv")) %>%
                          rename(abundance = est, Year = year) %>%
                          mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                          mutate(knots = 120, matsex = matsex2)
    
    Allmale.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Male_biomass_All_pre-1982_120_index.csv"),
                                read.csv("./BAIRDI/Output/Male_biomass_All_post-1982_120_index.csv")) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                          mutate(knots = 120, matsex = matsex2)
  
  ## Immature Females -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Spatial predictions
    # # 50 knots
    # Allimfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_All_pre-1982_50_spatialpreds.csv"),
    #                               read.csv("./BAIRDI/Output/Immature Female_abundance_All_post-1982_50_spatialpreds.csv")) %>% 
    #                              mutate(knots = 50, matsex = matsex2)
    # 
    # Allimfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_All_pre-1982_50_spatialpreds.csv"),
    #                             read.csv("./BAIRDI/Output/Immature Female_biomass_All_post-1982_50_spatialpreds.csv")) %>% 
    #                             mutate(knots = 50, matsex = matsex2)
    # 
    # # 90 knots
    # Allimfem.spat.abund90 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_All_pre-1982_90_spatialpreds.csv"),
    #                                read.csv("./BAIRDI/Output/Immature Female_abundance_All_post-1982_90_spatialpreds.csv")) %>% 
    #                               mutate(knots = 90, matsex = matsex2)
    # 
    # Allimfem.spat.bio90 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_All_pre-1982_90_spatialpreds.csv"),
    #                              read.csv("./BAIRDI/Output/Immature Female_biomass_All_post-1982_90_spatialpreds.csv")) %>% 
    #                             mutate(knots = 90, matsex = matsex2)
    # 
    # # 120 knots
    # Allimfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_All_pre-1982_120_spatialpreds.csv"),
    #                                read.csv("./BAIRDI/Output/Immature Female_abundance_All_post-1982_120_spatialpreds.csv"))%>% 
    #                               mutate(knots = 120, matsex = matsex2)
    # 
    # Allimfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_All_pre-1982_120_spatialpreds.csv"),
    #                              read.csv("./BAIRDI/Output/Immature Female_biomass_All_post-1982_120_spatialpreds.csv"))%>% 
    #                               mutate(knots = 120, matsex = matsex2)
    # 
  # Indices
    # 50
    Allimfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_All_pre-1982_50_index.csv"),
                                   read.csv("./BAIRDI/Output/Immature Female_abundance_All_post-1982_50_index.csv")) %>%
                                  rename(abundance = est, Year = year) %>%
                                  mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                                  mutate(knots = 50, matsex = matsex2)
    
    Allimfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_All_pre-1982_50_index.csv"),
                                 read.csv("./BAIRDI/Output/Immature Female_biomass_All_post-1982_50_index.csv")) %>%
                                  rename(biomass = est, Year = year) %>%
                                  mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                                  mutate(knots = 50, matsex = matsex2)
    
    # 90
    Allimfem.index.abund90 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_All_pre-1982_90_index.csv"),
                                    read.csv("./BAIRDI/Output/Immature Female_abundance_All_post-1982_90_index.csv")) %>%
                                    rename(abundance = est, Year = year) %>%
                                    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                                    mutate(knots = 90, matsex = matsex2)
    
    Allimfem.index.bio90 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_All_pre-1982_90_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_biomass_All_post-1982_90_index.csv")) %>%
                                  rename(biomass = est, Year = year) %>%
                                  mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                                  mutate(knots = 90, matsex = matsex2)
    
    # 120
    Allimfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_abundance_All_pre-1982_120_index.csv"),
                                    read.csv("./BAIRDI/Output/Immature Female_abundance_All_post-1982_120_index.csv")) %>%
                                    rename(abundance = est, Year = year) %>%
                                    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                                    mutate(knots = 120, matsex = matsex2)
    
    Allimfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Immature Female_biomass_All_pre-1982_120_index.csv"),
                                  read.csv("./BAIRDI/Output/Immature Female_biomass_All_post-1982_120_index.csv")) %>%
                                  rename(biomass = est, Year = year) %>%
                                  mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                                  mutate(knots = 120, matsex = matsex2)
    
  ## Mature Females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Spatial predictions
    # # 50
    # Allmatfem.spat.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_All_pre-1982_50_spatialpreds.csv"),
    #                                read.csv("./BAIRDI/Output/Mature Female_abundance_All_post-1982_50_spatialpreds.csv")) %>% 
    #                               mutate(knots = 50, matsex = matsex2)
    # 
    # Allmatfem.spat.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_All_pre-1982_50_spatialpreds.csv"),
    #                              read.csv("./BAIRDI/Output/Mature Female_biomass_All_post-1982_50_spatialpreds.csv")) %>% 
    #                               mutate(knots = 50, matsex = matsex2)
    # 
    # # 90
    # Allmatfem.spat.abund90 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_All_pre-1982_90_spatialpreds.csv"),
    #                                 read.csv("./BAIRDI/Output/Mature Female_abundance_All_post-1982_90_spatialpreds.csv")) %>% 
    #                                 mutate(knots = 90, matsex = matsex2)
    # 
    # Allmatfem.spat.bio90 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_All_pre-1982_90_spatialpreds.csv"),
    #                               read.csv("./BAIRDI/Output/Mature Female_biomass_All_post-1982_90_spatialpreds.csv")) %>% 
    #                               mutate(knots = 90, matsex = matsex2)
    # 
    # # 120
    # Allmatfem.spat.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_All_pre-1982_120_spatialpreds.csv"),
    #                                 read.csv("./BAIRDI/Output/Mature Female_abundance_All_post-1982_120_spatialpreds.csv"))%>% 
    #                                 mutate(knots = 120, matsex = matsex2)
    # 
    # Allmatfem.spat.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_All_pre-1982_120_spatialpreds.csv"),
    #                               read.csv("./BAIRDI/Output/Mature Female_biomass_All_post-1982_120_spatialpreds.csv"))%>% 
    #                               mutate(knots = 120, matsex = matsex2)
    
  # Indices
    # 50
    Allmatfem.index.abund50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_All_pre-1982_50_index.csv"),
                                    read.csv("./BAIRDI/Output/Mature Female_abundance_All_post-1982_50_index.csv")) %>%
                                    rename(abundance = est, Year = year) %>%
                                    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                                    mutate(knots = 50, matsex = matsex2)
    
    Allmatfem.index.bio50 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_All_pre-1982_50_index.csv"),
                                  read.csv("./BAIRDI/Output/Mature Female_biomass_All_post-1982_50_index.csv")) %>%
                                    rename(biomass = est, Year = year) %>%
                                    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                                    mutate(knots = 50, matsex = matsex2)
    
    # 90
    Allmatfem.index.abund90 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_All_pre-1982_90_index.csv"),
                                     read.csv("./BAIRDI/Output/Mature Female_abundance_All_post-1982_90_index.csv")) %>%
                                    rename(abundance = est, Year = year) %>%
                                    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                                    mutate(knots = 90, matsex = matsex2)
    
    Allmatfem.index.bio90 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_All_pre-1982_90_index.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_biomass_All_post-1982_90_index.csv")) %>%
                                    rename(biomass = est, Year = year) %>%
                                    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                                    mutate(knots = 90, matsex = matsex2)
    
    # 120
    Allmatfem.index.abund120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_abundance_All_pre-1982_120_index.csv"),
                                     read.csv("./BAIRDI/Output/Mature Female_abundance_All_post-1982_120_index.csv")) %>%
                                    rename(abundance = est, Year = year) %>%
                                    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                                    mutate(knots = 120, matsex = matsex2)
    
    Allmatfem.index.bio120 <- rbind(read.csv("./BAIRDI/Output/Mature Female_biomass_All_pre-1982_120_index.csv"),
                                   read.csv("./BAIRDI/Output/Mature Female_biomass_All_post-1982_120_index.csv")) %>%
                                    rename(biomass = est, Year = year) %>%
                                    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                                    mutate(knots = 120, matsex = matsex2)
  
  ### Plot timeseries
  
  # # Spatial model predictions
  # ggplot(Emale.spat.abund) +
  #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
  #   scale_fill_gradient2() + 
  #   labs(y = "Latitude",
  #        x = "Longitude",
  #        fill = "Log bairdi per sq.km") +
  #   theme_gray() + 
  #   facet_wrap(~year)+
  #   ggtitle("Male TannerW predicted abundance")+
  #   theme(axis.title = element_text(size = 10),
  #         legend.position = "bottom") -> abund_pred_plot
  # 
  # ggplot(Emale.spat.bio) +
  #   geom_tile(aes(y = lat, x = lon, fill = est)) + 
  #   scale_fill_gradient2() + 
  #   labs(y = "Latitude",
  #        x = "Longitude",
  #        fill = "Log bairdi biomass (kg) per sq.km") +
  #   theme_gray() + 
  #   facet_wrap(~year)+
  #   ggtitle("Male TannerW predicted biomass")+
  #   theme(axis.title = element_text(size = 10),
  #         legend.position = "bottom") -> bio_pred_plot
  # 
  # 
  # ggsave(plot = abund_pred_plot, "./BAIRDI/Figures/Male_abundance_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
  # ggsave(plot = bio_pred_plot, "./BAIRDI/Figures/Male_biomass_TannerE_predicted.png", height = 9, width = 8.5, units = "in")
  # 
 
  ## Join all data and plot -----
  rbind(Allmale.index.abund50, Allmale.index.abund90, Allmale.index.abund120, Allimfem.index.abund50,
        Allimfem.index.abund90, Allimfem.index.abund120, Allmatfem.index.abund50,
        Allmatfem.index.abund90, Allmatfem.index.abund120) -> All.abund.index
    
  All.abund.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, abundance = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.abund.index, dummy) -> All.abund.index
  
  rbind(Allmale.index.bio50, Allmale.index.bio90, Allmale.index.bio120, Allimfem.index.bio50, Allimfem.index.bio90,
        Allimfem.index.bio120, Allmatfem.index.bio50, Allmatfem.index.bio90, Allmatfem.index.bio120) -> All.bio.index
  
  All.bio.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, biomass = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.bio.index, dummy) -> All.bio.index
  
  
  # rbind(Allmale.spat.abund50, Allmale.spat.abund90, Allmale.spat.abund120, Allimfem.spat.abund50, Allimfem.spat.abund90,
  #       Allimfem.spat.abund120, Allmatfem.spat.abund50, Allmatfem.spat.abund90, Allmatfem.spat.abund120) -> All.abund.spatial
  # 
  # rbind(Allmale.spat.bio50, Allmale.spat.bio90, Allmale.spat.bio120, Allimfem.spat.bio50, Allimfem.spat.bio90,
  #       Allimfem.spat.bio120, Allmatfem.spat.bio50, Allmatfem.spat.bio90, Allmatfem.spat.bio120) -> All.bio.spatial
  
  tan.obs %>%
    group_by(Year, type, matsex) %>%
    reframe(value = sum(value),
            CI = sum(CI)) -> tan.obs3
  
  # Plot indices
  ggplot()+
    geom_ribbon(All.abund.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4)+
    geom_line(All.abund.index, mapping = aes(Year, abundance, color = as.factor(knots)))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> abund.ind.plot.EBS
  
  ggsave(plot = abund.ind.plot.EBS, "./BAIRDI/Figures/TannerEBS.abundance.index.png", height= 11, width = 8.5, units = "in")
  
  ggplot()+
    geom_ribbon(All.bio.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4) +
    geom_line(All.bio.index, mapping = aes(Year, biomass, color = as.factor(knots)))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Biomass (tons)")+
    ggtitle("EBS Tanner estimated biomass") +
    scale_color_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> bio.ind.plot.EBS
  
  ggsave(plot = bio.ind.plot.EBS, "./BAIRDI/Figures/TannerEBS.biomass.index.png", height= 11, width = 8.5, units = "in")
  