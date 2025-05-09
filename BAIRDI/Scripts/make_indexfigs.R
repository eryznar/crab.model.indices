### PURPOSE ----------------------------------------------------------------------
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Look at residuals
# 2) Add in scripts to load new survey data (CPUE, BIO/ABUND) and process each year (CPUE script is in TECHMEMONEW)

### LOAD LIBRARIES/FUNCTIONS/DATA --------------------------------------------------------
source("./BAIRDI/Scripts/load_libs_functions.R")
### Load VAST indices ---------------------------------
dir2 <- paste0(dir, "Data/VAST/")

# Tweedie
imfem.TW <- rbind(read.csv(paste0(dir2, "VAST_imfem_TW_abund50.csv")) %>%
                    mutate(knots = 50, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_TW_abund120.csv")) %>%
                    mutate(knots = 120, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_TW_bio50.csv")) %>%
                    mutate(knots = 50, type = "biomass"),
                  read.csv(paste0(dir2, "VAST_imfem_TW_bio120.csv")) %>%
                    mutate(knots = 120, type = "biomass")) %>%
  mutate(matsex = "Immature Female", family = "Tweedie")

matfem.TW <- rbind(read.csv(paste0(dir2, "VAST_matfem_TW_abund50.csv")) %>%
                     mutate(knots = 50, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_matfem_TW_abund120.csv")) %>%
                     mutate(knots = 120, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_matfem_TW_bio50.csv")) %>%
                     mutate(knots = 50, type = "biomass"),
                   read.csv(paste0(dir2, "VAST_matfem_TW_bio120.csv")) %>%
                     mutate(knots = 120, type = "biomass")) %>%
  mutate(matsex = "Mature Female", family = "Tweedie")


male.TW <- rbind(read.csv(paste0(dir2, "VAST_male_TW_abund50.csv")) %>%
                   mutate(knots = 50, type = "abundance"),
                 read.csv(paste0(dir2, "VAST_male_TW_abund120.csv")) %>%
                   mutate(knots = 120, type = "abundance"),
                 read.csv(paste0(dir2, "VAST_male_TW_bio50.csv")) %>%
                   mutate(knots = 50, type = "biomass"),
                 read.csv(paste0(dir2, "VAST_male_TW_bio120.csv")) %>%
                   mutate(knots = 120, type = "biomass")) %>%
  mutate(matsex = "Male", family = "Tweedie")

# Delta gamma
imfem.DG <- rbind(read.csv(paste0(dir2, "VAST_imfem_DG_abund50.csv")) %>%
                    mutate(knots = 50, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_abund90.csv")) %>%
                    mutate(knots = 90, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_abund120.csv")) %>%
                    mutate(knots = 120, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_abund750.csv")) %>%
                    mutate(knots = 750, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_bio50.csv")) %>%
                    mutate(knots = 50, type = "biomass"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_bio90.csv")) %>%
                    mutate(knots = 90, type = "biomass"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_bio120.csv")) %>%
                    mutate(knots = 120, type = "biomass"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_bio750.csv")) %>%
                    mutate(knots = 750, type = "biomass")) %>%
  mutate(matsex = "Immature Female", family = "Delta-gamma")

matfem.DG <- rbind(read.csv(paste0(dir2, "VAST_matfem_DG_abund50.csv")) %>%
                     mutate(knots = 50, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_abund90.csv")) %>%
                     mutate(knots = 90, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_abund120.csv")) %>%
                     mutate(knots = 120, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_abund750.csv")) %>%
                     mutate(knots = 750, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_bio50.csv")) %>%
                     mutate(knots = 50, type = "biomass"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_bio90.csv")) %>%
                     mutate(knots = 90, type = "biomass"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_bio120.csv")) %>%
                     mutate(knots = 120, type = "biomass"),
                   read.csv(paste0(dir2, "VAST_matfem_DG_bio750.csv")) %>%
                     mutate(knots = 750, type = "biomass")) %>%
  mutate(matsex = "Mature Female", family = "Delta-gamma")

male.DG <- rbind(read.csv(paste0(dir2, "VAST_male_DG_abund50.csv")) %>%
                   mutate(knots = 50, type = "abundance"),
                 read.csv(paste0(dir2, "VAST_matfem_DG_abund90.csv")) %>%
                   mutate(knots = 90, type = "abundance"),
                 read.csv(paste0(dir2, "VAST_male_DG_abund120.csv")) %>%
                   mutate(knots = 120, type = "abundance"),
                 read.csv(paste0(dir2, "VAST_male_DG_abund750.csv")) %>%
                   mutate(knots = 750, type = "abundance"),
                 read.csv(paste0(dir2, "VAST_male_DG_bio50.csv")) %>%
                   mutate(knots = 50, type = "biomass"),
                 read.csv(paste0(dir2, "VAST_male_DG_bio90.csv")) %>%
                   mutate(knots = 90, type = "biomass"),
                 read.csv(paste0(dir2, "VAST_male_DG_bio120.csv")) %>%
                   mutate(knots = 120, type = "biomass"),
                 read.csv(paste0(dir2, "VAST_male_DG_bio750.csv")) %>%
                   mutate(knots = 750, type = "biomass")) %>%
  mutate(matsex = "Male", family = "Delta-gamma")

# Bind all 
VAST.24 <- rbind(imfem.TW, matfem.TW, male.TW, imfem.DG, matfem.DG, male.DG) %>%
  rename(SE = "Std..Error.for.Estimate", Year = "Time") %>%
  mutate(Stratum = case_when((Stratum == "Stratum_1") ~ "EBS",
                             (Stratum == "Stratum_2") ~ "East",
                             TRUE ~ "West"),
         model = "VAST",
         CI = 1.96*SE)

dummy <- VAST.24 %>% filter(Year == 2020) %>%
  mutate(Year = 2020, Estimate = NA, SE = NA, Std..Error.for.ln.Estimate. = NA, CI = NA)

VAST.24 <- rbind(dummy, VAST.24 %>% filter(Year != 2020))

VAST.abund <- VAST.24 %>%
  filter(type == "abundance") %>%
  mutate(Estimate = Estimate/1e6, SE = SE/1e6, CI = CI/1e6) 

VAST.bio <- VAST.24 %>%
  filter(type == "biomass") 
# %>%
# mutate(Estimate = Estimate/100, SE = SE/100, CI = CI/100) # Delta gamma is not in kg

### Load survey observations ----------
tan.obs %>%
  group_by(Year, type, matsex) %>%
  reframe(value = sum(value),
          CI = sum(CI)) -> tan.obs3

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
  
  
  ## Males sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
 
  
  # Indices
    # 50
    Allmale.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_All_pre-1982_50-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Male_abundance_All_post-1982_50-Delta_gamma_index.csv"))) %>%
                            rename(abundance = est, Year = year) %>%
                            mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                            mutate(knots = 50, matsex = matsex2)
    
    Allmale.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_All_pre-1982_50-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Male_biomass_All_post-1982_50-Delta_gamma_index.csv"))) %>%
                            rename(biomass = est, Year = year) %>%
                            mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                            mutate(knots = 50, matsex = matsex2)
    
    # 90
    Allmale.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_All_pre-1982_90-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Male_abundance_All_post-1982_90-Delta_gamma_index.csv"))) %>%
                          rename(abundance = est, Year = year) %>%
                          mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                          mutate(knots = 90, matsex = matsex2)
    
    Allmale.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_All_pre-1982_90-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Male_biomass_All_post-1982_90-Delta_gamma_index.csv"))) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                          mutate(knots = 90, matsex = matsex2)
    
    # 120
    Allmale.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_All_pre-1982_120-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Male_abundance_All_post-1982_120-Delta_gamma_index.csv"))) %>%
                          rename(abundance = est, Year = year) %>%
                          mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
                          mutate(knots = 120, matsex = matsex2)
    
    Allmale.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_All_pre-1982_120-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Male_biomass_All_post-1982_120-Delta_gamma_index.csv"))) %>%
                          rename(biomass = est, Year = year) %>%
                          mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                          mutate(knots = 120, matsex = matsex2)
  
  ## Immature Females sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Indices
  # 50
  Allimfem.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_All_pre-1982_50-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Immature Female_abundance_All_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allimfem.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_All_pre-1982_50-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Immature Female_biomass_All_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allimfem.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_All_pre-1982_90-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Immature Female_abundance_All_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allimfem.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_All_pre-1982_90-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Immature Female_biomass_All_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allimfem.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_All_pre-1982_120-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Immature Female_abundance_All_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allimfem.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_All_pre-1982_120-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Immature Female_biomass_All_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  ## Mature Females sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Indices
  # 50
  Allmatfem.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_All_pre-1982_50-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Mature Female_abundance_All_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allmatfem.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_All_pre-1982_50-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Mature Female_biomass_All_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allmatfem.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_All_pre-1982_90-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Mature Female_abundance_All_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allmatfem.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_All_pre-1982_90-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Mature Female_biomass_All_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allmatfem.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_All_pre-1982_120-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Mature Female_abundance_All_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allmatfem.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_All_pre-1982_120-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Mature Female_biomass_All_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
 
  ## Join sdmTMB data and plot indices with obs -----
  rbind(Allmale.index.abund50, Allmale.index.abund90, Allmale.index.abund120, Allimfem.index.abund50,
        Allimfem.index.abund90, Allimfem.index.abund120, Allmatfem.index.abund50,
        Allmatfem.index.abund90, Allmatfem.index.abund120) -> All.abund.index
    
  All.abund.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, abundance = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.abund.index, dummy) %>%
    mutate(model = "sdmTMB") -> All.abund.index
  
  rbind(Allmale.index.bio50, Allmale.index.bio90, Allmale.index.bio120, Allimfem.index.bio50, Allimfem.index.bio90,
        Allimfem.index.bio120, Allmatfem.index.bio50, Allmatfem.index.bio90, Allmatfem.index.bio120) -> All.bio.index
  
  All.bio.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, biomass = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.bio.index, dummy) %>%
    mutate(model = "sdmTMB") -> All.bio.index
  
  write.csv(All.abund.index, paste0(dir, "Output/EBS_Delta_gamma_abundindex.csv"))
  write.csv(All.bio.index, paste0(dir, "Output/EBS_Delta_gamma_bioindex.csv"))
  
  
  
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
  
  ggsave(plot = abund.ind.plot.EBS, "./BAIRDI/Figures/TannerEBS.abundance.index.png", height= 9, width = 7.5, units = "in")
  
  All.abund.index %>%
    filter(matsex == "Male") -> aa
  
  tan.obs3 %>%
    filter(matsex == "Male") -> ss
  
  # Plot indices
  ggplot()+
    geom_ribbon(aa, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4)+
    geom_line(aa, mapping = aes(Year, abundance, color = as.factor(knots)), linewidth = 1)+
    geom_point(ss %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(ss %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          title = element_text(size = 16))
  
  ggsave("./BAIRDI/Figures/TannerEBS.maleabundance.index.png", height= 7, width = 12, units = "in")
  
  All.abund.index %>%
    filter(matsex == "Mature Female") -> aa
  
  tan.obs3 %>%
    filter(matsex == "Mature Female") -> ss
  
  # Plot indices
  # Plot indices
  ggplot()+
    geom_ribbon(aa, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4)+
    geom_line(aa, mapping = aes(Year, abundance, color = as.factor(knots)), linewidth = 1)+
    geom_point(ss %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(ss %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          title = element_text(size = 16))
  ggsave("./BAIRDI/Figures/TannerEBS.matfemabundance.index.png", height= 7, width = 12, units = "in")
  
  
  All.abund.index %>%
    filter(matsex == "Immature Female") -> aa
  
  tan.obs3 %>%
    filter(matsex == "Immature Female") -> ss
  
  # Plot indices
  ggplot()+
    geom_ribbon(aa, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4)+
    geom_line(aa, mapping = aes(Year, abundance, color = as.factor(knots)), linewidth = 1)+
    geom_point(ss %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(ss %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    scale_fill_manual(values = c("salmon", "turquoise", "violet"), labels = c("50", "90", "120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          title = element_text(size = 16))
  
  ggsave("./BAIRDI/Figures/TannerEBS.imfemabundance.index.png", height= 7, width = 12, units = "in")
  
  
  
  
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
  
  ggsave(plot = bio.ind.plot.EBS, "./BAIRDI/Figures/TannerEBS.biomass.index.png", height= 9, width = 7.5, units = "in")
  
  ## Compare sdmTMB and VAST ----
  VAST.abund <- VAST.abund %>% filter(family == "Delta-gamma")
  VAST.bio <- VAST.bio %>% filter(family == "Delta-gamma")
  
  
  ggplot()+
    geom_ribbon(All.abund.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index, mapping = aes(Year, abundance, color = model))+
    geom_ribbon(VAST.abund %>% filter(Stratum == "EBS", knots != 750), 
                  mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.abund %>% filter(Stratum == "EBS", knots !=750), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_grid(factor(knots, levels = c(50, 90, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot()+
    geom_ribbon(All.abund.index %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index %>% filter(knots == 50), mapping = aes(Year, abundance, color = model))+
    geom_ribbon(VAST.abund %>% filter(Stratum == "EBS", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.abund %>% filter(Stratum == "EBS", knots ==50), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> abund.compare
  
  ggsave(plot = abund.compare, "./BAIRDI/Figures/TannerEBS.abundance.sdmTMBVASTindex.png", height= 9, width = 7.5, units = "in")
  
  
  All.abund.index %>%
    filter(matsex == "Male") -> aa
  
  VAST.abund%>%
    filter(matsex == "Male") -> vv
  
  tan.obs3%>%
    filter(matsex == "Male") -> ss
  
  ggplot()+
    geom_ribbon(aa %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(aa %>% filter(knots == 50), mapping = aes(Year, abundance, color = model), linewidth = 1)+
    geom_ribbon(vv%>% filter(Stratum == "EBS", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(vv %>% filter(Stratum == "EBS", knots ==50), mapping = aes(Year, Estimate, color = model), linewidth= 1)+
    geom_point(ss %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(ss %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          title = element_text(size = 16))
  
  ggsave("./BAIRDI/Figures/TannerEBS.maleabundance.sdmTMBVASTindex.png", height= 7, width = 12, units = "in")
  
  All.abund.index %>%
    filter(matsex == "Mature Female") -> aa
  
  VAST.abund%>%
    filter(matsex == "Mature Female") -> vv
  
  tan.obs3%>%
    filter(matsex == "Mature Female") -> ss
  
  ggplot()+
    geom_ribbon(aa %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(aa %>% filter(knots == 50), mapping = aes(Year, abundance, color = model), linewidth = 1)+
    geom_ribbon(vv%>% filter(Stratum == "EBS", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(vv %>% filter(Stratum == "EBS", knots ==50), mapping = aes(Year, Estimate, color = model), linewidth= 1)+
    geom_point(ss %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(ss %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          title = element_text(size = 16))
  
  ggsave("./BAIRDI/Figures/TannerEBS.matfemabundance.sdmTMBVASTindex.png", height= 7, width = 12, units = "in")
  
  All.abund.index %>%
    filter(matsex == "Immature Female") -> aa
  
  VAST.abund%>%
    filter(matsex == "Immature Female") -> vv
  
  tan.obs3%>%
    filter(matsex == "Immature Female") -> ss
  
  ggplot()+
    geom_ribbon(aa %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(aa %>% filter(knots == 50), mapping = aes(Year, abundance, color = model), linewidth = 1)+
    geom_ribbon(vv%>% filter(Stratum == "EBS", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(vv %>% filter(Stratum == "EBS", knots ==50), mapping = aes(Year, Estimate, color = model), linewidth= 1)+
    geom_point(ss %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 1)+
    geom_errorbar(ss %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("EBS Tanner estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          title = element_text(size = 16))
  
  ggsave("./BAIRDI/Figures/TannerEBS.imfemabundance.sdmTMBVASTindex.png", height= 7, width = 12, units = "in")
  
  
  
  
  
  
  
  
  ggplot()+
    geom_ribbon(All.bio.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.bio.index, mapping = aes(Year, biomass, color = model))+
    geom_ribbon(VAST.bio %>% filter(Stratum == "EBS", knots != 750), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.bio %>% filter(Stratum == "EBS", knots !=750), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_grid(factor(knots, levels = c(50, 90, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
    theme_bw()+
    ylab("Biomass (tons)")+
    ggtitle("EBS Tanner estimated biomass") +
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot()+
    geom_ribbon(All.bio.index %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.bio.index %>% filter(knots == 50), mapping = aes(Year, biomass, color = model))+
    geom_ribbon(VAST.bio %>% filter(Stratum == "EBS", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.bio %>% filter(Stratum == "EBS", knots ==50), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    ylab("Biomass (tons)")+
    ggtitle("EBS Tanner estimated biomass") +
    # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> bio.compare
  
  ggsave(plot = bio.compare, "./BAIRDI/Figures/TannerEBS.biomass.sdmTMBVASTindex.png", height= 9, width = 7.5, units = "in")
  
  
### Tanner West ------------
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
  
  
  ## Males sdmTMB-----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "West"
  
  
  # Indices
  # 50
  Allmale.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_West_pre-1982_50-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Male_abundance_West_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allmale.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_West_pre-1982_50-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Male_biomass_West_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allmale.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_West_pre-1982_90-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Male_abundance_West_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allmale.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_West_pre-1982_90-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Male_biomass_West_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allmale.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_West_pre-1982_120-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Male_abundance_West_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allmale.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_West_pre-1982_120-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Male_biomass_West_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  ## Immature Females sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "West"
  
  # Indices
  # 50
  Allimfem.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_West_pre-1982_50-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Immature Female_abundance_West_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allimfem.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_West_pre-1982_50-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Immature Female_biomass_West_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allimfem.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_West_pre-1982_90-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Immature Female_abundance_West_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allimfem.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_West_pre-1982_90-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Immature Female_biomass_West_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allimfem.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_West_pre-1982_120-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Immature Female_abundance_West_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allimfem.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_West_pre-1982_120-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Immature Female_biomass_West_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  ## Mature Females sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "West"
  
  # Indices
  # 50
  Allmatfem.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_West_pre-1982_50-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Mature Female_abundance_West_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allmatfem.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_West_pre-1982_50-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Mature Female_biomass_West_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allmatfem.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_West_pre-1982_90-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Mature Female_abundance_West_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allmatfem.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_West_pre-1982_90-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Mature Female_biomass_West_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allmatfem.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_West_pre-1982_120-Delta_gamma_index.csv")),
                                    read.csv(paste0(dir, "Output/Mature Female_abundance_West_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allmatfem.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_West_pre-1982_120-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Mature Female_biomass_West_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  
  ## Join all data and plot -----
  rbind(Allmale.index.abund50, Allmale.index.abund90, Allmale.index.abund120, Allimfem.index.abund50,
        Allimfem.index.abund90, Allimfem.index.abund120, Allmatfem.index.abund50,
        Allmatfem.index.abund90, Allmatfem.index.abund120) -> All.abund.index
  
  All.abund.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, abundance = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.abund.index, dummy) %>%
    mutate(model = "sdmTMB") -> All.abund.index
  
  rbind(Allmale.index.bio50, Allmale.index.bio90, Allmale.index.bio120, Allimfem.index.bio50, Allimfem.index.bio90,
        Allimfem.index.bio120, Allmatfem.index.bio50, Allmatfem.index.bio90, Allmatfem.index.bio120) -> All.bio.index
  
  All.bio.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, biomass = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.bio.index, dummy) %>%
    mutate(model = "sdmTMB") -> All.bio.index
  
  write.csv(All.abund.index, paste0(dir, "Output/TannerWest_Delta_gamma_abundindex.csv"))
  write.csv(All.bio.index, paste0(dir, "Output/TannerWest_Delta_gamma_bioindex.csv"))
  
  # Survey observations
  tan.obs %>%
    filter(stock == "TannerW") %>%
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
    ggtitle("Tanner West estimated abundance") +
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
  
  ggsave(plot = abund.ind.plot.EBS, "./BAIRDI/Figures/TannerWest.abundance.index.png", height= 9, width = 7.5, units = "in")
  
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
    ggtitle("Tanner West estimated biomass") +
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
  
  ggsave(plot = bio.ind.plot.EBS, "./BAIRDI/Figures/TannerWest.biomass.index.png", height= 9, width = 7.5, units = "in")
  
  ## Compare sdmTMB and VAST ----
  ggplot()+
    geom_ribbon(All.abund.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index, mapping = aes(Year, abundance, color = model))+
    geom_ribbon(VAST.abund %>% filter(Stratum == "West", knots != 750), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.abund %>% filter(Stratum == "West", knots !=750), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_grid(factor(knots, levels = c(50, 90, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner West Tanner estimated abundance") +
    # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot()+
    geom_ribbon(All.abund.index %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index %>% filter(knots == 50), mapping = aes(Year, abundance, color = model))+
    geom_ribbon(VAST.abund %>% filter(Stratum == "West", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.abund %>% filter(Stratum == "West", knots ==50), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner West estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> abund.compare
  
  ggsave(plot = abund.compare, "./BAIRDI/Figures/TannerW.abundance.sdmTMBVASTindex.png", height= 9, width = 7.5, units = "in")
  
  
  ggplot()+
    geom_ribbon(All.bio.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.bio.index, mapping = aes(Year, biomass, color = model))+
    geom_ribbon(VAST.bio %>% filter(Stratum == "West", knots != 750), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.bio %>% filter(Stratum == "West", knots !=750), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_grid(factor(knots, levels = c(50, 90, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
    theme_bw()+
    ylab("Biomass (tons)")+
    ggtitle("Tanner West estimated biomass") +
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot()+
    geom_ribbon(All.bio.index %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.bio.index %>% filter(knots == 50), mapping = aes(Year, biomass, color = model))+
    geom_ribbon(VAST.bio %>% filter(Stratum == "West", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.bio %>% filter(Stratum == "West", knots ==50), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    ylab("Biomass (tons)")+
    ggtitle("Tanner West estimated biomass") +
    # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> bio.compare
  
  ggsave(plot = bio.compare, "./BAIRDI/Figures/TannerW.biomass.sdmTMBVASTindex.png", height= 9, width = 7.5, units = "in")
  
### Tanner East ------------
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
  
  
  ## Males sdm TMB-----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "East"
  
  
  # Indices
  # 50
  Allmale.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_East_pre-1982_50-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Male_abundance_East_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allmale.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_East_pre-1982_50-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Male_biomass_East_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allmale.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_East_pre-1982_90-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Male_abundance_East_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allmale.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_East_pre-1982_90-Delta_gamma_index.csv")),
                               read.csv(paste0(dir, "Output/Male_biomass_East_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allmale.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Male_abundance_East_pre-1982_120-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Male_abundance_East_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allmale.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Male_biomass_East_pre-1982_120-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Male_biomass_East_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  ## Immature Females sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "East"
  
  # Indices
  # 50
  Allimfem.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_East_pre-1982_50-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Immature Female_abundance_East_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allimfem.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_East_pre-1982_50-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Immature Female_biomass_East_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allimfem.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_East_pre-1982_90-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Immature Female_abundance_East_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allimfem.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_East_pre-1982_90-Delta_gamma_index.csv")),
                                read.csv(paste0(dir, "Output/Immature Female_biomass_East_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allimfem.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Immature Female_abundance_East_pre-1982_120-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Immature Female_abundance_East_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allimfem.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Immature Female_biomass_East_pre-1982_120-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Immature Female_biomass_East_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  ## Mature Females sdmTMB -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "East"
  
  # Indices
  # 50
  Allmatfem.index.abund50 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_East_pre-1982_50-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Mature Female_abundance_East_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  Allmatfem.index.bio50 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_East_pre-1982_50-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Mature Female_biomass_East_post-1982_50-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 50, matsex = matsex2)
  
  # 90
  Allmatfem.index.abund90 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_East_pre-1982_90-Delta_gamma_index.csv")),
                                   read.csv(paste0(dir, "Output/Mature Female_abundance_East_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  Allmatfem.index.bio90 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_East_pre-1982_90-Delta_gamma_index.csv")),
                                 read.csv(paste0(dir, "Output/Mature Female_biomass_East_post-1982_90-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 90, matsex = matsex2)
  
  # 120
  Allmatfem.index.abund120 <- rbind(read.csv(paste0(dir, "Output/Mature Female_abundance_East_pre-1982_120-Delta_gamma_index.csv")),
                                    read.csv(paste0(dir, "Output/Mature Female_abundance_East_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(abundance = est, Year = year) %>%
    mutate(abundance = abundance/1e6, lwr = lwr/1e6, upr = upr/1e6)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  Allmatfem.index.bio120 <- rbind(read.csv(paste0(dir, "Output/Mature Female_biomass_East_pre-1982_120-Delta_gamma_index.csv")),
                                  read.csv(paste0(dir, "Output/Mature Female_biomass_East_post-1982_120-Delta_gamma_index.csv"))) %>%
    rename(biomass = est, Year = year) %>%
    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
    mutate(knots = 120, matsex = matsex2)
  
  
  ## Join all data and plot -----
  rbind(Allmale.index.abund50, Allmale.index.abund90, Allmale.index.abund120, Allimfem.index.abund50,
        Allimfem.index.abund90, Allimfem.index.abund120, Allmatfem.index.abund50,
        Allmatfem.index.abund90, Allmatfem.index.abund120) -> All.abund.index
  
  All.abund.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, abundance = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.abund.index, dummy) %>%
    mutate(model = "sdmTMB")-> All.abund.index
  
  rbind(Allmale.index.bio50, Allmale.index.bio90, Allmale.index.bio120, Allimfem.index.bio50, Allimfem.index.bio90,
        Allimfem.index.bio120, Allmatfem.index.bio50, Allmatfem.index.bio90, Allmatfem.index.bio120) -> All.bio.index
  
  All.bio.index %>% 
    filter(Year == 2024) %>%
    mutate(Year = 2020, biomass = NA, lwr = NA, upr = NA, log_est = NA, se = NA) -> dummy
  
  rbind(All.bio.index, dummy) %>%
    mutate(model = "sdmTMB") -> All.bio.index
  
  write.csv(All.abund.index, paste0(dir, "Output/TannerEast_Delta_gamma_abundindex.csv"))
  write.csv(All.bio.index, paste0(dir, "Output/TannerEast_Delta_gamma_bioindex.csv"))
  
  # Survey observations
  tan.obs %>%
    filter(stock == "TannerE") %>%
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
    ggtitle("Tanner East estimated abundance") +
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
  
  ggsave(plot = abund.ind.plot.EBS, "./BAIRDI/Figures/TannerEast.abundance.index.png", height= 9, width = 7.5, units = "in")
  
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
    ggtitle("Tanner East estimated biomass") +
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
  
  ggsave(plot = bio.ind.plot.EBS, "./BAIRDI/Figures/TannerEast.biomass.index.png", height= 9, width = 7.5, units = "in")
  
  ## Compare sdmTMB and VAST ------------
  ggplot()+
    geom_ribbon(All.abund.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index, mapping = aes(Year, abundance, color = model))+
    geom_ribbon(VAST.abund %>% filter(Stratum == "East", knots != 750), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.abund %>% filter(Stratum == "East", knots !=750), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_grid(factor(knots, levels = c(50, 90, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner East Tanner estimated abundance") +
    # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot()+
    geom_ribbon(All.abund.index %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.abund.index %>% filter(knots == 50), mapping = aes(Year, abundance, color = model))+
    geom_ribbon(VAST.abund %>% filter(Stratum == "East", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.abund %>% filter(Stratum == "East", knots ==50), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "abundance"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner East estimated abundance") +
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> abund.compare
  
  ggsave(plot = abund.compare, "./BAIRDI/Figures/TannerE.abundance.sdmTMBVASTindex.png", height= 9, width = 7.5, units = "in")
  
  
  ggplot()+
    geom_ribbon(All.bio.index, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.bio.index, mapping = aes(Year, biomass, color = model))+
    geom_ribbon(VAST.bio %>% filter(Stratum == "East", knots != 750), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.bio %>% filter(Stratum == "East", knots !=750), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_grid(factor(knots, levels = c(50, 90, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
    theme_bw()+
    ylab("Biomass (tons)")+
    ggtitle("Tanner East estimated biomass") +
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot()+
    geom_ribbon(All.bio.index %>% filter(knots == 50), mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
    geom_line(All.bio.index %>% filter(knots == 50), mapping = aes(Year, biomass, color = model))+
    geom_ribbon(VAST.bio %>% filter(Stratum == "East", knots == 50), 
                mapping = aes(x = Year, ymin = Estimate -CI, ymax = Estimate + CI, fill = model), alpha = 0.4)+
    geom_line(VAST.bio %>% filter(Stratum == "East", knots ==50), mapping = aes(Year, Estimate, color = model))+
    geom_point(tan.obs3 %>% filter(type == "biomass"),
               mapping = aes(Year, value), color = "grey20", size = 0.75)+
    geom_errorbar(tan.obs3 %>% filter(type == "biomass"),
                  mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
    facet_wrap(~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y", nrow = 3)+
    theme_bw()+
    scale_color_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    scale_fill_manual(values = c("salmon", "turquoise"), labels = c("sdmTMB", "VAST"), name = "")+
    ylab("Biomass (tons)")+
    ggtitle("Tanner East estimated biomass") +
    # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
    # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16)) -> bio.compare
  
  ggsave(plot = bio.compare, "./BAIRDI/Figures/TannerE.biomass.sdmTMBVASTindex.png", height= 9, width = 7.5, units = "in")
  