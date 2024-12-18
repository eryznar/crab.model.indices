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
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerW_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
 
  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerW"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerW_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  ## Mature females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerW"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerW_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
 
### TANNER E ------------------------------------------------------------------
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "TannerE"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Male_TannerE_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  ## Immature females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "TannerE"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Immature Female_TannerE_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  ## Mature females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "TannerE"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_120_abundTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_120_abundTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_pre-1988_120_bioTMB.rda")
  post.model <- readRDS("./BAIRDI/Models/bairdi_Mature Female_TannerE_post-1988_120_bioTMB.rda")
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> out
  
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


## 50 knots ----
  ## Males -----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_50_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "TW", "IID", matsex2) -> ab.males
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_50_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, matsex2) -> bio.males
  
  ## Immature Females -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_50_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_50_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, matsex2) -> ab.imfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_50_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_50_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, matsex2) -> bio.imfem
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_50_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_50_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, matsex2) -> ab.matfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_50_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_50_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, matsex2) -> bio.matfem
  
  
  
  ## All abundance QQ plots ----
  rbind(ab.males[[3]], ab.imfem[[3]], ab.matfem[[3]]) %>%
    mutate(period = case_when((period == "<1988") ~ "<1982",
                              TRUE ~ "≥1982")) -> ab.resid
  
  ggplot()+
    theme_bw()+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    geom_point(ab.resid, mapping = aes(expected, observed), size = 1.5)+
    facet_grid(factor(matsex, levels = c("Male", "Mature Female", "Immature Female")) ~ period)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle("Tanner crab abundance DHARMa residual Q-Q plot (knots=50)")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  ggsave(filename = paste0("./BAIRDI/Figures/DHARMa_abundance_EBS_50_QQplot.png"), width=6, height=7, units="in")
  
  # All biomass QQ plots
  rbind(bio.males[[3]], bio.imfem[[3]], bio.matfem[[3]]) -> bio.resid
  
  ggplot()+
    theme_bw()+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    geom_point(ab.resid, mapping = aes(expected, observed), size = 1.5)+
    facet_grid(factor(matsex, levels = c("Male", "Mature Female", "Immature Female")) ~ period)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle("Tanner crab biomass DHARMa residual Q-Q plot (knots=50)")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  ggsave(filename = paste0("./BAIRDI/Figures/DHARMa_biomass_EBS_50_QQplot.png"), width=6, height=7, units="in")
  
  
## 90 knots ----
  ## Males -----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, matsex2) -> ab.males
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, matsex2) -> bio.males

  ## Immature Females -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_90_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_90_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, matsex2) -> ab.imfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_90_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_90_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, matsex2) -> bio.imfem
  
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_90_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_90_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, matsex2) -> ab.matfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_90_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_90_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, matsex2) -> bio.matfem
  
  
  
  ## All abundance QQ plots -----
 rbind(ab.males[[3]], ab.imfem[[3]], ab.matfem[[3]]) %>%
   mutate(period = case_when((period == "<1988") ~ "<1982",
                             TRUE ~ "≥1982")) -> ab.resid
 
 ggplot()+
   theme_bw()+
   geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
   geom_point(ab.resid, mapping = aes(expected, observed), size = 1.5)+
   facet_grid(factor(matsex, levels = c("Male", "Mature Female", "Immature Female")) ~ period)+
   ylab("Expected")+
   xlab("Observed")+
   ggtitle("Tanner crab abundance DHARMa residual Q-Q plot (knots=90)")+
   theme(axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))
 
 ggsave(filename = paste0("./BAIRDI/Figures/DHARMa_abundance_EBS_90_QQplot.png"), width=6, height=7, units="in")
 
 # All biomass QQ plots
 rbind(bio.males[[3]], bio.imfem[[3]], bio.matfem[[3]]) -> bio.resid
 
 ggplot()+
   theme_bw()+
   geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
   geom_point(ab.resid, mapping = aes(expected, observed), size = 1.5)+
   facet_grid(factor(matsex, levels = c("Male", "Mature Female", "Immature Female")) ~ period)+
   ylab("Expected")+
   xlab("Observed")+
   ggtitle("Tanner crab biomass DHARMa residual Q-Q plot (knots=90)")+
   theme(axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))
 
 ggsave(filename = paste0("./BAIRDI/Figures/DHARMa_biomass_EBS_90_QQplot.png"), width=6, height=7, units="in")
 
 
## 120 knots ----
 ## Males -----  
 data <- tan.cpue2
 matsex2 <- "Male"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_abundTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_120_abundTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> ab.males
 
 # Biomass
 type <- "biomass"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_bioTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_120_bioTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> bio.males
 
 ## Immature Females -----  
 data <- tan.cpue2
 matsex2 <- "Immature Female"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_120_abundTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_120_abundTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> ab.imfem
 
 # Biomass
 type <- "biomass"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_120_bioTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_120_bioTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> bio.imfem
 
 
 ## Mature Females -----  
 data <- tan.cpue2
 matsex2 <- "Mature Female"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_120_abundTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_120_abundTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> ab.matfem
 
 # Biomass
 type <- "biomass"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_120_bioTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_120_bioTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, matsex2) -> bio.matfem
 
 ## All abundance QQ plots -----
 rbind(ab.males[[3]], ab.imfem[[3]], ab.matfem[[3]]) %>%
   mutate(period = case_when((period == "<1988") ~ "<1982",
                             TRUE ~ "≥1982")) -> ab.resid
 
 ggplot()+
   theme_bw()+
   geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
   geom_point(ab.resid, mapping = aes(expected, observed), size = 1.5)+
   facet_grid(factor(matsex, levels = c("Male", "Mature Female", "Immature Female")) ~ period)+
   ylab("Expected")+
   xlab("Observed")+
   ggtitle("Tanner crab abundance DHARMa residual Q-Q plot (knots=120)")+
   theme(axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))
 
 ggsave(filename = paste0("./BAIRDI/Figures/DHARMa_abundance_EBS_120_QQplot.png"), width=6, height=7, units="in")
 
 # All biomass QQ plots
 rbind(bio.males[[3]], bio.imfem[[3]], bio.matfem[[3]]) -> bio.resid
 
 ggplot()+
   theme_bw()+
   geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
   geom_point(ab.resid, mapping = aes(expected, observed), size = 1.5)+
   facet_grid(factor(matsex, levels = c("Male", "Mature Female", "Immature Female")) ~ period)+
   ylab("Expected")+
   xlab("Observed")+
   ggtitle("Tanner crab biomass DHARMa residual Q-Q plot (knots=120)")+
   theme(axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))
 
 ggsave(filename = paste0("./BAIRDI/Figures/DHARMa_biomass_EBS_120_QQplot.png"), width=6, height=7, units="in")
 
 
 
### EVALUATE MESH ----------------------------------------------------------------
tan.cpue2 %>%
   mutate(period = case_when((year <1982) ~ "<1982",
                             TRUE ~ "≥1982")) %>%
   group_by(year, period, mat.sex) %>%
   reframe(N = n()) %>%
   group_by(period, mat.sex) %>%
   reframe(min.N = min(N)) 

    # Mesh is the same across matsex and abund/biomass, so just evaluating male models @ diff knots below
 
    # 120
    pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_abundTMB.rda"))
    post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_120_abundTMB.rda"))
    
    pre.model$spde$mesh$n # needs to be <136
    post.model$spde$mesh$n # needs to be <342
    
    ggplot(tan.cpue2 %>% filter(year %in% 1975:1981)) + 
      inlabru::gg(pre.model$spde$mesh) + 
      geom_point(aes(x = lon, y = lat)) +
      theme_bw() +
      ggtitle(paste0("<1982 mesh (knots=", pre.model$spde$mesh$n, ")"))+
      labs(x = "X", y = "Y") -> mesh.1
    
    ggplot(tan.cpue2 %>% filter(year %in% 1982:2024)) + 
      inlabru::gg(post.model$spde$mesh) + 
      geom_point(aes(x = lon, y = lat)) +
      theme_bw() +
      ggtitle(paste0("≥1982 mesh (knots=", post.model$spde$mesh$n, ")"))+
    labs(x = "X", y = "Y") -> mesh.2
    
    
    ggpubr::ggarrange(mesh.1, mesh.2, nrow = 2)
    ggsave("./BAIRDI/Figures/mesh120.png", height = 8, width = 6, units = "in")
    
    
    # 90
    pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abundTMB.rda"))
    post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abundTMB.rda"))
    
    pre.model$spde$mesh$n # needs to be <136
    post.model$spde$mesh$n # needs to be <342
    
    ggplot(tan.cpue2 %>% filter(year %in% 1975:1981)) + 
      inlabru::gg(pre.model$spde$mesh) + 
      geom_point(aes(x = lon, y = lat)) +
      theme_bw() +
      ggtitle(paste0("<1982 mesh (knots=", pre.model$spde$mesh$n, ")"))+
    labs(x = "X", y = "Y") -> mesh.1
    
    ggplot(tan.cpue2 %>% filter(year %in% 1982:2024)) + 
      inlabru::gg(post.model$spde$mesh) + 
      geom_point(aes(x = lon, y = lat)) +
      theme_bw() +
      ggtitle(paste0("≥1982 mesh (knots=", post.model$spde$mesh$n, ")"))+
    labs(x = "X", y = "Y") -> mesh.2
    
    ggpubr::ggarrange(mesh.1, mesh.2, nrow = 2)
    ggsave("./BAIRDI/Figures/mesh90.png", height = 8, width = 6, units = "in")
    
    # 50
    pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abundTMB.rda"))
    post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_50_abundTMB.rda"))
    
    pre.model$spde$mesh$n # needs to be <136
    post.model$spde$mesh$n # needs to be <342
  
    ggplot(tan.cpue2 %>% filter(year %in% 1975:1981)) + 
      inlabru::gg(pre.model$spde$mesh) + 
      geom_point(aes(x = lon, y = lat)) +
      theme_bw() +
      ggtitle(paste0("<1982 mesh (knots=", pre.model$spde$mesh$n, ")"))+
    labs(x = "X", y = "Y") -> mesh.1
    
    ggplot(tan.cpue2 %>% filter(year %in% 1982:2024)) + 
      inlabru::gg(post.model$spde$mesh) + 
      geom_point(aes(x = lon, y = lat)) +
      theme_bw() +
      ggtitle(paste0("≥1982 mesh (knots=", post.model$spde$mesh$n, ")"))+
    labs(x = "X", y = "Y") -> mesh.2
    
    ggpubr::ggarrange(mesh.1, mesh.2, nrow = 2)
    ggsave("./BAIRDI/Figures/mesh50.png", height = 8, width = 6, units = "in")
   