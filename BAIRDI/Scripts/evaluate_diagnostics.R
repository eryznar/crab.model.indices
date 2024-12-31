### PURPOSE ----------------------------------------------------------------------
# To generate model-based indices of abundance and biomass for EBS Tanner crab for all males, immature females, and 
# mature females. Minimum size is 25 mm CW and the time range is 1975-present. Try fitting separate models for 1975-1981 and 
# 1982+. 

# Author: Emily Ryznar

# TO DOs:
# 1) Look at residuals
# 2) Add in scripts to load new survey data (CPUE, BIO/ABUND) and process each year (CPUE script is in TECHMEMONEW)

### LOAD LIBRARIES/DATA --------------------------------------------------------
source("./BAIRDI/Scripts/load_libs_functions.R")

### LOAD FUNCTION --------------------------------------------------------------

evaluate_diagnostics <- function(data, pre.model, post.model, stock2, type, kn, family, method, matsex2){
  
  mod <- paste0(type, "-", kn, "-", family)
  
  # Run sanity check
  print("Pre-1982 model")
  sanity_check_pre <- sanity(pre.model) 
  print("Post-1982 model")
  sanity_check_post <- sanity(post.model) 
  
  # Calculate Dharma residuals
  resid1 <- simulate(pre.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(pre.model, plot = FALSE)
  
  ggplot()+
    theme_bw()+
    geom_point(resid1, mapping = aes(expected, observed), size = 2, fill = "black")+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle("<1982")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) -> r1.plot
  # 
  # png(filename = paste0("./BAIRDI/Figures/DHARMa_pre1988_", type, "_", stock2, "_", matsex2, "_", kn, ".png"), width=7, height=5, units="in", res=600)
  # 
  # plot(resid1, title= paste0("DHARMa residuals (Pre-1988, ", type, ",", stock2, " ", matsex2, ", kn=", kn, ")"))
  # 
  # dev.off()
  
  resid2 <- simulate(post.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(post.model, plot= FALSE)
  
  ggplot()+
    theme_bw()+
    geom_point(resid2, mapping = aes(expected, observed), size = 2, fill = "black")+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle("≥1982")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) -> r2.plot
  
  rbind(resid1 %>% mutate(period = "<1982"), resid2 %>% mutate(period = "≥1982")) %>%
    mutate(matsex = matsex2) -> all.resids
  
  # Test residuals
  resid1 <- simulate(pre.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(pre.model, return_DHARMa = TRUE)
  resid2 <- simulate(post.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(post.model, return_DHARMa = TRUE)
  
  round(DHARMa::testQuantiles(resid1, plot = FALSE)$p.value, 2) -> qq1
  round(DHARMa::testQuantiles(resid2, plot = FALSE)$p.value, 2) -> qq2
  
  round(DHARMa::testDispersion(resid1, plot = FALSE)$p.value, 2) -> dd1
  round(DHARMa::testDispersion(resid2, plot = FALSE)$p.value,2) -> dd2
  
  round(DHARMa::testOutliers(resid1, plot = FALSE)$p.value, 2) -> oo1
  round(DHARMa::testOutliers(resid2, plot = FALSE)$p.value, 2) -> oo2
  
  round(DHARMa::testZeroInflation(resid1, plot = FALSE)$p.value, 2) -> zz1
  round(DHARMa::testZeroInflation(resid2, plot = FALSE)$p.value, 2) -> zz2
  
  
  if(stock2 != "All"){
    data %>% 
      filter(mat.sex == matsex2, stck == stock2) %>%
      mutate(resids = c(resid1$scaledResiduals, resid2$scaledResiduals)) -> data2
  } else{
    data %>% 
      filter(mat.sex == matsex2) %>%
      mutate(resids = c(resid1$scaledResiduals, resid2$scaledResiduals)) -> data2
  }
  
  
  # # visualize residuals across the EBS
  # ggplot(data2) + 
  #   #geom_sf(data = shoreline) +
  #   geom_point(aes(y = lat, x = lon, color = resids), size = 1) +
  #   scale_color_gradient2(midpoint = 0) + 
  #   labs(y = "Latitude",
  #        x = "Longitude") +
  #   theme_gray() + 
  #   ggtitle(paste0(stock2," ", matsex2," ", type, " residuals (kn=", kn, ")"))+
  #   facet_wrap(~year)+
  #   theme(axis.title = element_text(size = 10),
  #         legend.position = "bottom") -> res_plot
  # 
  # ggsave(plot = res_plot, paste0("./BAIRDI/Figures/DHARMa", stock2, "_", matsex2, "_", mod, "_SPATIAL.png"), height = 9, width = 8.5, units = "in")
  
  # Calculate log-likelihood
  clust <- sample(seq_len(10), size = nrow(pre.model$data), replace = TRUE)
  
  pre.ll <- sdmTMB_cv(
    data = pre.model$data,
    formula = cpue_km ~ 0 + year_fac,
    spatial = "on",
    time = "year",
    mesh = pre.model$spde,
    spatiotemporal = "iid",
    silent = FALSE,
    anisotropy = TRUE,
    family = tweedie(link = "log"),
    fold_ids = clust
  )
  
  clust <- sample(seq_len(10), size = nrow(post.model$data), replace = TRUE)
  
  post.ll <- sdmTMB_cv(
    data = post.model$data,
    formula = cpue_km ~ 0 + year_fac,
    spatial = "on",
    time = "year",
    mesh = post.model$spde,
    spatiotemporal = "iid",
    silent = FALSE,
    anisotropy = TRUE,
    family = tweedie(link = "log"),
    fold_ids = clust
  )
  
  eval.df <- data.frame(matsex = rep(matsex2,2),
                        type = rep(type,2),
                        knots = rep(kn,2),
                        family = rep(family, 2),
                        method = rep(method,2),
                        period = c("<1982", "≥1982"),
                        loglik = c(pre.ll$sum_loglik, post.ll$sum_loglik),
                        quantiles = c(qq1, qq2),
                        dispersion = c(dd1, dd2),
                        outliers = c(oo1, oo2),
                        zeroinf = c(zz1, zz2))
  
  saveRDS(all.resids, paste0(dir, "Output/DHARMa", stock2, "_", matsex2, "_", mod, ".csv"))
  
  #return(list(sanity_check_pre, sanity_check_post, all.resids))
  
  return(list(sanity_check_pre, sanity_check_post, all.resids, eval.df))
}  
 


### EVALUATE MODELS EBS-wide ------------
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


## TWEEDIE IID 50 knots ----
  ## Males -----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_50_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "Tweedie", "IID", matsex2) -> ab.males
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_50_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "Tweedie", "IID", matsex2) -> bio.males
  
  ## Immature Females -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_50_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_50_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "Tweedie", "IID", matsex2) -> ab.imfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_50_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_50_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "Tweedie", "IID", matsex2) -> bio.imfem
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_50_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_50_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "Tweedie", "IID", matsex2) -> ab.matfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_50_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_50_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 50, "Tweedie", "IID", matsex2) -> bio.matfem
  
  ## Eval df -----
  eval.abund50 <- rbind(ab.males[[4]], ab.imfem[[4]], ab.matfem[[4]])
  eval.bio50 <- rbind(bio.males[[4]], bio.imfem[[4]], bio.matfem[[4]])
  
  write.csv(rbind(eval.abund50, eval.bio50), paste0(dir, "Output/eval.TW50.csv"))
  
  ## All abundance QQ plots ----
  rbind(ab.males[[3]], ab.imfem[[3]], ab.matfem[[3]]) %>%
    mutate(period = case_when((period == "<1982") ~ "<1982",
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
  
  
## TWEEDIE IID 90 knots ----
  ## Males -----  
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "IID", matsex2) -> ab.males
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "IID", matsex2) -> bio.males

  ## Immature Females -----  
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_90_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_90_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "IID", matsex2) -> ab.imfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_90_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_90_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "IID", matsex2) -> bio.imfem
  
  
  ## Mature Females -----  
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_90_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_90_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "IID", matsex2) -> ab.matfem
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_90_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_90_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "IID", matsex2) -> bio.matfem
  
  ## Eval df -----
  eval.abund90 <- rbind(ab.males[[4]], ab.imfem[[4]], ab.matfem[[4]])
  eval.bio90 <- rbind(bio.males[[4]], bio.imfem[[4]], bio.matfem[[4]])

  write.csv(rbind(eval.abund90, eval.bio90), paste0(dir, "Output/eval.TW90.csv"))
  
  ## All abundance QQ plots -----
 rbind(ab.males[[3]], ab.imfem[[3]], ab.matfem[[3]]) %>%
   mutate(period = case_when((period == "<1982") ~ "<1982",
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
 
 
## TWEEDIE IID 120 knots ----
 ## Males -----  
 data <- tan.cpue2
 matsex2 <- "Male"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_abundTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_120_abundTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, "Tweedie", "IID", matsex2) -> ab.males
 
 # Biomass
 type <- "biomass"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_bioTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_120_bioTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, "Tweedie", "IID", matsex2) -> bio.males
 
 ## Immature Females -----  
 data <- tan.cpue2
 matsex2 <- "Immature Female"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_120_abundTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_120_abundTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, "Tweedie", "IID", matsex2) -> ab.imfem
 
 # Biomass
 type <- "biomass"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_120_bioTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_120_bioTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, "Tweedie", "IID", matsex2) -> bio.imfem
 
 
 ## Mature Females -----  
 data <- tan.cpue2
 matsex2 <- "Mature Female"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_120_abundTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_120_abundTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, "TW", "IID", matsex2) -> ab.matfem
 
 # Biomass
 type <- "biomass"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_120_bioTMB.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_120_bioTMB.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 120, "Tweedie", "IID", matsex2) -> bio.matfem
 
 ## Eval df -----
 eval.abund120 <- rbind(ab.males[[4]], ab.imfem[[4]], ab.matfem[[4]])
 eval.bio120 <- rbind(bio.males[[4]], bio.imfem[[4]], bio.matfem[[4]])
 
 write.csv(rbind(eval.abund120, eval.bio120), paste0(dir, "Output/eval.TW120.csv"))
 
 
 ## All abundance QQ plots -----
 rbind(ab.males[[3]], ab.imfem[[3]], ab.matfem[[3]]) %>%
   mutate(period = case_when((period == "<1982") ~ "<1982",
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
 
 
 
 
## DELTA GAMMA IID 90 knots ----
 ## Males -----
 data <- tan.cpue2
 matsex2 <- "Male"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abund_DG_IID.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_DG_IID.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Delta_gamma", "IID", matsex2) -> ab.males.DG
 
## DELTA LOGNORMAL IID 90 knots ----
 ## Males -----
 data <- tan.cpue2
 matsex2 <- "Male"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abund_DLN_IID.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_DLN_IID.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Delta_lognormal", "IID", matsex2) -> ab.males.DLN
 
## TWEEDIE AR1 90 knots ----
 ## Males -----
 data <- tan.cpue2
 matsex2 <- "Male"
 stock2 <- "All"
 
 # Abundance
 type <- "abundance"
 pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_abund_T_AR1.rda"))
 post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_abund_T_AR1.rda"))
 
 evaluate_diagnostics(data, pre.model, post.model, stock2, type, knots = 90, "Tweedie", "AR1", matsex2) -> ab.males.TAR1
 
 
## Bind all evaluation dfs -----
eval.abund <- rbind(eval.abund50, eval.abund90, eval.abund120, ab.males.DG[[4]], ab.males.DLN[[4]], ab.males.TAR1[[4]])
eval.bio <- rbind(eval.bio50, eval.bio90, eval.bio120)

write.csv(eval.abund, paste0(dir, "Output/model_eval_abund.csv"))
write.csv(eval.bio, paste0(dir, "Output/model_eval_bio.csv"))

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
   
## DELTA GAMMA IID 90 kn ----
  ## Males biomass-----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_90_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_90_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, 90, "Delta_gamma", "IID", matsex2) -> bio.males.DG
  
  ## Immature Females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_90_DG_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_90_DG_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 90, "Delta_gamma", "IID", matsex2) -> ab.imfem.DG
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_90_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_90_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 90, "Delta_gamma", "IID", matsex2) -> bio.imfem.DG
  
  ## Mature Females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_90_DG_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_90_DG_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 90, "Delta_gamma", "IID", matsex2) -> ab.matfem.DG
  
  # Abundance
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_90_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_90_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 90, "Delta_gamma", "IID", matsex2) -> bio.matfem.DG
  
## DELTA GAMMA IID 50 kn ----
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_50_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, 50, "Delta_gamma", "IID", matsex2) -> bio.males.DG
  
  ## Immature Females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_50_DG_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_50_DG_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 50, "Delta_gamma", "IID", matsex2) -> ab.imfem.DG5
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_50_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_50_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 50, "Delta_gamma", "IID", matsex2) -> bio.imfem.DG5
  
  ## Mature Females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_50_DG_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_50_DG_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 50, "Delta_gamma", "IID", matsex2) -> ab.matfem.DG5
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_50_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_50_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 50, "Delta_gamma", "IID", matsex2) -> bio.matfem.DG5
  
## DELTA GAMMA IID 120 kn ----
  ## Males -----
  data <- tan.cpue2
  matsex2 <- "Male"
  stock2 <- "All"
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_post-1982_120_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, 120, "Delta_gamma", "IID", matsex2) -> bio.males.DG12
  
  ## Immature Females -----
  data <- tan.cpue2
  matsex2 <- "Immature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_120_DG_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_120_DG_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 120, "Delta_gamma", "IID", matsex2) -> ab.imfem.DG12
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_pre-1982_120_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Immature Female_All_post-1982_120_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 120, "Delta_gamma", "IID", matsex2) -> bio.imfem.DG12
  
  ## Mature Females -----
  data <- tan.cpue2
  matsex2 <- "Mature Female"
  stock2 <- "All"
  
  # Abundance
  type <- "abundance"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_120_DG_abundTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_120_DG_abundTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 120, "Delta_gamma", "IID", matsex2) -> ab.matfem.DG12
  
  # Biomass
  type <- "biomass"
  pre.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_pre-1982_120_DG_bioTMB.rda"))
  post.model <- readRDS(paste0(dir, "Models/bairdi_Mature Female_All_post-1982_120_DG_bioTMB.rda"))
  
  evaluate_diagnostics(data, pre.model, post.model, stock2, type, kn = 120, "Delta_gamma", "IID", matsex2) -> bio.matfem.DG12
  
  rbind(bio.males.DG[[4]], ab.imfem.DG[[4]], bio.imfem.DG[[4]], ab.matfem.DG[[4]], 
        bio.matfem.DG[[4]], ab.imfem.DG5[[4]], 
        bio.imfem.DG5[[4]], ab.matfem.DG5[[4]], bio.matfem.DG5[[4]], ab.imfem.DG12[[4]], 
        bio.imfem.DG12[[4]], ab.matfem.DG12[[4]], bio.matfem.DG12[[4]]) -> dg.eval
  
  write.csv(dg.eval, paste0(dir, "Output/delta_gammafem.eval.csv"))