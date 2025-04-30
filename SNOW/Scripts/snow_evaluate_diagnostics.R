### PURPOSE ----------------------------------------------------------------------
# To evaluate model-based indices of biomass for EBS snow crab for males >=95mm
# and mature females. Time range is 1980-present.

# Author: Emily Ryznar

# TO DOs:
# 1)

### LOAD LIBRARIES/DATA --------------------------------------------------------
source("./SNOW/Scripts/load_libs_functions.R")

### LOAD FUNCTION --------------------------------------------------------------
evaluate_diagnostics <- function(data, model, category, region, knots, dist){
  
  mod <- paste0(category, "-", region, "-", knots, "-", dist)
  
  # Run sanity check
  print("Sanity")
  sanity_check <- sanity(model) 
  
  # Calculate Dharma residuals
  resid <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, plot = FALSE)
  
  ggplot()+
    theme_bw()+
    geom_point(resid, mapping = aes(expected, observed), size = 2, fill = "black")+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle(mod)+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) -> r.plot
 
  
 resid %>%
    mutate(category = category) -> resids
  
  # Test residuals
  resid <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, return_DHARMa = TRUE)
  
  round(DHARMa::testQuantiles(resid, plot = FALSE)$p.value, 2) -> qq

  round(DHARMa::testDispersion(resid, plot = FALSE)$p.value, 2) -> dd

  round(DHARMa::testOutliers(resid, plot = FALSE)$p.value, 2) -> oo

  round(DHARMa::testZeroInflation(resid, plot = FALSE)$p.value, 2) -> zz

 if(region != "All"){
   data %>% 
     filter(region == region) %>%
     mutate(resids = c(resid$scaledResiduals)) -> data2
 } else{
   data %>% 
     mutate(resids = c(resid$scaledResiduals)) -> data2
 }
  
  # Set spatiotemporal estimator
  if(region != "EBS"){
    sptmp <- "ar1"
  }else{
    sptmp <- "iid"
  }
  
  # visualize residuals across the EBS
  ggplot(data2) +
    #geom_sf(data = shoreline) +
    geom_point(aes(y = lat, x = lon, color = resids), size = 1) +
    scale_color_gradient2(midpoint = 0) +
    labs(y = "Latitude",
         x = "Longitude") +
    theme_gray() +
    ggtitle(mod)+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom") -> res_plot

  ggsave(plot = res_plot, paste0("./SNOW/Figures/DHARMa", mod, "_SPATIAL.png"), height = 9, width = 8.5, units = "in")

  # Calculate log-likelihood
  clust <- sample(seq_len(10), size = nrow(model$data), replace = TRUE)
  
  if(dist == "DG"){
    ll <- sdmTMB_cv(
      data = model$data,
      formula = cpue_kg_km ~ 0 + year_fac,
      spatial = "on",
      time = "year",
      mesh = model$spde,
      spatiotemporal = sptmp,
      silent = FALSE,
      anisotropy = TRUE,
      family = delta_gamma(type = "poisson-link"),
      fold_ids = clust
    )
  } else if(dist == "TW"){
    ll <- sdmTMB_cv(
      data = model$data,
      formula = cpue_kg_km ~ 0 + year_fac,
      spatial = "on",
      time = "year",
      mesh = model$spde,
      spatiotemporal = sptmp,
      silent = FALSE,
      anisotropy = TRUE,
      family = tweedie(link = "log"),
      fold_ids = clust
    )
  }
  
  # Combine evaluation df
  eval.df <- data.frame(category = category,
                        region = region,
                        knots = knots,
                        dist = dist,
                        method = sptmp,
                        loglik = ll$sum_loglik,
                        quantiles = qq,
                        dispersion = dd,
                        outliers = oo,
                        zeroinf = zz)
  
  saveRDS(resids, paste0(dir, "Output/DHARMa_", mod, ".csv"))
  
  return(list(sanity_check_pre, sanity_check_post, resids, eval.df))
}  

make_diagnostic_plots <- function(data, pre.model, post.model, stock2, type, kn, family, method, matsex2){
  
  mod <- paste0(type, "-", kn, "-", family)
  
  # Calculate Dharma residuals
  resid1 <- simulate(pre.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(pre.model, plot = FALSE)
  
  resid2 <- simulate(post.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(post.model, plot= FALSE)
  
  
  # Bind all resids
  rbind(resid1 %>% mutate(period = "<1982"), resid2 %>% mutate(period = "≥1982")) %>%
    mutate(matsex = matsex2) -> all.resids
  
  # Test residuals
  resid1 <- simulate(pre.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(pre.model, return_DHARMa = TRUE)
  resid2 <- simulate(post.model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(post.model, return_DHARMa = TRUE)
  
  # Bind scale resids to pred vs. observed
  scaled <- c(resid1$scaledResiduals, resid2$scaledResiduals)
  
  data.frame(all.resids, scaled.resids = scaled) -> all.resids2
  
  saveRDS(all.resids2, paste0(dir, "Output/DHARMa", stock2, "_", matsex2, "_", mod, ".csv"))
  
  
  # Save DHARMa residual plots
  png(filename = paste0("./BAIRDI/Figures/DHARMa_", matsex2, "_", mod, "_pre1982resplot.png"), width=7, height=5, units="in", res=600)
  
  plot(resid1, title= paste0(matsex2, ", pre-1982, ", mod))
  
  dev.off()
  
  png(filename = paste0("./BAIRDI/Figures/DHARMa_", matsex2, "_", mod, "_post1982resplot.png"), width=7, height=5, units="in", res=600)
  
  plot(resid2, title= paste0(matsex2, ", post-1982, ", mod))
  
  
  dev.off()
  
  # Bind to lat/lon of data
  if(stock2 != "All"){
    data %>% 
      filter(mat.sex == matsex2, stck == stock2) %>%
      mutate(resids = c(resid1$scaledResiduals, resid2$scaledResiduals)) -> data2
  } else{
    data %>% 
      filter(mat.sex == matsex2) %>%
      mutate(resids = c(resid1$scaledResiduals, resid2$scaledResiduals)) -> data2
  }
  
  
  # visualize residuals across the EBS
  ggplot(data2) +
    #geom_sf(data = shoreline) +
    geom_point(aes(y = lat, x = lon, color = resids), size = 1) +
    scale_color_gradient2(midpoint = 0) +
    labs(y = "Latitude",
         x = "Longitude") +
    theme_gray() +
    ggtitle(paste0("EBS ", matsex2," residuals (", mod, ")"))+
    facet_wrap(~year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom") -> res_plot
  
  ggsave(plot = res_plot, paste0("./BAIRDI/Figures/DHARMa_", matsex2, "_", mod, "_SPATIAL.png"), height = 11, width = 10, units = "in")
  
  return(list(all.resids2))
}  

### EVALUATE DIAGNOSTICS ------------
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

evaluate_diagnostics(data, model, category, region, knots = 50, dist = "DG") -> ebs.mf.DG50


## Male EBS 50 knots Tweedie ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_TW_bioTMB.rda"))

evaluate_diagnostics(data, model, category, region, knots = 50, dist = "TW") -> ebs.m.TW50

## Male EBS 120 knots Tweedie ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_TW_bioTMB.rda"))

evaluate_diagnostics(data, model, category, region, knots = 120, dist = "TW") -> ebs.m.TW120

## Male EBS 50 knots Delta Gamma ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_DG_bioTMB.rda"))

evaluate_diagnostics(data, model, category, region, knots = 50, dist = "DG") -> ebs.m.DG50

## Male EBS 120 knots Delta Gamma ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_DG_bioTMB.rda"))

evaluate_diagnostics(data, model, category, region, knots = 120, dist = "DG") -> ebs.m.DG120


## Eval df -----
eval.abund50 <- rbind(ab.males[[4]], ab.imfem[[4]], ab.matfem[[4]])
eval.bio50 <- rbind(bio.males[[4]], bio.imfem[[4]], bio.matfem[[4]])

write.csv(rbind(eval.abund50, eval.bio50), paste0(dir, "Output/eval.TW50.csv"))

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


### EVALUATE MESH ----------------------------------------------------------------


# Mesh is the same across matsex and abund/biomass, so just evaluating male models @ diff knots below

# 50
model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_120_abundTMB.rda"))
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
pre.model <- readRDS(paste0(dir, "Models/bairdi_Male_All_pre-1982_50_abund_DG_IID.rda"))
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
