### PURPOSE ----------------------------------------------------------------------
# To evaluate model-based indices of biomass for EBS snow crab for males >=95mm
# and mature females. Time range is 1980-present.

# Author: Emily Ryznar

# TO DOs:
# 1)

### LOAD LIBRARIES/DATA ----------------------------- ---------------------------
source("./SNOW/Scripts/load_libs_functions.R")

### Set parallel processing for sdmTMB_cv()
parallelly::availableCores()
plan(multisession, workers =6)
options(future.globals.maxSize = 1000 * 1024^2)

### LOAD FUNCTION --------------------------------------------------------------
evaluate_diagnostics <- function(data, model, category, reg, knots, dist){
  
  mod <- paste0(category, "-", reg, "-", knots, "-", dist)
  
  # Run sanity check
  print("Sanity")
  sanity_check <- sanity(model) 
  
  # Calculate Dharma residuals
  resid <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, plot = FALSE)
  
  print("QQ plot")
  ggplot()+
    theme_bw()+
    geom_point(resid, mapping = aes(expected, observed), size = 2, fill = "black")+
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("Expected")+
    xlab("Observed")+
    ggtitle(mod)+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) -> r.plot
 
  ggsave(paste0("./SNOW/Figures/DHARMa_", mod, "_QQplot.png"), height = 8, width = 7, units = "in")
  

  
  print("Testing residuals")
  # Test residuals
  resid <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, return_DHARMa = TRUE)
  
  round(DHARMa::testQuantiles(resid, plot = FALSE)$p.value, 2) -> qq

  round(DHARMa::testDispersion(resid, plot = FALSE)$p.value, 2) -> dd

  round(DHARMa::testOutliers(resid, plot = FALSE)$p.value, 2) -> oo

  round(DHARMa::testZeroInflation(resid, plot = FALSE)$p.value, 2) -> zz

 if(reg != "All"){
   data %>% 
     filter(region == reg) %>%
     mutate(resids = resid$scaledResiduals) -> data2
 } else{
   data %>% 
     mutate(resids = resid$scaledResiduals) -> data2
 }
  
 if(reg != "EBS"){
   sptmp = "ar1"
 } else{
   sptmp = "iid"
 }
 
  
  print("Plotting spatial residuals")
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

  ggsave(plot = res_plot, paste0("./SNOW/Figures/DHARMa_", mod, "_SPATIAL.png"), height = 10, width = 10, units = "in")

  print("Calculating log-likelihood")
  # Calculate log-likelihood
  clust <- sample(seq_len(10), size = nrow(model$data), replace = TRUE)
  
  if(dist == "DG"){
    if(reg != "EBS"){
      ll <- sdmTMB_cv(
        data = model$data,
        formula = cpue_kg_km ~ 0 + year_fac,
        spatial = "on",
        time = "year",
        mesh = model$spde,
        spatiotemporal = "ar1",
        extra_time = c(2020),
        silent = FALSE,
        use_initial_fit = TRUE,
        anisotropy = TRUE,
        family = delta_gamma(type = "poisson-link"),
        fold_ids = clust,
        parallel = TRUE
      )
    } else{
      ll <- sdmTMB_cv(
        data = model$data,
        formula = cpue_kg_km ~ 0 + year_fac,
        spatial = "on",
        time = "year",
        mesh = model$spde,
        spatiotemporal = "iid",
        use_initial_fit = TRUE,
        silent = FALSE,
        anisotropy = TRUE,
        family = delta_gamma(type = "poisson-link"),
        fold_ids = clust,
        parallel = TRUE
        
      )
    }
    
  } else if(dist == "TW"){
    if(reg != "EBS"){
      ll <- sdmTMB_cv(
        data = model$data,
        formula = cpue_kg_km ~ 0 + year_fac,
        spatial = "on",
        time = "year",
        mesh = model$spde,
        spatiotemporal = "ar1",
        extra_time = c(2020),
        use_initial_fit = TRUE,
        silent = FALSE,
        anisotropy = TRUE,
        family = tweedie(link = "log"),
        fold_ids = clust,
        parallel = TRUE
        
      )
    } else{
      ll <- sdmTMB_cv(
        data = model$data,
        formula = cpue_kg_km ~ 0 + year_fac,
        spatial = "on",
        time = "year",
        mesh = model$spde,
        spatiotemporal = "iid",
        use_initial_fit = TRUE,
        silent = FALSE,
        anisotropy = TRUE,
        family = tweedie(link = "log"),
        fold_ids = clust,
        parallel = TRUE
        
      )
    }
    
  }
  
  print("Creating evaluation dataframe")
  # Combine evaluation df
  eval.df <- data.frame(category = category,
                        region = reg,
                        knots = knots,
                        dist = dist,
                        method = sptmp,
                        loglik = ll$sum_loglik,
                        quantiles = qq,
                        dispersion = dd,
                        outliers = oo,
                        zeroinf = zz)
  
  write.csv(eval.df, paste0(dir, "Output/DHARMa_", mod, "_modeleval.csv"))
  
 # saveRDS(resids, paste0(dir, "Output/DHARMa_", mod, ".csv"))
  
  return(list(eval.df))
} 


### EVALUATE DIAGNOSTICS ------------
years <- c(1980:2019, 2021:2024)

# Filter pred_grid (already in UTM), replicate by number of years
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", years) %>%
  rename(lon = X, lat = Y)

## Mature female EBS Delta gamma ----
data <- snow.matfem.cpue
category <- "Mature female"
region <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_120_DG_bioTMB.rda"))


evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "DG") -> ebs.mf.DG50
evaluate_diagnostics(data, model90, category, region, knots = 90, dist = "DG") -> ebs.mf.DG90
evaluate_diagnostics(data, model120, category, region, knots = 120, dist = "DG") -> ebs.mf.DG120

## Mature female EBS TW ----
data <- snow.matfem.cpue
category <- "Mature female"
region <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Mature female_50_TW_bioTMB.rda")) #90 and 120 knots does fit

evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "TW") -> ebs.mf.TW50

## Male95 EBS Delta gamma ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_DG_bioTMB.rda"))


evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "DG") -> ebs.m.DG50
evaluate_diagnostics(data, model90, category, region, knots = 90, dist = "DG") -> ebs.m.DG90
evaluate_diagnostics(data, model120, category, region, knots = 120, dist = "DG") -> ebs.m.DG120

## Male95 EBS Tweedie ----
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"
model50 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_TW_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_90_TW_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_TW_bioTMB.rda"))


evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "TW") -> ebs.m.TW0
evaluate_diagnostics(data, model90, category, region, knots = 90, dist = "TW") -> ebs.m.TW90
evaluate_diagnostics(data, model120, category, region, knots = 120, dist = "TW") -> ebs.m.TW120

## Mature female EBS-NBS Delta gamma ----
data <- snow.matfem.cpue
category <- "Mature female"
region <- "All"
model50 <- readRDS(paste0(dir, "Models/snow_All_Mature female_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_All_Mature female_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_All_Mature female_120_DG_bioTMB.rda"))


evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "DG") -> all.mf.DG50
evaluate_diagnostics(data, model90, category, region, knots = 90, dist = "DG") -> all.mf.DG90
evaluate_diagnostics(data, model120, category, region, knots = 120, dist = "DG") -> all.mf.DG120

## Mature female EBS-NBS TW ----
data <- snow.matfem.cpue
category <- "Mature female"
region <- "All"
model50 <- readRDS(paste0(dir, "Models/snow_All_Mature female_50_TW_bioTMB.rda")) 
model90 <- readRDS(paste0(dir, "Models/snow_All_Mature female_90_TW_bioTMB.rda")) 
model120 <- readRDS(paste0(dir, "Models/snow_All_Mature female_120_TW_bioTMB.rda")) 


evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "TW") ->all.mf.TW50
evaluate_diagnostics(data, model50, category, region, knots = 90, dist = "TW") -> all.mf.TW90
evaluate_diagnostics(data, model50, category, region, knots = 120, dist = "TW") -> all.mf.TW120

## Male95 EBS-NBS Delta gamma ----
data <- snow.male95.cpue
category <- "Male95"
region <- "All"
model50 <- readRDS(paste0(dir, "Models/snow_All_Male95_50_DG_bioTMB.rda"))
model90 <- readRDS(paste0(dir, "Models/snow_All_Male95_90_DG_bioTMB.rda"))
model120 <- readRDS(paste0(dir, "Models/snow_All_Male95_120_DG_bioTMB.rda"))


evaluate_diagnostics(data, model50, category, region, knots = 50, dist = "DG") -> all.m.DG50
evaluate_diagnostics(data, model90, category, region, knots = 90, dist = "DG") -> all.m.DG90
evaluate_diagnostics(data, model120, category, region, knots = 120, dist = "DG") -> all.m.DG120

## Male95 EBS-NBS Tweedie ----
data <- snow.male95.cpue
category <- "Male95"
region <- "All"
model90 <- readRDS(paste0(dir, "Models/snow_All_Male95_90_TW_bioTMB.rda")) # 50 doesn't fit!
model120 <- readRDS(paste0(dir, "Models/snow_All_Male95_120_TW_bioTMB.rda"))


evaluate_diagnostics(data, model90, category, region, knots = 90, dist = "TW") -> all.m.TW90
evaluate_diagnostics(data, model120, category, region, knots = 120, dist = "TW") -> all.m.TW120


### COMBINE EVALUATION OUTPUTS ----------------
files <- list.files(paste0(dir, "Output/"))

model_eval <- files[grep("_modeleval.csv", files)] %>%
  purrr::map_df(~read.csv(paste0(dir, "Output/", .))) %>%
  group_by(category) %>%
  arrange(desc(loglik), .by_group = TRUE)

### EVALUATE MESH --------------
snow.male95.cpue %>%
  filter(region == "EBS") %>%
  group_by(year) %>%
  reframe(N = n()) %>%
  filter(N == min(N)) # mesh vertices need to be less than 305

# EBS
mesh50 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_50_DG_bioTMB.rda"))$spde$mesh
mesh90 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_90_DG_bioTMB.rda"))$spde$mesh
mesh120 <- readRDS(paste0(dir, "Models/snow_EBS_Male95_120_DG_bioTMB.rda"))$spde$mesh

ggplot(snow.male95.cpue %>% filter(year %in% years, region == "EBS")) + 
  inlabru::gg(mesh50) + 
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  theme_bw() +
  ggtitle(paste0("Specified knots = 50, realized knots=", mesh50$n))+
  labs(x = "X", y = "Y")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> ebs.50

ggplot(snow.male95.cpue %>% filter(year %in% years, region == "EBS")) + 
  inlabru::gg(mesh90) + 
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  theme_bw() +
  ggtitle(paste0("Specified knots = 90, realized knots=", mesh90$n))+
  labs(x = "X", y = "Y")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> ebs.90

ggplot(snow.male95.cpue %>% filter(year %in% years, region == "EBS")) + 
  inlabru::gg(mesh120) + 
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  theme_bw() +
  ggtitle(paste0("Specified knots = 120, realized knots=", mesh120$n))+
  labs(x = "X", y = "Y") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> ebs.120

ggpubr::ggarrange(ebs.50, ebs.90, ebs.120, nrow = 3)
ggsave("./SNOW/Figures/snow_EBS_mesh.png", width = 5, height = 9)

# EBS-NBS
mesh50 <- readRDS(paste0(dir, "Models/snow_All_Male95_50_DG_bioTMB.rda"))$spde$mesh
mesh90 <- readRDS(paste0(dir, "Models/snow_All_Male95_90_DG_bioTMB.rda"))$spde$mesh
mesh120 <- readRDS(paste0(dir, "Models/snow_All_Male95_120_DG_bioTMB.rda"))$spde$mesh

ggplot(snow.male95.cpue %>% filter(year %in% years)) + 
  inlabru::gg(mesh50) + 
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  theme_bw() +
  ggtitle(paste0("Specified knots = 50, realized knots=", mesh50$n))+
  labs(x = "X", y = "Y")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> all.50

ggplot(snow.male95.cpue %>% filter(year %in% years)) + 
  inlabru::gg(mesh90) + 
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  theme_bw() +
  ggtitle(paste0("Specified knots = 90, realized knots=", mesh90$n))+
  labs(x = "X", y = "Y")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> all.90

ggplot(snow.male95.cpue %>% filter(year %in% years)) + 
  inlabru::gg(mesh120) + 
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  theme_bw() +
  ggtitle(paste0("Specified knots = 120, realized knots=", mesh120$n))+
  labs(x = "X", y = "Y") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> all.120

ggpubr::ggarrange(all.50, all.90, all.120, nrow = 3)
ggsave("./SNOW/Figures/snow_EBS-NBS_mesh.png", width = 5, height = 9)

