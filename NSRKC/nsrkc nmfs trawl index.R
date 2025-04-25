# ************************************************************************************************************
# Generating a spatiotemporal model-based index for Norton Sound red king crab in the NMFS trawl surveys
# (EBS and NBS)
# December 2024
# Caitlin Stern
# ************************************************************************************************************

# **************************************************************************************************************
# load libraries, set plot preferences ----
# **************************************************************************************************************

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sp)

cur_yr.ns <- 2024

plotdir.ns <- paste0(here::here(), "/NSRKC/plots")
modeldir.ns <- paste0(here::here(), "/NSRKC/models")

# set plot appearance
theme_sleek <- function(base_size = 12) {
  
  half_line <- base_size/2
  
  theme_light(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      #axis.text = element_text(colour = "grey30"),
      #axis.title = element_text(colour = "grey30"),
      #legend.title = element_text(colour = "grey30"),#, size = rel(0.9)
      panel.border = element_rect(fill = NA),#, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      #legend.text = element_text(size = rel(0.7)),#, colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)#,
      #plot.title = element_text(colour = "grey30"),#, size = rel(1)
      #plot.subtitle = element_text(colour = "grey30")#, size = rel(.85)
    )
  
}

theme_set(theme_sleek())
cbpalette <- colorRampPalette(colors = c("#009E73", "#0072B2","#E69F00" , "#56B4E9", "#D55E00", "#CC79A7","#F0E442", "black", "grey"))(9)

# **************************************************************************************************************
# read in grid used in VAST to use for predictions and convert to UTM ----
# **************************************************************************************************************

grid2km.ns <- readRDS(paste0(here::here(), "/NSRKC/data/NSRKC_Extrapolation_areas/NSRKC_2km_grid.rds")) %>%
  # get UTM zone
  mutate(zone = (floor((Lon + 180)/6) %% 60) + 1) 

unique(grid2km.ns$zone)

predgrid_utm_z2.ns <- grid2km.ns %>% filter(zone == 2)
predgrid_utm_z3.ns <- grid2km.ns %>% filter(zone == 3)
predgrid_utm_z4.ns <- grid2km.ns %>% filter(zone == 4)

u2.ns <- predgrid_utm_z2.ns
get_crs(u2.ns, c("Lon", "Lat"))
u2u.ns <- add_utm_columns(u2.ns, c("Lon", "Lat"))

u3.ns <- predgrid_utm_z3.ns
get_crs(u3.ns, c("Lon", "Lat"))
u3u.ns <- add_utm_columns(u3.ns, c("Lon", "Lat"))

u4.ns <- predgrid_utm_z4.ns
get_crs(u4.ns, c("Lon", "Lat"))
u4u.ns <- add_utm_columns(u4.ns, c("Lon", "Lat"))

predgrid_utm.ns <- rbind(u2u.ns, u3u.ns, u4u.ns) %>%
  mutate(X = X / 1000, Y = Y / 1000) # convert UTM coordinates from meter to kilometers to reduce computing needs


# **************************************************************************************************************
# read in and process survey data ----
# **************************************************************************************************************

# using the data sets that are used in the assessment

# read in NMFS haul data
# comment in csv: "Stations without ADFG station indicate outside of NS area" - filter out rows with missing ADFG station?
haul_rkc_nmfs <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NMFS_1976_1991_Haul_data.csv'), skip = 6)

# read in NBS haul data
haul_rkc_nbs2010 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NOAA_NBS_Haul_2010.csv'), skip = 6)
haul_rkc_nbs2017 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NOAA_NBS_Haul_2017.csv'), skip = 6)
haul_rkc_nbs2019 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NOAA_NBS_Haul_2019.csv'), skip = 6)
haul_rkc_nbs2021 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NOAA_NBS_Haul_2021.csv'), skip = 6)
haul_rkc_nbs2022 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NOAA_NBS_Haul_2022.csv'), skip = 6)
haul_rkc_nbs2023 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/NOAA_NBS_Haul_2023.csv'), skip = 6)

# read in ADFG haul data
haul_rkc_adfg1996 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_1996.csv'), skip = 6)
haul_rkc_adfg1999 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_1999.csv'), skip = 6)
haul_rkc_adfg2002 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2002.csv'), skip = 6)
haul_rkc_adfg2006 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2006.csv'), skip = 6)
haul_rkc_adfg2008 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2008.csv'), skip = 6)
haul_rkc_adfg2011 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2011.csv'), skip = 6)
haul_rkc_adfg2014 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2014.csv'), skip = 6)
haul_rkc_adfg2017 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2017.csv'), skip = 6)
haul_rkc_adfg2018 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2018.csv'), skip = 6)
haul_rkc_adfg2019 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2019.csv'), skip = 6)
haul_rkc_adfg2020 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2020.csv'), skip = 6)
haul_rkc_adfg2021 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2021.csv'), skip = 6)
haul_rkc_adfg2023 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2023.csv'), skip = 6)
haul_rkc_adfg2024 <- read.csv(paste0(here::here(), '/NSRKC/data/data from assessment/ADFG_Haul_2024.csv'), skip = 6)



# **************************************************************************************************************
# make SPDE mesh ----
# **************************************************************************************************************

NS_spde <- make_mesh(nsrkc_kgkm, xy_cols = c("X","Y"), n_knots = 120, type = "kmeans")
NS_spde_50kn <- make_mesh(nsrkc_kgkm, xy_cols = c("X","Y"), n_knots = 50, type = "kmeans")

plot(NS_spde)
NS_spde$mesh$n

plot(NS_spde_50kn)
NS_spde_50kn$mesh$n


# **************************************************************************************************************
# fit and check models ----
# **************************************************************************************************************

# fit a GLMM with spatiotemporal = AR1 and 120 knots
m.nsrkc.ar1 <- sdmTMB(
  data = rkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = NS_spde, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 120 knots
m.nsrkc.iid <- sdmTMB(
  data = rkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = NS_spde, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = AR1 and 50 knots
m.nsrkc.50kn <- sdmTMB(
  data = rkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = NS_spde_50kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# save the fitted models
saveRDS(m.nsrkc.ar1, file = paste0(modeldir.ns, "/m_nsrkc_ar1_2024_12_02.RDS"))
saveRDS(m.nsrkc.iid, file = paste0(modeldir.ns, "/m_nsrkc_iid_2024_12_02.RDS"))
saveRDS(m.nsrkc.50kn, file = paste0(modeldir.ns, "/m_nsrkc_50kn_2024_12_02.RDS"))

# print the model fit
m.nsrkc.ar1
m.nsrkc.ar1$sd_report

# view parameters
tidy.nsrkc.ar1 <- tidy(m.nsrkc.ar1, conf.int = TRUE)
tidy.nsrkc.ran.ar1 <- tidy(m.nsrkc.ar1 , effects = "ran_pars", conf.int = TRUE)

# run sanity checks
sanity(m.nsrkc.ar1)
sanity(m.nsrkc.iid)
sanity(m.nsrkc.50kn)


# **************************************************************************************************************
# model residuals ----
# **************************************************************************************************************

# randomized quantile residuals
rkc_kgkm_ar1$resids <- residuals(m.nsrkc.ar1) # randomized quantile residuals

hist(rkc_kgkm_ar1$resids)
qqnorm(rkc_kgkm_ar1$resids)
abline(a=0, b=1)

# DHARMa residuals
m.nsrkc.ar1.resid <- simulate(update(m.nsrkc.ar1, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.ar1, return_DHARMa = TRUE)
plot(m.nsrkc.ar1.resid)

m.nsrkc.iid.resid <- simulate(update(m.nsrkc.iid, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.iid, return_DHARMa = TRUE)
plot(m.nsrkc.iid.resid)

m.nsrkc.50kn.resid <- simulate(update(m.nsrkc.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.50kn, return_DHARMa = TRUE)
plot(m.nsrkc.50kn.resid)



# **************************************************************************************************************
# generate predictions ----
# **************************************************************************************************************

# replicate grid across all years
pred_grid.ns <- replicate_df(predgrid_utm.ns, "year_f", unique(rkc_kgkm_utm$year_f))
pred_grid.ns$year <- as.integer(as.character(factor(pred_grid.ns$year_f)))
dplyr::glimpse(pred_grid.ns)

pred_grid.ns <- replicate_df(predgrid_utm.ns, "SURVEY_YEAR", unique(rkc_kgkm_utm$SURVEY_YEAR))
pred_grid.ns$year_f <- factor(pred_grid.ns$SURVEY_YEAR)
dplyr::glimpse(pred_grid.ns)

# predictions on new data
predictions.ar1.ns <- predict(m.nsrkc.ar1, newdata = pred_grid.ns, return_tmb_object = T)
predictions.50kn.ns <- predict(m.nsrkc.50kn, newdata = pred_grid.ns, return_tmb_object = T)
predictions.iid.ns <- predict(m.nsrkc.iid, newdata = pred_grid.ns, return_tmb_object = T)

# save the predictions
saveRDS(predictions.ar1.ns, file = paste0(modeldir.ns, "/m_nsrkc_predictions_ar1_2024_12_02.RDS"))
saveRDS(predictions.50kn.ns, file = paste0(modeldir.ns, "/m_nsrkc_predictions_50kn_2024_12_02.RDS"))
saveRDS(predictions.iid.ns, file = paste0(modeldir.ns, "/m_nsrkc_predictions_iid_2024_12_02.RDS"))

# make a map function
plot_map_utm <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~SURVEY_YEAR) +
    coord_fixed()
}

# types of predictions
## predictions that incorporate all fixef and ranef
map.fe.re.utm <- plot_map_utm(predictions_noland4km_utm$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

map.fe.re.yrs <- plot_map(predictions$data %>% filter(SURVEY_YEAR > 2003), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

ggsave(file.path(plotdir.ns, "map_fixedeff_randomeff.png"), plot = map.fe.re.yrs, height = 4.2, width = 7, units = "in")

# just look at fixef
plot_map(predictions$data, exp(est_non_rf)) +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

# look at spatial ranef that are not accounted for by our fixef
## these represent consistent biotic and abiotic factors that are affecting biomass density but are not accounted for in the model
plot_map(predictions$data, omega_s) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# look at spatiotemporal ranefs that rep deviation from the fexef pred and the spatial ranef
## these represent biotic and abiotic factors that are changing through time and are not accounted for in the model
plot_map(predictions$data, epsilon_st) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()


# **************************************************************************************************************
# generate index ----
# **************************************************************************************************************

#use get_index() function to extract the total biomass calculations and standard errors
# not sure what argument to use for area. Area should be equal to area of grid cells that you're predicting over.
# check grid set up after converting to UTM
# when CPUE is catch/km^2
# can change aggregation function - if response variable is density, makes sense to sum. But could make more sense
# to take average over the area. 
# UTM is in meters

index.ar1.ns <- get_index(predictions.ar1.ns, area = predictions.ar1.ns$data$Area_km2, bias_correct = TRUE)
index.50kn.ns <- get_index(predictions.50kn.ns, area = predictions.50kn.ns$data$Area_km2, bias_correct = TRUE)
index.iid.ns <- get_index(predictions.iid.ns, area = predictions.iid.ns$data$Area_km2, bias_correct = TRUE)

# index values in metric tons rather than kg
index.ar1.ns.t <- index.ar1.ns %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.50kn.ns.t <- index.50kn.ns %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.ns.t <- index.iid.ns %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# save index values
write.csv(index.ar1.ns.t, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_ar1_2024_12_02.csv"))
write.csv(index.50kn.ns.t, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_50kn_2024_12_02.csv"))
write.csv(index.iid.ns.t, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_iid_2024_12_02.csv"))

# plot index
ggplot(index.ar1.ns.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.50kn.ns.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.iid.ns.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

# biomass estimates
mutate(index.ar1.ns, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))


# **************************************************************************************************************
# compare predicted index to observations ----
# **************************************************************************************************************


# ******************************************************************************