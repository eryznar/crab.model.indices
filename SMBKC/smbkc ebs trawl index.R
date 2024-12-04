# ************************************************************************************************************
# Generating a spatiotemporal model-based index for St Matthew Island blue king crab in the NMFS trawl survey
# September 2024
# Caitlin Stern
# ************************************************************************************************************

# **************************************************************************************************************
# load libraries, set plot preferences ----
# **************************************************************************************************************

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sp)
#library(sdmTMBextra)

cur_yr <- 2024

plotdir <- paste0(here::here(), "/SMBKC/plots")
modeldir <- paste0(here::here(), "/SMBKC/models")

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

grid4km <- readRDS(paste0(here::here(), "/SMBKC/data/SMBKC_Extrapolation_regions/SMBKC_Extrapolation_region/SMBKC_4kmgrid_No_Land.rds")) %>%
  # get UTM zone
  mutate(zone = (floor((Lon + 180)/6) %% 60) + 1) 

predgrid_utm_z1 <- grid4km %>% filter(zone == 1)
predgrid_utm_z2 <- grid4km %>% filter(zone == 2)

u1 <- predgrid_utm_z1
get_crs(u1, c("Lon", "Lat"))
u1u <- add_utm_columns(u1, c("Lon", "Lat"))

u2 <- predgrid_utm_z2
get_crs(u2, c("Lon", "Lat"))
u2u <- add_utm_columns(u2, c("Lon", "Lat"))

predgrid_utm <- rbind(u1u, u2u) %>%
  mutate(X = X / 1000, Y = Y / 1000) # convert UTM coordinates from meters to kilometers to reduce computing needs

# **************************************************************************************************************
# read in and process survey data ----
# **************************************************************************************************************

haul_bkc <- read.csv(paste0(here::here(), '/SMBKC/data/trawl_survey/EBSCrab_Haul_', cur_yr, '.csv'), skip = 5)

# error check
haul_bkc %>% filter(SAMPLING_FACTOR >= 1 & is.na(LENGTH) == TRUE)

# prepare data
bkc_kgkm <- haul_bkc %>% 
  filter(AKFIN_SURVEY_YEAR >= 1978, MID_LATITUDE > 58.5) %>% 
  dplyr::select(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, AREA_SWEPT, SPECIES_NAME, SEX, LENGTH, SAMPLING_FACTOR) %>% 
  # calculate weight in kg for males >= 90mm. Mean weight by stage: 0.7 kg for Stage-1, 1.2 kg for Stage-2, and 1.9 kg for Stage-3
  mutate(wt.kg1 = case_when(
    SEX == 1 & LENGTH >= 90 & LENGTH < 104 ~ 0.7,
    SEX == 1 & LENGTH >= 104 & LENGTH < 120 ~ 1.2,
    SEX == 1 & LENGTH >= 120 ~ 1.9,
    TRUE ~ 0
  )) %>%
  # make sampling factor column for males >= 90mm only; for all others sampling factor = 0
  mutate(sampling.factor = case_when(
    SEX == 1 & LENGTH >= 90 & LENGTH < 104 ~ SAMPLING_FACTOR,
    SEX == 1 & LENGTH >= 104 & LENGTH < 120 ~ SAMPLING_FACTOR,
    SEX == 1 & LENGTH >= 120 ~ SAMPLING_FACTOR,
    TRUE ~ 0
  )) %>%
  mutate(wt.kg = wt.kg1 * sampling.factor) %>%
  group_by(AKFIN_SURVEY_YEAR, GIS_STATION) %>% 
  mutate(total.wt.kg = sum(wt.kg)) %>% 
  mutate(AREA_SWEPT_km2 = AREA_SWEPT/0.29155335) %>% # convert from square nautical miles to square km
  mutate(kg.km = total.wt.kg / AREA_SWEPT_km2) %>%
  select(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, AREA_SWEPT, SPECIES_NAME, total.wt.kg, kg.km) %>%
  distinct() %>%
  arrange(AKFIN_SURVEY_YEAR, GIS_STATION) %>%
  rename("SURVEY_YEAR" = "AKFIN_SURVEY_YEAR", "LATITUDE" = "MID_LATITUDE", "LONGITUDE" = "MID_LONGITUDE") %>%
  mutate("year_f" = factor(SURVEY_YEAR))
  
bkc_kgkm_utm_zone <- bkc_kgkm %>%  
  # get UTM zone
  mutate(zone = (floor((LONGITUDE + 180)/6) %% 60) + 1) 

unique(bkc_kgkm_utm_zone$zone)

bkc_kgkm_utm_z1 <- bkc_kgkm_utm_zone %>% filter(zone == 1)
bkc_kgkm_utm_z2 <- bkc_kgkm_utm_zone %>% filter(zone == 2)
bkc_kgkm_utm_z3 <- bkc_kgkm_utm_zone %>% filter(zone == 3)

d1 <- bkc_kgkm_utm_z1
get_crs(d1, c("LONGITUDE", "LATITUDE"))
d1u <- add_utm_columns(d1, c("LONGITUDE", "LATITUDE"))

d2 <- bkc_kgkm_utm_z2
get_crs(d2, c("LONGITUDE", "LATITUDE"))
d2u <- add_utm_columns(d2, c("LONGITUDE", "LATITUDE"))

d3 <- bkc_kgkm_utm_z3
get_crs(d3, c("LONGITUDE", "LATITUDE"))
d3u <- add_utm_columns(d3, c("LONGITUDE", "LATITUDE"))

bkc_kgkm_utm <- rbind(d1u, d2u, d3u) %>%
  mutate(X = X / 1000, Y = Y / 1000) # convert UTM coordinates from meter to kilometers to reduce computing needs

# number of unique stations
length(unique(bkc_kgkm_utm$GIS_STATION))

# **************************************************************************************************************
# make SPDE mesh ----
# **************************************************************************************************************

#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), cutoff="10") 
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 100)
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 50)
BK_spde <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 120, type = "kmeans")
BK_spde_110kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 110, type = "kmeans")
BK_spde_100kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 100, type = "kmeans")
BK_spde_90kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 90, type = "kmeans")
BK_spde_75kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 75, type = "kmeans")
BK_spde_50kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 50, type = "kmeans")
BK_spde_25kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 25, type = "kmeans")
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 120)
plot(BK_spde)
BK_spde$mesh$n

#BK_spde.2 <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 10)

# look at spatial range of the model (Matern range - want 2 knots per range distance)


# **************************************************************************************************************
# fit and check models ----
# **************************************************************************************************************

# fit a GLMM with spatiotemporal = AR1 and 120 knots
m.smbkc.utm <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 120 knots
m.smbkc.iid <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = AR1 and 50 knots
m.smbkc.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 110 knots
m.smbkc.iid.110kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_110kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 100 knots
m.smbkc.iid.100kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_100kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 90 knots
m.smbkc.iid.90kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 75 knots
m.smbkc.iid.75kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_75kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 50 knots
m.smbkc.iid.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 25 knots
m.smbkc.iid.25kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_25kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# save the fitted models
saveRDS(m.smbkc.utm, file = paste0(modeldir, "/m_smbkc_2024_11_27.RDS"))
saveRDS(m.smbkc.iid, file = paste0(modeldir, "/m_smbkc_iid_2024_12_01.RDS"))
saveRDS(m.smbkc.50kn, file = paste0(modeldir, "/m_smbkc_50kn_2024_12_01.RDS"))
saveRDS(m.smbkc.iid.110kn, file = paste0(modeldir, "/m_smbkc_iid_110kn_2024_12_01.RDS"))
saveRDS(m.smbkc.iid.90kn, file = paste0(modeldir, "/m_smbkc_iid_90kn_2024_12_01.RDS"))
saveRDS(m.smbkc.iid.75kn, file = paste0(modeldir, "/m_smbkc_iid_75kn_2024_12_01.RDS"))
saveRDS(m.smbkc.iid.50kn, file = paste0(modeldir, "/m_smbkc_iid_50kn_2024_12_01.RDS"))
saveRDS(m.smbkc.iid.25kn, file = paste0(modeldir, "/m_smbkc_iid_25kn_2024_12_01.RDS"))
saveRDS(m.smbkc.iid.100kn, file = paste0(modeldir, "/m_smbkc_iid_100kn_2024_12_01.RDS"))

# print the model fit
m.smbkc.utm
m.smbkc.utm$sd_report

# view parameters
tidy.smbkc.utm <- tidy(m.smbkc.utm, conf.int = TRUE)
tidy.smbkc.ran.utm <- tidy(m.smbkc.utm , effects = "ran_pars", conf.int = TRUE)

# run sanity checks
sanity(m.smbkc.utm)
sanity(m.smbkc.iid)
sanity(m.smbkc.50kn)
sanity(m.smbkc.iid.110kn)
sanity(m.smbkc.iid.100kn)
sanity(m.smbkc.iid.90kn)
sanity(m.smbkc.iid.75kn)
sanity(m.smbkc.iid.50kn)
sanity(m.smbkc.iid.25kn)

# **************************************************************************************************************
# model residuals ----
# **************************************************************************************************************

# randomized quantile residuals
bkc_kgkm_utm$resids <- residuals(m.smbkc.utm) # randomized quantile residuals

hist(bkc_kgkm_utm$resids)
qqnorm(bkc_kgkm_utm$resids)
abline(a=0, b=1)

# DHARMa residuals
bkc_kgkm_utm$dharma_resids <- dharma_residuals(m.smbkc.utm)

fit_dl <- update(m.smbkc.utm, family = delta_lognormal())
simulate(fit_dl, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.utm)

ret <- simulate(fit_dl, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.utm, return_DHARMa = TRUE)
plot(ret)

m.smbkc.utm.resid <- simulate(update(m.smbkc.utm, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.utm, return_DHARMa = TRUE)
plot(m.smbkc.utm.resid)

m.smbkc.iid.resid <- simulate(update(m.smbkc.iid, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid, return_DHARMa = TRUE)
plot(m.smbkc.iid.resid)

m.smbkc.iid.100kn.resid <- simulate(update(m.smbkc.iid.100kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.100kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.100kn.resid)

m.smbkc.iid.90kn.resid <- simulate(update(m.smbkc.iid.90kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.90kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.90kn.resid)

m.smbkc.iid.75kn.resid <- simulate(update(m.smbkc.iid.75kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.75kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.75kn.resid)

m.smbkc.iid.50kn.resid <- simulate(update(m.smbkc.iid.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.50kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.50kn.resid)

m.smbkc.iid.25kn.resid <- simulate(update(m.smbkc.iid.25kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.25kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.25kn.resid)

m.smbkc.50kn.resid <- simulate(update(m.smbkc.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.50kn, return_DHARMa = TRUE)
plot(m.smbkc.50kn.resid)



ggplot(bkc_kgkm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()


ggplot(bkc_kgkm_utm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()


# **************************************************************************************************************
# generate predictions ----
# **************************************************************************************************************

# replicate grid across all years
pred_grid <- replicate_df(predgrid_utm, "year_f", unique(bkc_kgkm_utm$year_f))
pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
dplyr::glimpse(pred_grid)

pred_grid <- replicate_df(predgrid_utm, "SURVEY_YEAR", unique(bkc_kgkm_utm$SURVEY_YEAR))
pred_grid$year_f <- factor(pred_grid$SURVEY_YEAR)
dplyr::glimpse(pred_grid)

# predictions on new data
predictions <- predict(m.smbkc.utm, newdata = pred_grid, return_tmb_object = T)
predictions.50kn <- predict(m.smbkc.50kn, newdata = pred_grid, return_tmb_object = T)
predictions.iid <- predict(m.smbkc.iid, newdata = pred_grid, return_tmb_object = T)
predictions.iid.110kn <- predict(m.smbkc.iid.110kn, newdata = pred_grid, return_tmb_object = T)
predictions.iid.100kn <- predict(m.smbkc.iid.100kn, newdata = pred_grid, return_tmb_object = T)
predictions.iid.90kn <- predict(m.smbkc.iid.90kn, newdata = pred_grid, return_tmb_object = T)
predictions.iid.75kn <- predict(m.smbkc.iid.75kn, newdata = pred_grid, return_tmb_object = T)
predictions.iid.50kn <- predict(m.smbkc.iid.50kn, newdata = pred_grid, return_tmb_object = T)
predictions.iid.25kn <- predict(m.smbkc.iid.25kn, newdata = pred_grid, return_tmb_object = T)

# save the predictions
saveRDS(predictions, file = paste0(modeldir, "/m_smbkc_predictions_2024_11_27.RDS"))
saveRDS(predictions.50kn, file = paste0(modeldir, "/m_smbkc_predictions_50kn_2024_12_01.RDS"))
saveRDS(predictions.iid, file = paste0(modeldir, "/m_smbkc_predictions_iid_2024_12_01.RDS"))
saveRDS(predictions.iid.110kn, file = paste0(modeldir, "/m_smbkc_predictions_iid_110kn_2024_12_02.RDS"))
saveRDS(predictions.iid.100kn, file = paste0(modeldir, "/m_smbkc_predictions_iid_100kn_2024_12_02.RDS"))
saveRDS(predictions.iid.90kn, file = paste0(modeldir, "/m_smbkc_predictions_iid_90kn_2024_12_02.RDS"))
saveRDS(predictions.iid.75kn, file = paste0(modeldir, "/m_smbkc_predictions_iid_75kn_2024_12_02.RDS"))
saveRDS(predictions.iid.50kn, file = paste0(modeldir, "/m_smbkc_predictions_iid_50kn_2024_12_02.RDS"))
saveRDS(predictions.iid.25kn, file = paste0(modeldir, "/m_smbkc_predictions_iid_25kn_2024_12_02.RDS"))

# make a map function
plot_map_utm <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~SURVEY_YEAR) +
    coord_fixed()
}
 
# types of predictions
## predictions that incorporate all fixef and ranef
map.fe.re <- plot_map(predictions$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

map.fe.re.utm <- plot_map_utm(predictions_noland4km_utm$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

map.fe.re.yrs <- plot_map(predictions$data %>% filter(SURVEY_YEAR > 2003), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

ggsave(file.path(plotdir, "map_fixedeff_randomeff.png"), plot = map.fe.re.yrs, height = 4.2, width = 7, units = "in")

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

View(predictions)

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

index.ar1 <- get_index(predictions, area = predictions$data$Area_km2, bias_correct = TRUE)
index.50kn <- get_index(predictions.50kn, area = predictions.50kn$data$Area_km2, bias_correct = TRUE)
index.iid <- get_index(predictions.iid, area = predictions.iid$data$Area_km2, bias_correct = TRUE)
index.iid.110kn <- get_index(predictions.iid.110kn, area = predictions.iid.110kn$data$Area_km2, bias_correct = TRUE)
index.iid.100kn <- get_index(predictions.iid.100kn, area = predictions.iid.100kn$data$Area_km2, bias_correct = TRUE)
index.iid.90kn <- get_index(predictions.iid.90kn, area = predictions.iid.90kn$data$Area_km2, bias_correct = TRUE)
index.iid.75kn <- get_index(predictions.iid.75kn, area = predictions.iid.75kn$data$Area_km2, bias_correct = TRUE)
index.iid.50kn <- get_index(predictions.iid.50kn, area = predictions.iid.50kn$data$Area_km2, bias_correct = TRUE)
index.iid.25kn <- get_index(predictions.iid.25kn, area = predictions.iid.25kn$data$Area_km2, bias_correct = TRUE)

# index values in metric tons rather than kg
index.t <- index %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.50kn.t <- index.50kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.t <- index.iid %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.110kn.t <- index.iid.110kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.100kn.t <- index.iid.100kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.90kn.t <- index.iid.90kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.75kn.t <- index.iid.75kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.50kn.t <- index.iid.50kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.25kn.t <- index.iid.25kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# save index values
write.csv(index.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_2024_11_27.csv"))
write.csv(index.ar1.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_2024_12_02.csv"))
write.csv(index.50kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_50kn_2024_12_02.csv"))
write.csv(index.iid.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_2024_12_02.csv"))
write.csv(index.iid.110kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_110kn_2024_12_02.csv"))
write.csv(index.iid.100kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_100kn_2024_12_02.csv"))
write.csv(index.iid.90kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_90kn_2024_12_02.csv"))
write.csv(index.iid.75kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_75kn_2024_12_02.csv"))
write.csv(index.iid.50kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_50kn_2024_12_02.csv"))
write.csv(index.iid.25kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_25kn_2024_12_02.csv"))

# plot index
ggplot(index.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.50kn.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.iid.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

mod.compare <- rbind(index.t %>% mutate(model = "AR1, 120 knots"), index.50kn.t %>% mutate(model = "AR1, 50 knots")) %>%
  rbind(index.iid.t %>% mutate(model = "IID, 120 knots"))

mod.compare.iid <- rbind(index.iid.t %>% mutate(model = "120 knots"), index.iid.50kn.t %>% mutate(model = "50 knots")) %>%
  rbind(index.iid.110kn.t %>% mutate(model = "110 knots")) %>%
  #rbind(index.iid.100kn.t %>% mutate(model = "100 knots")) %>%
  rbind(index.iid.90kn.t %>% mutate(model = "90 knots")) %>%
  rbind(index.iid.75kn.t %>% mutate(model = "75 knots")) %>%
  rbind(index.iid.25kn.t %>% mutate(model = "25 knots"))

ggplot(mod.compare, aes(SURVEY_YEAR, est.t, group = model)) + 
  geom_line(aes(SURVEY_YEAR, est.t, color = model)) +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(mod.compare.iid, aes(SURVEY_YEAR, est.t, group = model)) + 
  geom_line(aes(SURVEY_YEAR, est.t, color = model)) +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

# biomass estimates
mutate(index, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))


# **************************************************************************************************************
# compare predicted index to observations ----
# **************************************************************************************************************

# plot predicted vs. observed index

obs.biom <- read.csv(paste0(here::here(), '/SMBKC/data/trawl_survey/survey_biomass_mt2.csv'))

obs.pred <- left_join(index.ar1.t, obs.biom, by = "SURVEY_YEAR") %>% filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
       obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

ggplot(obs.pred)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t), color = cbpalette[2])+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t), alpha = 0.2, fill = cbpalette[2]) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette)
  #scale_x_discrete(labels = yraxis$labels, breaks = yraxis$breaks)

# compare different models to each other and observations
obs.pred.comp <- left_join(mod.compare.iid, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

ggplot(obs.pred.comp)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette)

# plots with observations for models individually

ggplot(obs.pred.comp %>% filter(model == "AR1, 120 knots"))+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette)

ggplot(obs.pred.comp %>% filter(model == "AR1, 50 knots"))+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette)









ggplot(index.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

# ******************************************************************************
