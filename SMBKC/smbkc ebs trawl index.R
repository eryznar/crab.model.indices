# *************************************************************************************************
# Generating a spatiotemporal model-based index for St Matthew Island blue king crab in the NMFS trawl survey
# September 2024
# Caitlin Stern
# *************************************************************************************************

# *************************************************************************************************
# load libraries, set plot preferences ----
# *************************************************************************************************

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sp)
#library(sdmTMBextra)

cur_yr <- 2024

plotdir.sm <- paste0(here::here(), "/SMBKC/plots")
modeldir.sm <- paste0(here::here(), "/SMBKC/models")

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

# *************************************************************************************************
# read in grid used in VAST to use for predictions and convert to UTM ----
# *************************************************************************************************

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

# plot prediction grid
predgrid4kmplot <- ggplot(predgrid_utm, aes(x = Lon, y = Lat, color = Area_km2)) +
  geom_point() + 
  theme(legend.position="none") +
  labs(x = "Longitude", y = "Latitude")

ggsave(file.path(plotdir.sm, "prediction_grid_4km.png"), plot = predgrid4kmplot, height = 4.2, width = 7, units = "in")

# *************************************************************************************************
# read in and process survey data ----
# *************************************************************************************************

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

# *************************************************************************************************
# make SPDE mesh ----
# *************************************************************************************************

#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), cutoff="10") 
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 100)
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 50)
BK_spde_120kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 120, type = "kmeans")
BK_spde_110kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 110, type = "kmeans")
BK_spde_100kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 100, type = "kmeans")
BK_spde_90kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 90, type = "kmeans")
BK_spde_75kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 75, type = "kmeans")
BK_spde_50kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 50, type = "kmeans")
BK_spde_25kn <- make_mesh(bkc_kgkm_utm, xy_cols = c("X","Y"), n_knots = 25, type = "kmeans")
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 120)
plot(BK_spde_120kn)
BK_spde_120kn$mesh$n

#BK_spde.2 <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 10)

# look at spatial range of the model (Matern range - want 2 knots per range distance)

# plot mesh with data points ----

mesh.plot.120kn.sm <- ggplot(bkc_kgkm_utm) + 
  inlabru::gg(BK_spde_120kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

mesh.plot.90kn.sm <- ggplot(bkc_kgkm_utm) + 
  inlabru::gg(BK_spde_90kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

mesh.plot.50kn.sm <- ggplot(bkc_kgkm_utm) + 
  inlabru::gg(BK_spde_50kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

# export mesh plots

ggsave(file.path(plotdir.sm, "mesh_120kn_sm.png"), plot = mesh.plot.120kn.sm, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "mesh_90kn_sm.png"), plot = mesh.plot.90kn.sm, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "mesh_50kn_sm.png"), plot = mesh.plot.50kn.sm, height = 4.2, width = 7, units = "in")


# *************************************************************************************************
# fit and check models ----
# *************************************************************************************************

# fit IID/Tweedie models ----

# fit a GLMM with spatiotemporal = IID and 120 knots
m.smbkc.iid.tw.120kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 110 knots
m.smbkc.iid.tw.110kn <- sdmTMB(
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
m.smbkc.iid.tw.100kn <- sdmTMB(
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
m.smbkc.iid.tw.90kn <- sdmTMB(
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
m.smbkc.iid.tw.75kn <- sdmTMB(
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
m.smbkc.iid.tw.50kn <- sdmTMB(
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
m.smbkc.iid.tw.25kn <- sdmTMB(
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

# fit IID/delta gamma models ----

# IID, delta gamma, 120 knots
m.smbkc.iid.dg.120kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# IID, delta gamma, 90 knots
m.smbkc.iid.dg.90kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# IID, delta gamma, 50 knots
m.smbkc.iid.dg.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# fit IID/delta lognormal models ----

# IID, delta lognormal, 120 knots
m.smbkc.iid.dl.120kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# IID, delta lognormal, 90 knots
m.smbkc.iid.dl.90kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# IID, delta lognormal, 50 knots
m.smbkc.iid.dl.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# fit random walk/Tweedie models ----

# Tweedie, RW, 120 knots
m.smbkc.iid.rw.120kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# Tweedie, RW, 90 knots
m.smbkc.iid.rw.90kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# Tweedie, RW, 50 knots
m.smbkc.iid.rw.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit random walk/delta gamma models ----

# RW, delta gamma, 50 knots
m.smbkc.rw.dg.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# fit random walk/delta lognormal models ----

# RW, delta lognormal, 50 knots
m.smbkc.rw.dl.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# fit AR1/Tweedie models ----

# AR1, Tweedie, 120 knots
m.smbkc.ar1.tw.120kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# AR1, Tweedie, 90 knots
m.smbkc.ar1.tw.90kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# AR1, Tweedie, 50 knots
m.smbkc.ar1.tw.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit AR1/delta gamma models ----

# AR1, delta gamma, 50 knots
m.smbkc.ar1.dg.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  delta_gamma(type = "standard"))

# fit AR1/delta lognormal models ----

# AR1, delta lognormal, 50 knots
m.smbkc.ar1.dl.50kn <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  delta_lognormal(type = "standard"))



# save the fitted models ----

# IID/Tweedie models
saveRDS(m.smbkc.iid.tw.120kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_120kn.RDS"))
saveRDS(m.smbkc.iid.tw.110kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_110kn.RDS"))
saveRDS(m.smbkc.iid.tw.100kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_100kn.RDS"))
saveRDS(m.smbkc.iid.tw.90kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_90kn.RDS"))
saveRDS(m.smbkc.iid.tw.75kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_75kn.RDS"))
saveRDS(m.smbkc.iid.tw.50kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_50kn.RDS"))
saveRDS(m.smbkc.iid.tw.25kn, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_25kn.RDS"))

# IID/delta gamma models
saveRDS(m.smbkc.iid.dg.120kn, file = paste0(modeldir.sm, "/m_smbkc_iid_dg_120kn.RDS"))
saveRDS(m.smbkc.iid.dg.90kn, file = paste0(modeldir.sm, "/m_smbkc_iid_dg_90kn.RDS"))
saveRDS(m.smbkc.iid.dg.50kn, file = paste0(modeldir.sm, "/m_smbkc_iid_dg_50kn.RDS"))

# IID/delta lognormal models
saveRDS(m.smbkc.iid.dl.120kn, file = paste0(modeldir.sm, "/m_smbkc_iid_dl_120kn.RDS"))
saveRDS(m.smbkc.iid.dl.90kn, file = paste0(modeldir.sm, "/m_smbkc_iid_dl_90kn.RDS"))
saveRDS(m.smbkc.iid.dl.50kn, file = paste0(modeldir.sm, "/m_smbkc_iid_dl_50kn.RDS"))

# RW/Tweedie models
saveRDS(m.smbkc.iid.rw.120kn, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_120kn.RDS"))
saveRDS(m.smbkc.iid.rw.90kn, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_90kn.RDS"))
saveRDS(m.smbkc.iid.rw.50kn, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_50kn.RDS"))

# RW/dg models
saveRDS(m.smbkc.rw.dg.50kn, file = paste0(modeldir.sm, "/m_smbkc_rw_dg_50kn.RDS"))

# RW/dl models
saveRDS(m.smbkc.rw.dl.50kn, file = paste0(modeldir.sm, "/m_smbkc_rw_dl_50kn.RDS"))

# AR1/Tweedie models
saveRDS(m.smbkc.ar1.tw.120kn, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_120kn.RDS"))
saveRDS(m.smbkc.ar1.tw.90kn, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_90kn.RDS"))
saveRDS(m.smbkc.ar1.tw.50kn, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_50kn.RDS"))

# AR1/dg models
saveRDS(m.smbkc.ar1.dg.50kn, file = paste0(modeldir.sm, "/m_smbkc_ar1_dg_50kn.RDS"))

# AR1/dl models
saveRDS(m.smbkc.ar1.dl.50kn, file = paste0(modeldir.sm, "/m_smbkc_ar1_dl_50kn.RDS"))

# print the model fit
m.smbkc.utm
m.smbkc.utm$sd_report

# view parameters
tidy.smbkc.utm <- tidy(m.smbkc.utm, conf.int = TRUE)
tidy.smbkc.ran.utm <- tidy(m.smbkc.utm , effects = "ran_pars", conf.int = TRUE)

# run sanity checks ----

# IID/Tweedie models
sanity(m.smbkc.iid.tw.120kn)
sanity(m.smbkc.iid.tw.110kn)
sanity(m.smbkc.iid.tw.100kn)
sanity(m.smbkc.iid.tw.90kn)
sanity(m.smbkc.iid.tw.75kn)
sanity(m.smbkc.iid.tw.50kn)
sanity(m.smbkc.iid.tw.25kn)

# IID/delta gamma models
sanity(m.smbkc.iid.dg.120kn) # doesn't pass
sanity(m.smbkc.iid.dg.90kn) # doesn't pass
sanity(m.smbkc.iid.dg.50kn) # doesn't pass

# IID/delta lognormal models
sanity(m.smbkc.iid.dl.120kn) # doesn't pass
sanity(m.smbkc.iid.dl.90kn) # doesn't pass
sanity(m.smbkc.iid.dl.50kn) # model won't fit - error

# RW/Tweedie models
sanity(m.smbkc.iid.rw.120kn) 
sanity(m.smbkc.iid.rw.90kn) 
sanity(m.smbkc.iid.rw.50kn) 

# RW/dg models
sanity(m.smbkc.rw.dg.50kn) # doesn't pass

# RW/dl models
sanity(m.smbkc.rw.dl.50kn) # doesn't pass

# AR1/Tweedie models
sanity(m.smbkc.ar1.tw.120kn) 
sanity(m.smbkc.ar1.tw.90kn) 
sanity(m.smbkc.ar1.tw.50kn) 

# AR1/dg models
sanity(m.smbkc.ar1.dg.50kn) # doesn't pass

# AR1/dl models
sanity(m.smbkc.ar1.dl.50kn) # doesn't pass

# ************************************************************************************************
# compare models using cross-validation ----
# ************************************************************************************************

#clust <- sample(seq_len(2), size = nrow(nsrkc_utm), replace = TRUE)
#clust <- kmeans(nsrkc_utm[, c("X", "Y")], 20)$cluster
#clust <- sample(1:2, size = nrow(nsrkc_utm), replace = TRUE, prob = c(0.1, 0.9))

# cross validation for IID/Tweedie models ----

#clust.iid.tw.120kn <- sample(seq_len(10), size = nrow(m.smbkc.iid.tw.120kn$data), replace = TRUE)

m.smbkc.iid.tw.120kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.iid.tw.120kn)

#clust.iid.tw.90kn <- sample(seq_len(10), size = nrow(m.smbkc.iid.tw.90kn$data), replace = TRUE)

m.smbkc.iid.tw.90kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.iid.tw.90kn)

#clust.iid.tw.50kn <- sample(seq_len(10), size = nrow(m.smbkc.iid.tw.50kn$data), replace = TRUE)

m.smbkc.iid.tw.50kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "iid",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.iid.tw.50kn)


# cross validation for RW/Tweedie models ----

#clust.rw.tw.120kn <- sample(seq_len(10), size = nrow(m.smbkc.iid.rw.120kn$data), replace = TRUE)

m.smbkc.rw.tw.120kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.rw.tw.120kn)

#clust.rw.tw.90kn <- sample(seq_len(10), size = nrow(m.smbkc.iid.rw.90kn$data), replace = TRUE)

m.smbkc.rw.tw.90kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.rw.tw.90kn)

#clust.rw.tw.50kn <- sample(seq_len(10), size = nrow(m.smbkc.iid.rw.50kn$data), replace = TRUE)

m.smbkc.rw.tw.50kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "rw",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.rw.tw.50kn)

# cross validation for AR1/Tweedie models ----

#clust.ar1.tw.120kn <- sample(seq_len(10), size = nrow(m.smbkc.ar1.tw.120kn$data), replace = TRUE)

m.smbkc.ar1.tw.120kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_120kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.ar1.tw.120kn)

#clust.ar1.tw.90kn <- sample(seq_len(10), size = nrow(m.smbkc.ar1.tw.90kn$data), replace = TRUE)

m.smbkc.ar1.tw.90kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f, 
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_90kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.ar1.tw.90kn)

#clust.ar1.tw.50kn <- sample(seq_len(10), size = nrow(m.smbkc.ar1.tw.50kn$data), replace = TRUE)

m.smbkc.ar1.tw.50kn.cv <- sdmTMB_cv(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + year_f,
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde_50kn, 
  spatiotemporal = "ar1",
  extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)
  #fold_ids = clust.ar1.tw.50kn)


# save cross-validation models ----

# IID/Tweedie cross-validation models
saveRDS(m.smbkc.iid.tw.120kn.cv, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_120kn_cv.RDS"))
saveRDS(m.smbkc.iid.tw.90kn.cv, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_90kn_cv.RDS"))
saveRDS(m.smbkc.iid.tw.50kn.cv, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_50kn_cv.RDS"))

# RW/Tweedie cross-validation models
saveRDS(m.smbkc.rw.tw.120kn.cv, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_120kn_cv.RDS"))
saveRDS(m.smbkc.rw.tw.90kn.cv, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_90kn_cv.RDS"))
saveRDS(m.smbkc.rw.tw.50kn.cv, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_50kn_cv.RDS"))

# AR1/Tweedie cross-validation models
saveRDS(m.smbkc.ar1.tw.120kn.cv, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_120kn_cv.RDS"))
saveRDS(m.smbkc.ar1.tw.90kn.cv, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_90kn_cv.RDS"))
saveRDS(m.smbkc.ar1.tw.50kn.cv, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_50kn_cv.RDS"))

# table of loglik values

# make table with loglikelihood results ----
m.smbkc.compare.loglik <- data.frame(family=character(),
                                     estimation=character(), 
                                     knots=numeric(), 
                                     loglik=numeric(),
                                     stringsAsFactors=FALSE) %>%
  # IID models
  add_row(family = "Tweedie", estimation = "IID", knots = 120, loglik = m.smbkc.iid.120kn.tw.cv$sum_loglik) %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 90, loglik = m.smbkc.iid.90kn.tw.cv$sum_loglik) %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 50, loglik = m.smbkc.iid.50kn.tw.cv$sum_loglik) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 120, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 90, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 50, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 120, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 90, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 50, loglik = NA) %>%
  # RW models
  add_row(family = "Tweedie", estimation = "RW", knots = 120, loglik = m.smbkc.rw.tw.120kn.cv$sum_loglik) %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 90, loglik = m.smbkc.rw.tw.90kn.cv$sum_loglik) %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 50, loglik = m.smbkc.rw.tw.50kn.cv$sum_loglik) %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 50, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 50, loglik = NA) %>%
# AR1 models
  add_row(family = "Tweedie", estimation = "AR1", knots = 120, loglik = m.smbkc.ar1.tw.120kn.cv$sum_loglik) %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 90, loglik = m.smbkc.ar1.tw.90kn.cv$sum_loglik) %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 50, loglik = m.smbkc.ar1.tw.50kn.cv$sum_loglik) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 50, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 50, loglik = NA)
  
# save table with predictive skill values
write.csv(m.smbkc.compare.loglik, paste0(here::here(), "/SMBKC/output/smbkc_compare_loglik.csv"))

# ************************************************************************************************
# model residuals ----
# *************************************************************************************************

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

# IID/Tweedie DHARMa residuals ----

m.smbkc.iid.tw.120kn.resid <- simulate(update(m.smbkc.iid.tw.120kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.tw.120kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.tw.120kn.resid)

m.smbkc.iid.tw.90kn.resid <- simulate(update(m.smbkc.iid.tw.90kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.tw.90kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.tw.90kn.resid)

m.smbkc.iid.tw.50kn.resid <- simulate(update(m.smbkc.iid.tw.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.tw.50kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.tw.50kn.resid)

# IID/delta gamma DHARMa residuals ----

m.smbkc.iid.dg.120kn.resid <- simulate(update(m.smbkc.iid.dg.120kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.dg.120kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.dg.120kn.resid)

m.smbkc.iid.dg.90kn.resid <- simulate(update(m.smbkc.iid.dg.90kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.dg.90kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.dg.90kn.resid)

m.smbkc.iid.dg.50kn.resid <- simulate(update(m.smbkc.iid.dg.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.dg.50kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.dg.50kn.resid)

# IID/delta lognormal DHARMa residuals ----

m.smbkc.iid.dl.120kn.resid <- simulate(update(m.smbkc.iid.dl.120kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.dl.120kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.dl.120kn.resid)

m.smbkc.iid.dl.90kn.resid <- simulate(update(m.smbkc.iid.dl.90kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.dl.90kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.dl.90kn.resid)

m.smbkc.iid.dl.50kn.resid <- simulate(update(m.smbkc.iid.dl.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.dl.50kn, return_DHARMa = TRUE)
plot(m.smbkc.iid.dl.50kn.resid)

# RW/Tweedie DHARMa residuals ----

m.smbkc.rw.tw.120kn.resid <- simulate(update(m.smbkc.iid.rw.120kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.rw.120kn, return_DHARMa = TRUE)
plot(m.smbkc.rw.tw.120kn.resid)

m.smbkc.rw.tw.90kn.resid <- simulate(update(m.smbkc.iid.rw.90kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.rw.90kn, return_DHARMa = TRUE)
plot(m.smbkc.rw.tw.90kn.resid)

m.smbkc.rw.tw.50kn.resid1 <- simulate(update(m.smbkc.iid.rw.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  #dharma_residuals(m.smbkc.iid.rw.50kn, return_DHARMa = TRUE)
  dharma_residuals(m.smbkc.iid.rw.50kn, plot = FALSE)

m.smbkc.rw.tw.50kn.resid2 <- simulate(update(m.smbkc.iid.rw.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.iid.rw.50kn, return_DHARMa = TRUE)
plot(m.smbkc.rw.tw.50kn.resid)

ggplot()+
  theme_bw()+
  geom_point(m.smbkc.rw.tw.50kn.resid1, mapping = aes(expected, observed), size = 2, fill = "black")+
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("Expected")+
  xlab("Observed")+
  ggtitle("SMBKC random walk, Tweedie, 50 knots")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave()


# add all of the scaled residuals to one data frame and then make spatial residual plots ----

bkc_kgkm_utm_resids <- cbind(bkc_kgkm_utm, 
                             data.frame(m.smbkc.iid.tw.120kn.resid$scaledResiduals),
                             data.frame(m.smbkc.iid.tw.90kn.resid$scaledResiduals),
                             data.frame(m.smbkc.rw.tw.120kn.resid$scaledResiduals),
                             data.frame(m.smbkc.rw.tw.90kn.resid$scaledResiduals),
                             data.frame(m.smbkc.rw.tw.50kn.resid$scaledResiduals),
                             data.frame(m.smbkc.ar1.tw.120kn.resid$scaledResiduals),
                             data.frame(m.smbkc.ar1.tw.90kn.resid$scaledResiduals),
                             data.frame(m.smbkc.ar1.tw.50kn.resid$scaledResiduals))


spat.resid.iid.tw.120kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.iid.tw.120kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC IID Tweedie model residuals (120 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.iid.tw.90kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.iid.tw.90kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC IID Tweedie model residuals (90 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.rw.tw.120kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.rw.tw.120kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC random walk Tweedie model residuals (120 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.rw.tw.90kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.rw.tw.90kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC random walk Tweedie model residuals (90 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.rw.tw.50kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.rw.tw.50kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC random walk Tweedie model residuals (50 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.ar1.tw.120kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.ar1.tw.120kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC AR1 Tweedie model residuals (120 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.ar1.tw.90kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.ar1.tw.90kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC AR1 Tweedie model residuals (90 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

spat.resid.ar1.tw.50kn <- ggplot(bkc_kgkm_utm_resids) + 
  geom_point(aes(y = LATITUDE, x = LONGITUDE, color = m.smbkc.ar1.tw.50kn.resid.scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("SMBKC AR1 Tweedie model residuals (50 knots)")+
  facet_wrap(~SURVEY_YEAR)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

ggsave(paste0(plotdir.sm, "/spat_resid_iid_tw_120kn.png"), spat.resid.iid.tw.120kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_iid_tw_90kn.png"), spat.resid.iid.tw.90kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_rw_tw_120kn.png"), spat.resid.rw.tw.120kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_rw_tw_90kn.png"), spat.resid.rw.tw.90kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_rw_tw_50kn.png"), spat.resid.rw.tw.50kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_ar1_tw_120kn.png"), spat.resid.ar1.tw.120kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_ar1_tw_90kn.png"), spat.resid.ar1.tw.90kn, height = 7, width = 7, units = "in")
ggsave(paste0(plotdir.sm, "/spat_resid_ar1_tw_50kn.png"), spat.resid.ar1.tw.50kn, height = 7, width = 7, units = "in")

# AR1/Tweedie DHARMa residuals ----

m.smbkc.ar1.tw.120kn.resid <- simulate(update(m.smbkc.ar1.tw.120kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.ar1.tw.120kn, return_DHARMa = TRUE)
plot(m.smbkc.ar1.tw.120kn.resid)

m.smbkc.ar1.tw.90kn.resid <- simulate(update(m.smbkc.ar1.tw.90kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.ar1.tw.90kn, return_DHARMa = TRUE)
plot(m.smbkc.ar1.tw.90kn.resid)

m.smbkc.ar1.tw.50kn.resid <- simulate(update(m.smbkc.ar1.tw.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.ar1.tw.50kn, return_DHARMa = TRUE)
plot(m.smbkc.ar1.tw.50kn.resid)

# save DHARMa residuals
saveRDS(m.smbkc.iid.tw.120kn.resid, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_120kn_resid.RDS"))
saveRDS(m.smbkc.iid.tw.90kn.resid, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_90kn_resid.RDS"))
saveRDS(m.smbkc.iid.tw.50kn.resid, file = paste0(modeldir.sm, "/m_smbkc_iid_tw_50kn_resid.RDS"))
saveRDS(m.smbkc.rw.tw.120kn.resid, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_120kn_resid.RDS"))
saveRDS(m.smbkc.rw.tw.90kn.resid, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_90kn_resid.RDS"))
saveRDS(m.smbkc.rw.tw.50kn.resid, file = paste0(modeldir.sm, "/m_smbkc_rw_tw_50kn_resid.RDS"))
saveRDS(m.smbkc.ar1.tw.120kn.resid, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_120kn_resid.RDS"))
saveRDS(m.smbkc.ar1.tw.90kn.resid, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_90kn_resid.RDS"))
saveRDS(m.smbkc.ar1.tw.50kn.resid, file = paste0(modeldir.sm, "/m_smbkc_ar1_tw_50kn_resid.RDS"))

# bring in saved residuals
m.smbkc.iid.tw.120kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_iid_tw_120kn_resid.RDS"))
m.smbkc.iid.tw.90kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_iid_tw_90kn_resid.RDS"))
m.smbkc.iid.tw.50kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_iid_tw_50kn_resid.RDS"))
m.smbkc.rw.tw.120kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_rw_tw_120kn_resid.RDS"))
m.smbkc.rw.tw.90kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_rw_tw_90kn_resid.RDS"))
m.smbkc.rw.tw.50kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_rw_tw_50kn_resid.RDS"))
m.smbkc.ar1.tw.120kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_ar1_tw_120kn_resid.RDS"))
m.smbkc.ar1.tw.90kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_ar1_tw_90kn_resid.RDS"))
m.smbkc.ar1.tw.50kn.resid <- readRDS(paste0(modeldir.sm, "/m_smbkc_ar1_tw_50kn_resid.RDS"))

# plot residuals

ggplot(bkc_kgkm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()


ggplot(bkc_kgkm_utm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()

# DHARMa residuals tests ----

# IID/Tweedie models
DHARMa::testOutliers(m.smbkc.iid.tw.120kn.resid)
DHARMa::testOutliers(m.smbkc.iid.tw.90kn.resid)
DHARMa::testOutliers(m.smbkc.iid.tw.50kn.resid)

DHARMa::testQuantiles(m.smbkc.iid.tw.120kn.resid)
DHARMa::testQuantiles(m.smbkc.iid.tw.90kn.resid)
DHARMa::testQuantiles(m.smbkc.iid.tw.50kn.resid)

DHARMa::testDispersion(m.smbkc.iid.tw.120kn.resid)
DHARMa::testDispersion(m.smbkc.iid.tw.90kn.resid)
DHARMa::testDispersion(m.smbkc.iid.tw.50kn.resid)

DHARMa::testResiduals(m.smbkc.iid.tw.120kn.resid)
DHARMa::testResiduals(m.smbkc.iid.tw.90kn.resid)
DHARMa::testResiduals(m.smbkc.iid.tw.50kn.resid)

DHARMa::testSpatialAutocorrelation(m.smbkc.iid.120kn.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)
DHARMa::testSpatialAutocorrelation(m.smbkc.iid.90kn.dgs.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)

DHARMa::testZeroInflation(m.smbkc.iid.tw.120kn.resid)
DHARMa::testZeroInflation(m.smbkc.iid.tw.90kn.resid)
DHARMa::testZeroInflation(m.smbkc.iid.tw.50kn.resid)

# RW/Tweedie models
DHARMa::testOutliers(m.smbkc.rw.tw.120kn.resid)
DHARMa::testOutliers(m.smbkc.rw.tw.90kn.resid)
DHARMa::testOutliers(m.smbkc.rw.tw.50kn.resid)

DHARMa::testQuantiles(m.smbkc.rw.tw.120kn.resid)
DHARMa::testQuantiles(m.smbkc.rw.tw.90kn.resid)
DHARMa::testQuantiles(m.smbkc.rw.tw.50kn.resid)

DHARMa::testDispersion(m.smbkc.rw.tw.120kn.resid)
DHARMa::testDispersion(m.smbkc.rw.tw.90kn.resid)
DHARMa::testDispersion(m.smbkc.rw.tw.50kn.resid)

DHARMa::testResiduals(m.smbkc.rw.tw.120kn.resid)
DHARMa::testResiduals(m.smbkc.rw.tw.90kn.resid)
DHARMa::testResiduals(m.smbkc.rw.tw.50kn.resid)

DHARMa::testSpatialAutocorrelation(m.smbkc.rw.120kn.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)
DHARMa::testSpatialAutocorrelation(m.smbkc.rw.90kn.dgs.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)

DHARMa::testZeroInflation(m.smbkc.rw.tw.120kn.resid)
DHARMa::testZeroInflation(m.smbkc.rw.tw.90kn.resid)
DHARMa::testZeroInflation(m.smbkc.rw.tw.50kn.resid)

# AR1/Tweedie models

DHARMa::testOutliers(m.smbkc.ar1.tw.120kn.resid)
DHARMa::testOutliers(m.smbkc.ar1.tw.90kn.resid)
DHARMa::testOutliers(m.smbkc.ar1.tw.50kn.resid)

DHARMa::testQuantiles(m.smbkc.ar1.tw.120kn.resid)
DHARMa::testQuantiles(m.smbkc.ar1.tw.90kn.resid)
DHARMa::testQuantiles(m.smbkc.ar1.tw.50kn.resid)

DHARMa::testDispersion(m.smbkc.ar1.tw.120kn.resid)
DHARMa::testDispersion(m.smbkc.ar1.tw.90kn.resid)
DHARMa::testDispersion(m.smbkc.ar1.tw.50kn.resid)

DHARMa::testResiduals(m.smbkc.ar1.tw.120kn.resid)
DHARMa::testResiduals(m.smbkc.ar1.tw.90kn.resid)
DHARMa::testResiduals(m.smbkc.ar1.tw.50kn.resid)

DHARMa::testSpatialAutocorrelation(m.smbkc.ar1.tw.120kn.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)
DHARMa::testSpatialAutocorrelation(m.smbkc.ar1.tw.90kn.dgs.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)

DHARMa::testZeroInflation(m.smbkc.ar1.tw.120kn.resid)
DHARMa::testZeroInflation(m.smbkc.ar1.tw.90kn.resid)
DHARMa::testZeroInflation(m.smbkc.ar1.tw.50kn.resid)


# make table with residual analysis results ----
m.smbkc.compare.resid <- data.frame(family=character(),
                                     estimation=character(), 
                                     knots=numeric(), 
                                     quantile=character(),
                                     dispersion=character(),
                                     outliers=character(),
                                     zero_inf=character(),
                                     stringsAsFactors=FALSE) %>%
  # IID models
  add_row(family = "Tweedie", estimation = "IID", knots = 120, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "sig.") %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 90, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 120, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 90, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 120, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 90, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  # RW models
  add_row(family = "Tweedie", estimation = "RW", knots = 120, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "sig.") %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 90, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 50, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  # AR1 models
  add_row(family = "Tweedie", estimation = "AR1", knots = 120, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 90, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 50, quantile = "sig.", dispersion = "under", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA)

# save table with residuals results
write.csv(m.smbkc.compare.resid, paste0(here::here(), "/SMBKC/output/smbkc_compare_resid.csv"))

# ************************************************************************************************
# generate predictions ----
# *************************************************************************************************

# replicate grid across all years
pred_grid <- replicate_df(predgrid_utm, "year_f", unique(bkc_kgkm_utm$year_f))
pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
dplyr::glimpse(pred_grid)

predgrid_utm <- replicate_df(predgrid_utm, "SURVEY_YEAR", unique(bkc_kgkm_utm$SURVEY_YEAR))
predgrid_utm$year_f <- factor(pred_grid$SURVEY_YEAR)
dplyr::glimpse(predgrid_utm)

# predictions on new data ----

# IID/Tweedie predictions ----
pred.iid.tw.120kn <- predict(m.smbkc.iid.tw.120kn, newdata = predgrid_utm, return_tmb_object = T)
pred.iid.tw.110kn <- predict(m.smbkc.iid.tw.110kn, newdata = predgrid_utm, return_tmb_object = T)
pred.iid.tw.100kn <- predict(m.smbkc.iid.tw.100kn, newdata = predgrid_utm, return_tmb_object = T)
pred.iid.tw.90kn <- predict(m.smbkc.iid.tw.90kn, newdata = predgrid_utm, return_tmb_object = T)
pred.iid.tw.75kn <- predict(m.smbkc.iid.tw.75kn, newdata = predgrid_utm, return_tmb_object = T)
pred.iid.tw.50kn <- predict(m.smbkc.iid.tw.50kn, newdata = predgrid_utm, return_tmb_object = T)
pred.iid.tw.25kn <- predict(m.smbkc.iid.tw.25kn, newdata = predgrid_utm, return_tmb_object = T)

# IID/delta gamma predictions ----
pred.iid.dg.120kn <- predict(m.smbkc.iid.dg.120kn, newdata = pred_grid, return_tmb_object = T)
pred.iid.dg.90kn <- predict(m.smbkc.iid.dg.90kn, newdata = pred_grid, return_tmb_object = T)
pred.iid.dg.50kn <- predict(m.smbkc.iid.dg.50kn, newdata = pred_grid, return_tmb_object = T)

# IID/delta lognormal predictions ----
pred.iid.dl.120kn <- predict(m.smbkc.iid.dl.120kn, newdata = pred_grid, return_tmb_object = T)
pred.iid.dl.90kn <- predict(m.smbkc.iid.dl.90kn, newdata = pred_grid, return_tmb_object = T)
pred.iid.dl.50kn <- predict(m.smbkc.iid.dl.50kn, newdata = pred_grid, return_tmb_object = T)

# RW/Tweedie predictions
pred.rw.tw.120kn <- predict(m.smbkc.iid.rw.120kn, newdata = pred_grid, return_tmb_object = T)
pred.rw.tw.90kn <- predict(m.smbkc.iid.rw.90kn, newdata = pred_grid, return_tmb_object = T)
pred.rw.tw.50kn <- predict(m.smbkc.iid.rw.50kn, newdata = pred_grid, return_tmb_object = T)

# AR1/Tweedie predictions
pred.ar1.tw.120kn <- predict(m.smbkc.ar1.tw.120kn, newdata = pred_grid, return_tmb_object = T)
pred.ar1.tw.90kn <- predict(m.smbkc.ar1.tw.90kn, newdata = pred_grid, return_tmb_object = T)
pred.ar1.tw.50kn <- predict(m.smbkc.ar1.tw.50kn, newdata = pred_grid, return_tmb_object = T)

# save the predictions
saveRDS(pred.iid.tw.120kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_120kn.RDS"))
saveRDS(pred.iid.tw.110kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_110kn.RDS"))
saveRDS(pred.iid.tw.100kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_100kn.RDS"))
saveRDS(pred.iid.tw.90kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_90kn.RDS"))
saveRDS(pred.iid.tw.75kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_75kn.RDS"))
saveRDS(pred.iid.tw.50kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_50kn.RDS"))
saveRDS(pred.iid.tw.25kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_tw_25kn.RDS"))

saveRDS(pred.iid.dg.120kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_dg_120kn.RDS"))
saveRDS(pred.iid.dg.90kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_dg_90kn.RDS"))
saveRDS(pred.iid.dg.50kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_dg_50kn.RDS"))

saveRDS(pred.iid.dl.120kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_dl_120kn.RDS"))
saveRDS(pred.iid.dl.90kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_dl_90kn.RDS"))
saveRDS(pred.iid.dl.50kn, file = paste0(modeldir.sm, "/m_smbkc_pred_iid_dl_50kn.RDS"))

saveRDS(pred.rw.tw.120kn, file = paste0(modeldir.sm, "/m_smbkc_pred_rw_tw_120kn.RDS"))
saveRDS(pred.rw.tw.90kn, file = paste0(modeldir.sm, "/m_smbkc_pred_rw_tw_90kn.RDS"))
saveRDS(pred.rw.tw.50kn, file = paste0(modeldir.sm, "/m_smbkc_pred_rw_tw_50kn.RDS"))

saveRDS(pred.ar1.tw.120kn, file = paste0(modeldir.sm, "/m_smbkc_pred_ar1_tw_120kn.RDS"))
saveRDS(pred.ar1.tw.90kn, file = paste0(modeldir.sm, "/m_smbkc_pred_ar1_tw_90kn.RDS"))
saveRDS(pred.ar1.tw.50kn, file = paste0(modeldir.sm, "/m_smbkc_pred_ar1_tw_50kn.RDS"))

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

ggsave(file.path(plotdir.sm, "map_fixedeff_randomeff.png"), plot = map.fe.re.yrs, height = 4.2, width = 7, units = "in")

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

# bring in saved predictions
pred.iid.tw.120kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_iid_tw_120kn.rds"))
pred.iid.tw.90kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_iid_tw_90kn.rds"))
pred.iid.tw.50kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_iid_tw_50kn.rds"))
pred.rw.tw.120kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_rw_tw_120kn.rds"))
pred.rw.tw.90kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_rw_tw_90kn.rds"))
pred.rw.tw.50kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_rw_tw_50kn.rds"))
pred.ar1.tw.120kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_ar1_tw_120kn.rds"))
pred.ar1.tw.90kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_ar1_tw_90kn.rds"))
pred.ar1.tw.50kn <- readRDS(paste0(here::here(), "/SMBKC/models/m_smbkc_pred_ar1_tw_50kn.rds"))


# spatial biomass prediction plots ----

midpoint <- pred.iid.tw.50kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.iid.tw.50kn <- ggplot(pred.iid.tw.50kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC IID, Tweedie, 50 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.iid.tw.90kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.iid.tw.90kn <- ggplot(pred.iid.tw.90kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC IID, Tweedie, 90 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.iid.tw.120kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.iid.tw.120kn <- ggplot(pred.iid.tw.120kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC IID, Tweedie, 120 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.rw.tw.50kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.rw.tw.50kn <- ggplot(pred.rw.tw.50kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC RW, Tweedie, 50 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.rw.tw.90kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.rw.tw.90kn <- ggplot(pred.rw.tw.90kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC RW, Tweedie, 90 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.rw.tw.120kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.rw.tw.120kn <- ggplot(pred.rw.tw.120kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC RW, Tweedie, 120 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.ar1.tw.50kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.ar1.tw.50kn <- ggplot(pred.ar1.tw.50kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC AR1, Tweedie, 50 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.ar1.tw.90kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.ar1.tw.90kn <- ggplot(pred.ar1.tw.90kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC AR1, Tweedie, 90 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

midpoint <- pred.ar1.tw.120kn$data %>% filter(SURVEY_YEAR > 2011) %>% mutate(est.exp = exp(est)/1000) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.ar1.tw.120kn <- ggplot(pred.ar1.tw.120kn$data %>% filter(SURVEY_YEAR > 2011)) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)/1000), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted biomass") +
  theme_gray() + 
  facet_wrap(~SURVEY_YEAR)+
  ggtitle("SMBKC AR1, Tweedie, 120 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

# save spatial biomass prediction plots

ggsave(file.path(plotdir.sm, "spatial_pred_iid_tw_120kn.png"), plot = sppred.iid.tw.120kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_iid_tw_90kn.png"), plot = sppred.iid.tw.90kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_iid_tw_50kn.png"), plot = sppred.iid.tw.50kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_rw_tw_120kn.png"), plot = sppred.rw.tw.120kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_rw_tw_90kn.png"), plot = sppred.rw.tw.90kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_rw_tw_50kn.png"), plot = sppred.rw.tw.50kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_ar1_tw_120kn.png"), plot = sppred.ar1.tw.120kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_ar1_tw_90kn.png"), plot = sppred.ar1.tw.90kn, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.sm, "spatial_pred_ar1_tw_50kn.png"), plot = sppred.ar1.tw.50kn, height = 5, width = 7, units = "in")

# *************************************************************************************************
# generate index ----
# *************************************************************************************************

#use get_index() function to extract the total biomass calculations and standard errors
# not sure what argument to use for area. Area should be equal to area of grid cells that you're predicting over.
# check grid set up after converting to UTM
# when CPUE is catch/km^2
# can change aggregation function - if response variable is density, makes sense to sum. But could make more sense
# to take average over the area. 
# UTM is in meters

# IID/Tweedie indices
index.iid.tw.120kn <- get_index(pred.iid.tw.120kn, area = pred.iid.tw.120kn$data$Area_km2, bias_correct = TRUE)
index.iid.110kn <- get_index(pred.iid.110kn, area = pred.iid.110kn$data$Area_km2, bias_correct = TRUE)
index.iid.tw.100kn <- get_index(pred.iid.100kn, area = pred.iid.100kn$data$Area_km2, bias_correct = TRUE)
index.iid.tw.90kn <- get_index(pred.iid.tw.90kn, area = pred.iid.tw.90kn$data$Area_km2, bias_correct = TRUE)
index.iid.75kn <- get_index(pred.iid.75kn, area = pred.iid.75kn$data$Area_km2, bias_correct = TRUE)
index.iid.tw.50kn <- get_index(pred.iid.tw.50kn, area = pred.iid.tw.50kn$data$Area_km2, bias_correct = TRUE)
index.iid.25kn <- get_index(pred.iid.25kn, area = pred.iid.25kn$data$Area_km2, bias_correct = TRUE)

# RW/Tweedie indices
index.rw.tw.120kn <- get_index(pred.rw.tw.120kn, area = pred.rw.tw.120kn$data$Area_km2, bias_correct = TRUE)
index.rw.tw.90kn <- get_index(pred.rw.tw.90kn, area = pred.rw.tw.90kn$data$Area_km2, bias_correct = TRUE)
index.rw.tw.50kn <- get_index(pred.rw.tw.50kn, area = pred.rw.tw.50kn$data$Area_km2, bias_correct = TRUE)

# AR1/Tweedie indices
index.ar1.tw.120kn <- get_index(pred.ar1.tw.120kn, area = pred.ar1.tw.120kn$data$Area_km2, bias_correct = TRUE)
index.ar1.tw.90kn <- get_index(pred.ar1.tw.90kn, area = pred.ar1.tw.90kn$data$Area_km2, bias_correct = TRUE)
index.ar1.tw.50kn <- get_index(pred.ar1.tw.50kn, area = pred.ar1.tw.50kn$data$Area_km2, bias_correct = TRUE)


# IID/Tweedie index values in metric tons rather than kg
index.t <- index %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.50kn.t <- index.50kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.tw.120kn.t <- index.iid.tw.120kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.110kn.t <- index.iid.110kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.100kn.t <- index.iid.100kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.tw.90kn.t <- index.iid.tw.90kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.75kn.t <- index.iid.75kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.tw.50kn.t <- index.iid.tw.50kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.iid.25kn.t <- index.iid.25kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# RW/Tweedie index values in metric tons rather than kg
index.rw.tw.120kn.t <- index.rw.tw.120kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.rw.tw.90kn.t <- index.rw.tw.90kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.rw.tw.50kn.t <- index.rw.tw.50kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# AR1/Tweedie index values in metric tons rather than kg
index.ar1.tw.120kn.t <- index.ar1.tw.120kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.ar1.tw.90kn.t <- index.ar1.tw.90kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)
index.ar1.tw.50kn.t <- index.ar1.tw.50kn %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# save IID/Tweedie index values
write.csv(index.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_2024_11_27.csv"))
write.csv(index.ar1.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_2024_12_02.csv"))
write.csv(index.50kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_50kn_2024_12_02.csv"))
write.csv(index.iid.tw.120kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_tw_120kn.csv"))
write.csv(index.iid.110kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_110kn_2024_12_02.csv"))
write.csv(index.iid.100kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_100kn_2024_12_02.csv"))
write.csv(index.iid.tw.90kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_tw_90kn.csv"))
write.csv(index.iid.75kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_75kn_2024_12_02.csv"))
write.csv(index.iid.tw.50kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_tw_50kn.csv"))
write.csv(index.iid.25kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_25kn_2024_12_02.csv"))

# save RW/Tweedie index values
write.csv(index.rw.tw.120kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_rw_tw_120kn.csv"))
write.csv(index.rw.tw.90kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_rw_tw_90kn.csv"))
write.csv(index.rw.tw.50kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_rw_tw_50kn.csv"))

# save AR1/Tweedie index values
write.csv(index.ar1.tw.120kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_tw_120kn.csv"))
write.csv(index.ar1.tw.90kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_tw_90kn.csv"))
write.csv(index.ar1.tw.50kn.t, paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_tw_50kn.csv"))

# bring in saved indices
index.iid.tw.120kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_tw_120kn.csv"))
index.iid.tw.90kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_tw_90kn.csv"))
index.iid.tw.50kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_iid_tw_50kn.csv"))
index.rw.tw.120kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_rw_tw_120kn.csv"))
index.rw.tw.90kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_rw_tw_90kn.csv"))
index.rw.tw.50kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_rw_tw_50kn.csv"))
index.ar1.tw.120kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_tw_120kn.csv"))
index.ar1.tw.90kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_tw_90kn.csv"))
index.ar1.tw.50kn.t <- read.csv(paste0(here::here(), "/SMBKC/output/smbkc_ebs_trawl_index_ar1_tw_50kn.csv"))

# plot IID/Tweedie indices
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

# combine IID/Tweedie indices in one data frame for plotting
mod.compare <- rbind(index.t %>% mutate(model = "AR1, 120 knots"), index.50kn.t %>% mutate(model = "AR1, 50 knots")) %>%
  rbind(index.iid.t %>% mutate(model = "IID, 120 knots"))

mod.compare.iid <- rbind(index.iid.tw.120kn.t %>% mutate(model = "IID, Tweedie, 120 knots"), index.iid.tw.90kn.t %>% mutate(model = "IID, Tweedie, 90 knots")) %>%
  #rbind(index.iid.110kn.t %>% mutate(model = "110 knots")) %>%
  #rbind(index.iid.100kn.t %>% mutate(model = "100 knots")) %>%
  rbind(index.iid.tw.50kn.t %>% mutate(model = "IID, Tweedie, 50 knots")) #%>%
  #rbind(index.iid.75kn.t %>% mutate(model = "75 knots")) %>%
  #rbind(index.iid.25kn.t %>% mutate(model = "25 knots"))

# combine RW/Tweedie indices in one data frame for plotting
mod.compare.rw <- rbind(index.rw.tw.120kn.t %>% mutate(model = "RW, Tweedie, 120 knots"), index.rw.tw.90kn.t %>% mutate(model = "RW, Tweedie, 90 knots")) %>%
  rbind(index.rw.tw.50kn.t %>% mutate(model = "RW, Tweedie, 50 knots")) 

# combine AR1/Tweedie indices in one data frame for plotting
mod.compare.ar1 <- rbind(index.ar1.tw.120kn.t %>% mutate(model = "AR1, Tweedie, 120 knots"), index.ar1.tw.90kn.t %>% mutate(model = "AR1, Tweedie, 90 knots")) %>%
  rbind(index.ar1.tw.50kn.t %>% mutate(model = "AR1, Tweedie, 50 knots")) 

# combine model with best predictive skill from each estimation method into one data frame
mod.compare.50kn <- rbind(index.iid.tw.50kn.t %>% mutate(model = "IID, Tweedie, 50 knots"), index.rw.tw.50kn.t %>% mutate(model = "RW, Tweedie, 50 knots")) %>%
  rbind(index.ar1.tw.50kn.t %>% mutate(model = "AR1, Tweedie, 50 knots")) 

# plot IID/Tweedie indices
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




# ************************************************************************************************
# compare predicted index to observations ----
# ************************************************************************************************

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
obs.pred.iid.tw <- 
  left_join(mod.compare.iid, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

obs.pred.rw.tw <- 
  left_join(mod.compare.rw, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

obs.pred.ar1.tw <- 
  left_join(mod.compare.ar1, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

obs.pred.50kn.tw <- 
  left_join(mod.compare.50kn, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

sm.fit.iid.tw <- ggplot(obs.pred.iid.tw)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") 

ggsave(file.path(plotdir.sm, "smbkc_iid_tw_fit.png"), plot = sm.fit.iid.tw, height = 5, width = 7, units = "in")

sm.fit.rw.tw <- ggplot(obs.pred.rw.tw)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") 

ggsave(file.path(plotdir.sm, "smbkc_rw_tw_fit.png"), plot = sm.fit.rw.tw, height = 5, width = 7, units = "in")

sm.fit.ar1.tw <- ggplot(obs.pred.ar1.tw)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") 

ggsave(file.path(plotdir.sm, "smbkc_ar1_tw_fit.png"), plot = sm.fit.ar1.tw, height = 5, width = 7, units = "in")

sm.fit.50kn.tw <- ggplot(obs.pred.50kn.tw)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank())

ggsave(file.path(plotdir.sm, "smbkc_50kn_tw_fit.png"), plot = sm.fit.50kn.tw, height = 5, width = 7, units = "in")

sm.fit.50kn.tw.noci <- ggplot(obs.pred.50kn.tw)+
  geom_line(aes(x = SURVEY_YEAR, y = est.t, color = model))+
  #geom_ribbon(aes(x = SURVEY_YEAR, y = est.t, ymin = lwr.t, ymax = upr.t, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank())

ggsave(file.path(plotdir.sm, "smbkc_50kn_tw_fit_noci.png"), plot = sm.fit.50kn.tw, height = 5, width = 7, units = "in")

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
# plot VAST indices ----
# ******************************************************************************

# VAST IID/delta gamma

index.vast.iid.dg.120kn <- read.csv(paste0(here::here(), "/SMBKC/VAST output/SMBKC_GE90/SMBKC_GE90_ObsModel_2_1_120kts_V2025_Epsilon2_Omega2_disabled/Index.csv")) %>% 
  mutate(SURVEY_YEAR = Time, vast.est = Estimate, vast.est.se = Std..Error.for.Estimate) %>%
  mutate(vast.lwr = vast.est - 1.96*vast.est.se, vast.upr = vast.est + 1.96*vast.est.se, model = "IID, delta gamma, 120 knots") 
index.vast.iid.dg.90kn <- read.csv(paste0(here::here(), "/SMBKC/VAST output/SMBKC_GE90/SMBKC_GE90_ObsModel_2_1_90kts_V2025_Epsilon2_Omega2_disabled/Index.csv")) %>% 
  mutate(SURVEY_YEAR = Time, vast.est = Estimate, vast.est.se = Std..Error.for.Estimate) %>%
  mutate(vast.lwr = vast.est - 1.96*vast.est.se, vast.upr = vast.est + 1.96*vast.est.se, model = "IID, delta gamma, 90 knots") 
index.vast.iid.dg.50kn <- read.csv(paste0(here::here(), "/SMBKC/VAST output/SMBKC_GE90/SMBKC_GE90_ObsModel_2_1_50kts_V2025_Epsilon2_Omega2_disabled/Index.csv")) %>% 
  mutate(SURVEY_YEAR = Time, vast.est = Estimate, vast.est.se = Std..Error.for.Estimate) %>%
  mutate(vast.lwr = vast.est - 1.96*vast.est.se, vast.upr = vast.est + 1.96*vast.est.se, model = "IID, delta gamma, 50 knots") 

index.vast.iid.dg <- rbind(index.vast.iid.dg.120kn, index.vast.iid.dg.90kn, index.vast.iid.dg.50kn)

obs.pred.vast.iid.dg <- 
  left_join(index.vast.iid.dg, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

sm.fit.vast.iid.dg <- ggplot(obs.pred.vast.iid.dg)+
  geom_line(aes(x = SURVEY_YEAR, y = vast.est, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = vast.est, ymin = vast.lwr, ymax = vast.upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  ggtitle("VAST")

ggsave(file.path(plotdir.sm, "smbkc_vast_iid_dg_fit.png"), plot = sm.fit.vast.iid.dg, height = 5, width = 7, units = "in")

# VAST IID/Tweedie 

index.vast.iid.tw.120kn <- read.csv(paste0(here::here(), "/SMBKC/VAST output/SMBKC_GE90/SMBKC_GE90_ObsModel_10_2_120kts_V2025/Index.csv")) %>% 
  mutate(SURVEY_YEAR = Time, vast.est = Estimate, vast.est.se = Std..Error.for.Estimate) %>%
  mutate(vast.lwr = vast.est - 1.96*vast.est.se, vast.upr = vast.est + 1.96*vast.est.se, model = "IID, Tweedie, 120 knots") 
index.vast.iid.tw.90kn <- read.csv(paste0(here::here(), "/SMBKC/VAST output/SMBKC_GE90/SMBKC_GE90_ObsModel_10_2_90kts_V2025/Index.csv")) %>% 
  mutate(SURVEY_YEAR = Time, vast.est = Estimate, vast.est.se = Std..Error.for.Estimate) %>%
  mutate(vast.lwr = vast.est - 1.96*vast.est.se, vast.upr = vast.est + 1.96*vast.est.se, model = "IID, Tweedie, 90 knots") 
index.vast.iid.tw.50kn <- read.csv(paste0(here::here(), "/SMBKC/VAST output/SMBKC_GE90/SMBKC_GE90_ObsModel_10_2_50kts_V2025/Index.csv")) %>% 
  mutate(SURVEY_YEAR = Time, vast.est = Estimate, vast.est.se = Std..Error.for.Estimate) %>%
  mutate(vast.lwr = vast.est - 1.96*vast.est.se, vast.upr = vast.est + 1.96*vast.est.se, model = "IID, Tweedie, 50 knots") 

index.vast.iid.tw <- rbind(index.vast.iid.tw.120kn, index.vast.iid.tw.90kn, index.vast.iid.tw.50kn)

obs.pred.vast.iid.tw <- 
  left_join(index.vast.iid.tw, obs.biom, by = "SURVEY_YEAR") %>% 
  filter(SURVEY_YEAR >= 1978) %>%
  mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
         obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))

sm.fit.vast.iid.tw <- ggplot(obs.pred.vast.iid.tw)+
  geom_line(aes(x = SURVEY_YEAR, y = vast.est, color = model))+
  geom_ribbon(aes(x = SURVEY_YEAR, y = vast.est, ymin = vast.lwr, ymax = vast.upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = SURVEY_YEAR, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = SURVEY_YEAR, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  ggtitle("VAST")

ggsave(file.path(plotdir.sm, "smbkc_vast_iid_tw_fit.png"), plot = sm.fit.vast.iid.tw, height = 5, width = 7, units = "in")

# from Jon:
# 10_2 IID/Tweedie
# delta gamma 2_1
# no AR1 (crashed), no random walk (not sure how to implement)


# ******************************************************************************
# compare sdmTMB and VAST predicted indices ----
# ******************************************************************************

vast.sdmTMB.index.smbkc <- rbind(
  index.vast.iid.tw.50kn %>% 
    mutate(model = "VAST, IID, Tweedie, 50 kn") %>%
    mutate(Year = Time, est = vast.est, lwr = vast.lwr, upr = vast.upr) %>%
    select(c(Year, est, lwr, upr, model)),
  index.vast.iid.tw.90kn %>% 
    mutate(model = "VAST, IID, Tweedie, 90 kn") %>%
    mutate(Year = Time, est = vast.est, lwr = vast.lwr, upr = vast.upr) %>%
    select(c(Year, est, lwr, upr, model)),
  index.rw.tw.50kn.t %>% 
    mutate(Year = SURVEY_YEAR, model = "sdmTMB, RW, Tweedie, 50 knots") %>%
    mutate(est = est/1000, lwr = lwr/1000, upr = upr/1000) %>%
    select(c(Year, est, lwr, upr, model)),
  index.rw.tw.90kn.t %>% 
    mutate(Year = SURVEY_YEAR, model = "sdmTMB, RW, Tweedie, 90 knots") %>%
    mutate(est = est/1000, lwr = lwr/1000, upr = upr/1000) %>%
    select(c(Year, est, lwr, upr, model))
)

vast.sdmTMB.obs.pred.smbkc <- left_join(vast.sdmTMB.index.smbkc, 
                                        obs.biom %>% 
                                          rename(Year = SURVEY_YEAR) %>%
                                          mutate(obs_l95 = BIOMASS_MT * exp(-1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))),
                                                 obs_u95 = BIOMASS_MT * exp(1.96 * sqrt(log(1 + BIOMASS_MT_CV^2))))) %>%
  filter(Year >= 1978)


vast.sdmTMB.smbkc.fit.50kn <- ggplot(vast.sdmTMB.obs.pred.smbkc %>%
                                       filter(model %in% c("VAST, IID, Tweedie, 50 kn", "sdmTMB, RW, Tweedie, 50 knots")))+
  geom_line(aes(x = Year, y = est, color = model))+
  geom_ribbon(aes(x = Year, y = est, ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = Year, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom")

ggsave(file.path(plotdir.sm, "vast_sdmTMB_smbkc_fit_50kn.png"), plot = vast.sdmTMB.smbkc.fit.50kn, height = 5, width = 7, units = "in")

vast.sdmTMB.smbkc.fit.90kn <- ggplot(vast.sdmTMB.obs.pred.smbkc %>%
                                       filter(model %in% c("VAST, IID, Tweedie, 90 kn", "sdmTMB, RW, Tweedie, 90 knots")))+
  geom_line(aes(x = Year, y = est, color = model))+
  geom_ribbon(aes(x = Year, y = est, ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = Year, y = BIOMASS_MT), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Blue king crab biomass estimate (t)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom")

ggsave(file.path(plotdir.sm, "vast_sdmTMB_smbkc_fit_90kn.png"), plot = vast.sdmTMB.smbkc.fit.90kn, height = 5, width = 7, units = "in")



