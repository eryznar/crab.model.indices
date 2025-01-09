# ************************************************************************************************
# Generating a spatiotemporal model-based index for Norton Sound red king crab 
# in the Northern Bering Sea survey only (for comparison with VAST estimates)
# December 2024
# Caitlin Stern
# ************************************************************************************************

# ************************************************************************************************
# load libraries, set plot preferences ----
# ************************************************************************************************

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sp)
library(pROC)

cur_yr.ns <- 2024

plotdir.ns <- paste0(here::here(), "/NSRKC/plots")
modeldir.ns <- paste0(here::here(), "/NSRKC/models")
outdir.ns <- paste0(here::here(), "/NSRKC/output")

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

grid2km.ns <- readRDS(paste0(here::here(), "/NSRKC/data/NSRKC_Extrapolation_areas/NSRKC_5km_grid.rds")) %>%
  # get UTM zone
  mutate(zone = (floor((Lon + 180)/6) %% 60) + 1) 

unique(grid2km.ns$zone)

predgrid_utm_z3.ns <- grid2km.ns %>% filter(zone == 3)
predgrid_utm_z4.ns <- grid2km.ns %>% filter(zone == 4)

u3.ns <- predgrid_utm_z3.ns
get_crs(u3.ns, c("Lon", "Lat"))
u3u.ns <- add_utm_columns(u3.ns, c("Lon", "Lat"))

u4.ns <- predgrid_utm_z4.ns
get_crs(u4.ns, c("Lon", "Lat"))
u4u.ns <- add_utm_columns(u4.ns, c("Lon", "Lat"))

predgrid_utm.ns <- rbind(u3u.ns, u4u.ns) %>%
  mutate(X = X / 1000, Y = Y / 1000) # convert UTM coordinates from meter to kilometers to reduce computing needs

# plot prediction grid
predgrid5kmplot <- ggplot(predgrid_utm.ns, aes(x = Lon, y = Lat, color = Area_km2)) +
  geom_point() + 
  theme(legend.position="none") +
  labs(x = "Longitude", y = "Latitude")

ggsave(file.path(plotdir.ns, "prediction_grid_5km.png"), plot = predgrid5kmplot, height = 4.2, width = 7, units = "in")

# *************************************************************************************************
# read in and process survey data ----
# *************************************************************************************************

nsrkc.survey <- read.csv(paste0(here::here(), "/NSRKC/data/NSRKC_trawl_survey_abundance.csv"))

nsrkc.dt <- nsrkc.survey %>%
  #select(Year, ADFG_Station, Latitude, Longitude, crab.km2) %>%
  arrange(Year, ADFG_Station) %>%
  mutate("year_f" = factor(Year)) %>%
  select(-X)

# check that all Lat/Lon combinations are valid. Norton Sound includes 
# 160 < longitude < 168
# 61.49 < latitude < 66 (source: https://www.adfg.alaska.gov/static/applications/dcfnewsrelease/750733004.pdf)
filter(nsrkc.dt, Latitude > 66 | Latitude < 61.49)
filter(nsrkc.dt, Longitude > -160 | Longitude < -168)

nsrkc.dt2 <- nsrkc.dt %>%
  filter(Latitude <= 65.9) %>%
  filter(Latitude >= 61.49) %>%
  filter(Longitude <= -160) %>%
  filter(Longitude >= -168)

nsrkc_utm_zone <- nsrkc.dt2 %>%  
  # get UTM zone
  mutate(zone = (floor((Longitude + 180)/6) %% 60) + 1) 

unique(nsrkc_utm_zone$zone)

nsrkc_utm_z3 <- nsrkc_utm_zone %>% filter(zone == 3)
nsrkc_utm_z4 <- nsrkc_utm_zone %>% filter(zone == 4)

d3 <- nsrkc_utm_z3
get_crs(d3, c("Longitude", "Latitude"))
d3u <- add_utm_columns(d3, c("Longitude", "Latitude"))

d4 <- nsrkc_utm_z4
get_crs(d4, c("Longitude", "Latitude"))
d4u <- add_utm_columns(d4, c("Longitude", "Latitude"))

nsrkc_utm <- rbind(d3u, d4u) %>%
  mutate(X = X / 1000, Y = Y / 1000) # convert UTM coordinates from meter to kilometers to reduce computing needs

# include only data from the Northern Bering Sea survey
nsrkc_utm_nbs <- nsrkc_utm %>%
  filter(Agent == "NBS")

# number of unique stations
length(unique(nsrkc_utm_nbs$ADFG_Station))

# check number of observations per year (this should be larger than the number of vertices in the mesh)
nsrkc_utm_nbs %>%
  group_by(Year) %>%
  count() %>%
  print(n = 50)

# plot survey data
ggplot(nsrkc_utm_nbs) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

ggplot(nsrkc_utm_nbs) + 
  geom_point(aes(x = Longitude, y = Latitude)) +
  theme_bw()


# total abundance

abund.crab.nbs <- nsrkc_utm_nbs %>%
  #filter(Agent == "NBS") %>%
  group_by(Agent, Year) %>%
  mutate(crab.ab = crab.km2 * area.km2) %>%
  summarise(total.ab = sum(crab.ab),
            mean.crab.ab = mean(crab.ab), 
            sd.crab.ab = sd(crab.ab, na.rm = TRUE),
            n.station = length(ADFG_Station)) %>%
  mutate(se.crab.ab = sd.crab.ab / sqrt(n.station)) %>%
  mutate(obs.mean.lwr95 = mean.crab.ab - 1.96*se.crab.ab,
         obs.mean.upr95 = mean.crab.ab + 1.96*se.crab.ab,
         obs.lwr95 = obs.mean.lwr95 * n.station,
         obs.upr95 = obs.mean.upr95 * n.station)


# *************************************************************************************************
# make SPDE mesh ----
# *************************************************************************************************

nbs_mesh_100kn <- make_mesh(nsrkc_utm_nbs, xy_cols = c("X","Y"), n_knots = 100, type = "kmeans")

plot(nbs_mesh_100kn)
nbs_mesh_100kn$mesh$n

nbs_mesh_50kn <- make_mesh(nsrkc_utm_nbs, xy_cols = c("X","Y"), n_knots = 50, type = "kmeans")

plot(nbs_mesh_50kn)
nbs_mesh_50kn$mesh$n

nbs_mesh_30kn <- make_mesh(nsrkc_utm_nbs, xy_cols = c("X","Y"), n_knots = 30, type = "kmeans")

plot(nbs_mesh_30kn)
nbs_mesh_30kn$mesh$n


max.edge = diff(range(st_coordinates(nsrkc_utm_nbs)[,1]))/(3*5)
inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(nsrkc_utm$X, nsrkc_utm$Y), # coordinates
  max.edge = c(1,2)*max.edge, # max triangle edge length; inner and outer meshes
  #offset = c(5, 25),  # inner and outer border widths
  #cutoff = 5 # minimum triangle edge length
)

mesh <- make_mesh(nsrkc_utm_nbs, c("X", "Y"), mesh = inla_mesh)
plot(mesh)

# plot mesh with data points

mesh.plot.100kn.nbs <- ggplot(nsrkc_utm_nbs) + 
  inlabru::gg(nbs_mesh_100kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

mesh.plot.50kn.nbs <- ggplot(nsrkc_utm_nbs) + 
  inlabru::gg(nbs_mesh_50kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

mesh.plot.30kn.nbs <- ggplot(nsrkc_utm_nbs) + 
  inlabru::gg(nbs_mesh_30kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

# export mesh plots

ggsave(file.path(plotdir.ns, "mesh_100kn_nbs.png"), plot = mesh.plot.100kn.nbs, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "mesh_50kn_nbs.png"), plot = mesh.plot.50kn.nbs, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "mesh_30kn_nbs.png"), plot = mesh.plot.30kn.nbs, height = 4.2, width = 7, units = "in")

# *************************************************************************************************
# fit and check models ----
# *************************************************************************************************

# IID/Tweedie models ----

# 100 knots
m.nbs.iid.tw.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# 50 knots
m.nbs.iid.tw.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# 30 knots
m.nbs.iid.tw.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# IID/delta gamma models ----

# 100 knots
m.nbs.iid.dg.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# 50 knots
m.nbs.iid.dg.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# 30 knots
m.nbs.iid.dg.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# IID/delta lognormal models ----

# 100 knots
m.nbs.iid.dl.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# 50 knots
m.nbs.iid.dl.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# 30 knots
m.nbs.iid.dl.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "iid",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# RW/Tweedie models ----

# 100 knots
m.nbs.rw.tw.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# 50 knots
m.nbs.rw.tw.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# 30 knots
m.nbs.rw.tw.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# RW/delta gamma models ----

# 100 knots
m.nbs.rw.dg.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# 50 knots
m.nbs.rw.dg.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# 30 knots
m.nbs.rw.dg.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# Random walk/delta lognormal models ----

# 100 knots
m.nbs.rw.dl.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# 50 knots
m.nbs.rw.dl.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# 30 knots
m.nbs.rw.dl.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# AR1/Tweedie models ----

# 100 knots
m.nbs.ar1.tw.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# 50 knots
m.nbs.ar1.tw.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# 30 knots
m.nbs.ar1.tw.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))


# AR1/delta gamma models ----

# 100 knots
m.nbs.ar1.dg.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# 50 knots
m.nbs.ar1.dg.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))

# 30 knots
m.nbs.ar1.dg.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_gamma(type = "standard"))


# AR1/delta lognormal models ----

# 100 knots
m.nbs.ar1.dl.100kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# 50 knots
m.nbs.ar1.dl.50kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# 30 knots
m.nbs.ar1.dl.30kn <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"))

# save fitted models ----

# save the fitted IID/Tweedie models
saveRDS(m.nbs.iid.tw.100kn, file = paste0(modeldir.ns, "/m_nbs_iid_tw_100kn.RDS"))
saveRDS(m.nbs.iid.tw.50kn, file = paste0(modeldir.ns, "/m_nbs_iid_tw_50kn.RDS"))
saveRDS(m.nbs.iid.tw.30kn, file = paste0(modeldir.ns, "/m_nbs_iid_tw_30kn.RDS"))

# save the fitted IID/delta gamma models
saveRDS(m.nbs.iid.dg.100kn, file = paste0(modeldir.ns, "/m_nbs_iid_dg_100kn.RDS"))
saveRDS(m.nbs.iid.dg.50kn, file = paste0(modeldir.ns, "/m_nbs_iid_dg_50kn.RDS"))
saveRDS(m.nbs.iid.dg.30kn, file = paste0(modeldir.ns, "/m_nbs_iid_dg_30kn.RDS"))

# save the fitted IID/delta lognormal models
saveRDS(m.nbs.iid.dl.100kn, file = paste0(modeldir.ns, "/m_nbs_iid_dl_100kn.RDS"))
saveRDS(m.nbs.iid.dl.50kn, file = paste0(modeldir.ns, "/m_nbs_iid_dl_50kn.RDS"))
saveRDS(m.nbs.iid.dl.30kn, file = paste0(modeldir.ns, "/m_nbs_iid_dl_30kn.RDS"))

# save the fitted RW/Tweedie models
saveRDS(m.nbs.rw.tw.100kn, file = paste0(modeldir.ns, "/m_nbs_rw_tw_100kn.RDS"))
saveRDS(m.nbs.rw.tw.50kn, file = paste0(modeldir.ns, "/m_nbs_rw_tw_50kn.RDS"))
saveRDS(m.nbs.rw.tw.30kn, file = paste0(modeldir.ns, "/m_nbs_rw_tw_30kn.RDS"))

# save the fitted RW/delta gamma models
saveRDS(m.nbs.rw.dg.100kn, file = paste0(modeldir.ns, "/m_nbs_rw_dg_100kn.RDS"))
saveRDS(m.nbs.rw.dg.50kn, file = paste0(modeldir.ns, "/m_nbs_rw_dg_50kn.RDS"))
saveRDS(m.nbs.rw.dg.30kn, file = paste0(modeldir.ns, "/m_nbs_rw_dg_30kn.RDS"))

# save the fitted RW/delta lognormal models
saveRDS(m.nbs.rw.dl.100kn, file = paste0(modeldir.ns, "/m_nbs_rw_dl_100kn.RDS"))
saveRDS(m.nbs.rw.dl.50kn, file = paste0(modeldir.ns, "/m_nbs_rw_dl_50kn.RDS"))
saveRDS(m.nbs.rw.dl.30kn, file = paste0(modeldir.ns, "/m_nbs_rw_dl_30kn.RDS"))

# save the fitted AR1/Tweedie models
saveRDS(m.nbs.ar1.tw.100kn, file = paste0(modeldir.ns, "/m_nbs_ar1_tw_100kn.RDS"))
saveRDS(m.nbs.ar1.tw.50kn, file = paste0(modeldir.ns, "/m_nbs_ar1_tw_50kn.RDS"))
saveRDS(m.nbs.ar1.tw.30kn, file = paste0(modeldir.ns, "/m_nbs_ar1_tw_30kn.RDS"))

# save the fitted AR1/delta gamma models
saveRDS(m.nbs.ar1.dg.100kn, file = paste0(modeldir.ns, "/m_nbs_ar1_dg_100kn.RDS"))
saveRDS(m.nbs.ar1.dg.50kn, file = paste0(modeldir.ns, "/m_nbs_ar1_dg_50kn.RDS"))
saveRDS(m.nbs.ar1.dg.30kn, file = paste0(modeldir.ns, "/m_nbs_ar1_dg_30kn.RDS"))

# save the fitted AR1/delta lognormal models
saveRDS(m.nbs.ar1.dl.100kn, file = paste0(modeldir.ns, "/m_nbs_ar1_dl_100kn.RDS"))
saveRDS(m.nbs.ar1.dl.50kn, file = paste0(modeldir.ns, "/m_nbs_ar1_dl_50kn.RDS"))
saveRDS(m.nbs.ar1.dl.30kn, file = paste0(modeldir.ns, "/m_nbs_ar1_dl_30kn.RDS"))

# run sanity checks ----

# IID models
sanity(m.nbs.iid.tw.100kn) # doesn't pass
sanity(m.nbs.iid.tw.50kn) # doesn't pass
sanity(m.nbs.iid.tw.30kn) # doesn't pass
sanity(m.nbs.iid.dg.100kn) # doesn't pass
sanity(m.nbs.iid.dg.50kn) # doesn't pass
sanity(m.nbs.iid.dg.30kn) # doesn't pass
sanity(m.nbs.iid.dl.100kn) # doesn't pass
sanity(m.nbs.iid.dl.50kn) # doesn't pass
sanity(m.nbs.iid.dl.30kn) # doesn't pass

# RW models
sanity(m.nbs.rw.tw.100kn) # doesn't pass but almost
sanity(m.nbs.rw.tw.50kn) # doesn't pass but almost
sanity(m.nbs.rw.tw.30kn) # passes
sanity(m.nbs.rw.dg.100kn) # doesn't pass
sanity(m.nbs.rw.dg.50kn) # doesn't pass
sanity(m.nbs.rw.dg.30kn) # doesn't pass
sanity(m.nbs.rw.dl.100kn) # doesn't pass
sanity(m.nbs.rw.dl.50kn) # doesn't pass but almost
sanity(m.nbs.rw.dl.30kn) # doesn't pass

# AR1 models
sanity(m.nbs.ar1.tw.100kn) # doesn't pass
sanity(m.nbs.ar1.tw.50kn) # doesn't pass
sanity(m.nbs.ar1.tw.30kn) # doesn't pass but almost
sanity(m.nbs.ar1.dg.100kn) # doesn't pass
sanity(m.nbs.ar1.dg.50kn) # doesn't pass
sanity(m.nbs.ar1.dg.30kn) # doesn't pass
sanity(m.nbs.ar1.dl.100kn) # doesn't pass
sanity(m.nbs.ar1.dl.50kn) # doesn't pass
sanity(m.nbs.ar1.dl.30kn) # doesn't pass


# ************************************************************************************************
# compare models using cross-validation ----
# ************************************************************************************************

# cross-validation for models that passed or almost passed sanity checks

# RW, Tweedie, 100 knots
m.nbs.rw.tw.100kn.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_100kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)

# RW, Tweedie, 50 knots
m.nbs.rw.tw.50kn.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)

# RW, Tweedie, 30 knots
m.nbs.rw.tw.30kn.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)

# RW, delta lognormal, 50 knots
m.nbs.rw.dl.50kn.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_50kn, 
  spatiotemporal = "rw",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = delta_lognormal(type = "standard"),
  k_folds = 3)

# AR1, Tweedie, 30 knots
m.nbs.ar1.tw.30kn.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = nbs_mesh_30kn, 
  spatiotemporal = "ar1",
  # NBS survey took place in 2010, 2017, 2019, 2021, 2022, 2023
  extra_time = c(2011, 2012, 2013, 2014, 2015, 2016, 2018, 2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"),
  k_folds = 3)

m.nbs.rw.tw.100kn.cv$sum_loglik
m.nbs.rw.tw.50kn.cv$sum_loglik
m.nbs.rw.tw.30kn.cv$sum_loglik


# save cross-validation models ----
saveRDS(m.nbs.rw.tw.100kn.cv, file = paste0(modeldir.ns, "/m_nbs_rw_tw_100kn_cv.RDS"))
saveRDS(m.nbs.rw.tw.50kn.cv, file = paste0(modeldir.ns, "/m_nbs_rw_tw_50kn_cv.RDS"))
saveRDS(m.nbs.rw.tw.30kn.cv, file = paste0(modeldir.ns, "/m_nbs_rw_tw_30kn_cv.RDS"))

# make table with loglikelihood results ----
m.nsrkc.compare.loglik <- data.frame(family=character(),
                                     estimation=character(), 
                                     knots=numeric(), 
                                     loglik=numeric(),
                                     stringsAsFactors=FALSE) %>%
  # IID models
  add_row(family = "Tweedie", estimation = "IID", knots = 100, loglik = NA) %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 50, loglik = NA) %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 30, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 100, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 50, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 30, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 100, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 50, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 30, loglik = NA) %>%
  # RW models
  add_row(family = "Tweedie", estimation = "RW", knots = 100, loglik = NA) %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 50, loglik = NA) %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 30, loglik = m.nbs.rw.tw.30kn.cv$sum_loglik) %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 100, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 50, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 30, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 100, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 50, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 30, loglik = NA) %>%
  # AR1 models
  add_row(family = "Tweedie", estimation = "AR1", knots = 100, loglik = NA) %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 50, loglik = NA) %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 30, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 100, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 50, loglik = NA) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 30, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 100, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 50, loglik = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 30, loglik = NA)

# save table with predictive skill values
write.csv(m.nsrkc.compare.loglik, paste0(here::here(), "/NSRKC/output/nsrkc_compare_loglik.csv"))

# ************************************************************************************************
# model residuals ----
# *************************************************************************************************

m.nbs.rw.tw.100kn.resid <- simulate(update(m.nbs.rw.tw.100kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nbs.rw.tw.100kn, return_DHARMa = TRUE)

plot(m.nbs.rw.tw.100kn.resid)

m.nbs.rw.tw.50kn.resid <- simulate(update(m.nbs.rw.tw.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nbs.rw.tw.50kn, return_DHARMa = TRUE)

plot(m.nbs.rw.tw.50kn.resid)

m.nbs.rw.tw.30kn.resid <- simulate(update(m.nbs.rw.tw.30kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nbs.rw.tw.30kn, return_DHARMa = TRUE)

plot(m.nbs.rw.tw.30kn.resid)

# save residuals
saveRDS(m.nbs.rw.tw.100kn.resid, file = paste0(modeldir.ns, "/m_nbs_resid_rw_tw_100kn.RDS"))
saveRDS(m.nbs.rw.tw.50kn.resid, file = paste0(modeldir.ns, "/m_nbs_resid_rw_tw_50kn.RDS"))
saveRDS(m.nbs.rw.tw.30kn.resid, file = paste0(modeldir.ns, "/m_nbs_resid_rw_tw_30kn.RDS"))

# residuals tests
# 100 knots
DHARMa::testOutliers(m.nbs.rw.tw.100kn.resid) 
DHARMa::testQuantiles(m.nbs.rw.tw.100kn.resid) 
DHARMa::testDispersion(m.nbs.rw.tw.100kn.resid) 
DHARMa::testResiduals(m.nbs.rw.tw.100kn.resid)
DHARMa::testZeroInflation(m.nbs.rw.tw.100kn.resid) 
# 50 knots
DHARMa::testOutliers(m.nbs.rw.tw.50kn.resid) 
DHARMa::testQuantiles(m.nbs.rw.tw.50kn.resid) 
DHARMa::testDispersion(m.nbs.rw.tw.50kn.resid) 
DHARMa::testResiduals(m.nbs.rw.tw.50kn.resid)
DHARMa::testZeroInflation(m.nbs.rw.tw.50kn.resid) 
# 30 knots
DHARMa::testOutliers(m.nbs.rw.tw.30kn.resid) # not significant
DHARMa::testQuantiles(m.nbs.rw.tw.30kn.resid) # quantile deviations detected
DHARMa::testDispersion(m.nbs.rw.tw.30kn.resid) # not significant
DHARMa::testResiduals(m.nbs.rw.tw.30kn.resid)
DHARMa::testZeroInflation(m.nbs.rw.tw.30kn.resid) # not significant

# add the scaled residuals to one data frame and then make spatial residual plots ----

nsrkc_utm_nbs_resids <- cbind(nsrkc_utm_nbs, 
                             data.frame(m.nbs.rw.tw.30kn.resid$scaledResiduals),
                             data.frame(m.nbs.rw.tw.50kn.resid$scaledResiduals),
                             data.frame(m.nbs.rw.tw.100kn.resid$scaledResiduals))


spat.resid.nbs.rw.tw.30kn <- ggplot(nsrkc_utm_nbs_resids) + 
  geom_point(aes(y = Latitude, x = Longitude, color = m.nbs.rw.tw.30kn.resid$scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("NSRKC RW Tweedie model residuals (30 knots)")+
  facet_wrap(~Year)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

ggsave(paste0(plotdir.ns, "/spat_resid_nbs_rw_tw_30kn.png"), spat.resid.nbs.rw.tw.30kn, height = 5, width = 7, units = "in")

spat.resid.nbs.rw.tw.50kn <- ggplot(nsrkc_utm_nbs_resids) + 
  geom_point(aes(y = Latitude, x = Longitude, color = m.nbs.rw.tw.50kn.resid$scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("NSRKC RW Tweedie model residuals (50 knots)")+
  facet_wrap(~Year)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

ggsave(paste0(plotdir.ns, "/spat_resid_nbs_rw_tw_50kn.png"), spat.resid.nbs.rw.tw.50kn, height = 5, width = 7, units = "in")

spat.resid.nbs.rw.tw.100kn <- ggplot(nsrkc_utm_nbs_resids) + 
  geom_point(aes(y = Latitude, x = Longitude, color = m.nbs.rw.tw.100kn.resid$scaledResiduals), size = 1) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  ggtitle("NSRKC RW Tweedie model residuals (30 knots)")+
  facet_wrap(~Year)+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom") +
  labs(color = "scaled residuals")

ggsave(paste0(plotdir.ns, "/spat_resid_nbs_rw_tw_100kn.png"), spat.resid.nbs.rw.tw.100kn, height = 5, width = 7, units = "in")

# make table with residual analysis results ----
m.nsrkc.compare.resid <- data.frame(family=character(),
                                    estimation=character(), 
                                    knots=numeric(), 
                                    quantile=character(),
                                    dispersion=character(),
                                    outliers=character(),
                                    zero_inf=character(),
                                    stringsAsFactors=FALSE) %>%
  # IID models
  add_row(family = "Tweedie", estimation = "IID", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "Tweedie", estimation = "IID", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "IID", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "IID", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  # RW models
  add_row(family = "Tweedie", estimation = "RW", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "Tweedie", estimation = "RW", knots = 30, quantile = "sig.", dispersion = "n.s.", outliers = "n.s.", zero_inf = "n.s.") %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "RW", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "RW", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  # AR1 models
  add_row(family = "Tweedie", estimation = "AR1", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "Tweedie", estimation = "AR1", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta gamma", estimation = "AR1", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 100, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
  add_row(family = "delta lognormal", estimation = "AR1", knots = 50, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA) %>%
add_row(family = "delta lognormal", estimation = "AR1", knots = 30, quantile = NA, dispersion = NA, outliers = NA, zero_inf = NA)

# save table with residuals results
write.csv(m.nsrkc.compare.resid, paste0(here::here(), "/NSRKC/output/nsrkc_compare_resid.csv"))

# ************************************************************************************************
# generate predictions ----
# *************************************************************************************************

# replicate grid across all years
pred.grid.nbs <- replicate_df(predgrid_utm.ns, "year_f", unique(nsrkc_utm_nbs$year_f))
pred.grid.nbs$year <- as.integer(as.character(factor(pred.grid.nbs$year_f)))
dplyr::glimpse(pred.grid.nbs)

predgrid.utm.nbs <- replicate_df(predgrid_utm.ns, "Year", unique(nsrkc_utm_nbs$Year))
predgrid.utm.nbs$year_f <- factor(predgrid.utm.nbs$Year)
dplyr::glimpse(predgrid.utm.nbs)

# predictions on new data ----
pred.nbs.rw.tw.100kn <- predict(m.nbs.rw.tw.100kn, newdata = predgrid.utm.nbs, return_tmb_object = T)
pred.nbs.rw.tw.50kn <- predict(m.nbs.rw.tw.50kn, newdata = predgrid.utm.nbs, return_tmb_object = T)
pred.nbs.rw.tw.30kn <- predict(m.nbs.rw.tw.30kn, newdata = predgrid.utm.nbs, return_tmb_object = T)

# save predictions
saveRDS(pred.nbs.rw.tw.100kn, file = paste0(modeldir.ns, "/m_nbs_pred_rw_tw_100kn.RDS"))
saveRDS(pred.nbs.rw.tw.50kn, file = paste0(modeldir.ns, "/m_nbs_pred_rw_tw_50kn.RDS"))
saveRDS(pred.nbs.rw.tw.30kn, file = paste0(modeldir.ns, "/m_nbs_pred_rw_tw_30kn.RDS"))

# bring in saved predictions
pred.nbs.rw.tw.100kn <- readRDS(paste0(modeldir.ns, "/m_nbs_pred_rw_tw_100kn.RDS"))
pred.nbs.rw.tw.50kn <- readRDS(paste0(modeldir.ns, "/m_nbs_pred_rw_tw_50kn.RDS"))
pred.nbs.rw.tw.30kn <- readRDS(paste0(modeldir.ns, "/m_nbs_pred_rw_tw_30kn.RDS"))

# spatial biomass prediction plots ----
midpoint <- pred.nbs.rw.tw.100kn$data %>% mutate(est.exp = exp(est)) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.nbs.rw.tw.100kn <- ggplot(pred.nbs.rw.tw.100kn$data) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted abundance") +
  theme_gray() + 
  facet_wrap(~Year)+
  ggtitle("NSRKC RW, Tweedie, 100 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(plotdir.ns, "spatial_pred_rw_tw_100kn.png"), plot = sppred.nbs.rw.tw.100kn, height = 5, width = 7, units = "in")

midpoint <- pred.nbs.rw.tw.50kn$data %>% mutate(est.exp = exp(est)) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.nbs.rw.tw.50kn <- ggplot(pred.nbs.rw.tw.50kn$data) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted abundance") +
  theme_gray() + 
  facet_wrap(~Year)+
  ggtitle("NSRKC RW, Tweedie, 50 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(plotdir.ns, "spatial_pred_rw_tw_50kn.png"), plot = sppred.nbs.rw.tw.50kn, height = 5, width = 7, units = "in")


midpoint <- pred.nbs.rw.tw.30kn$data %>% mutate(est.exp = exp(est)) %>%
  summarise(est.exp.med = median(est.exp)) %>% pull(est.exp.med)
sppred.nbs.rw.tw.30kn <- ggplot(pred.nbs.rw.tw.30kn$data) +
  #geom_tile(aes(y = Lat, x = Lon, fill = est)) + 
  geom_point(aes(y = Lat, x = Lon, color = exp(est)), size = 1, alpha = 0.3) +
  scale_color_gradient2(low="#0072B2",  mid = "#0072B2", high="#F0E442", midpoint = midpoint) + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Predicted abundance") +
  theme_gray() + 
  facet_wrap(~Year)+
  ggtitle("NSRKC RW, Tweedie, 30 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(plotdir.ns, "spatial_pred_rw_tw_30kn.png"), plot = sppred.nbs.rw.tw.30kn, height = 5, width = 7, units = "in")


# *************************************************************************************************
# generate index ----
# *************************************************************************************************

index.nbs.rw.tw.100kn <- get_index(pred.nbs.rw.tw.100kn, area = pred.nbs.rw.tw.100kn$data$Area_km2, bias_correct = TRUE)
index.nbs.rw.tw.50kn <- get_index(pred.nbs.rw.tw.50kn, area = pred.nbs.rw.tw.50kn$data$Area_km2, bias_correct = TRUE)
index.nbs.rw.tw.30kn <- get_index(pred.nbs.rw.tw.30kn, area = pred.nbs.rw.tw.30kn$data$Area_km2, bias_correct = TRUE)

# index values in metric tons rather than kg
#index.nbs.rw.tw.30kn.t <- index.nbs.rw.tw.30kn %>% 
  #mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# save index values
write.csv(index.nbs.rw.tw.100kn, paste0(here::here(), "/NSRKC/output/nbs_index_rw_tw_100kn.csv"))
write.csv(index.nbs.rw.tw.50kn, paste0(here::here(), "/NSRKC/output/nbs_index_rw_tw_50kn.csv"))
write.csv(index.nbs.rw.tw.30kn, paste0(here::here(), "/NSRKC/output/nbs_index_rw_tw_30kn.csv"))

# plot index
ggplot(index.nbs.rw.tw.30kn, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.nbs.rw.tw.50kn, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.nbs.rw.tw.100kn, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

index.nbs.rw.tw <- rbind(
  index.nbs.rw.tw.30kn %>% mutate(model = "RW, Tweedie, 30 kn"),
  index.nbs.rw.tw.50kn %>% mutate(model = "RW, Tweedie, 50 kn"),
  index.nbs.rw.tw.100kn %>% mutate(model = "RW, Tweedie, 100 kn")
)

# ************************************************************************************************
# compare predicted index to observations ----
# ************************************************************************************************

obs.pred.nbs <- abund.crab.nbs %>%
  left_join(index.nbs.rw.tw) 

nbs.fit.rw.tw <- ggplot(obs.pred.nbs)+
  geom_line(aes(x = Year, y = est, color = model))+
  geom_ribbon(aes(x = Year, y = est, ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = Year, y = total.ab), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs.lwr95, ymax = obs.upr95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Red king crab abundance estimate (numbers)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") 

ggsave(file.path(plotdir.ns, "nsrkc_rw_tw_fit.png"), plot = nbs.fit.rw.tw, height = 5, width = 7, units = "in")


# ************************************************************************************************
# plot VAST index ----
# ************************************************************************************************

vast.index.nbs.iid.dg.30kn <- read.csv(paste0(here::here(), "/NSRKC/VAST output/NSRKC/NSRKC_MALE_GE64_30kts_V2025_Epsilon2_Omega1_Omega2_disabled/Index.csv")) %>%
  mutate(lwr = Estimate - 1.96*Std..Error.for.Estimate, upr = Estimate + 1.96*Std..Error.for.Estimate) %>%
  rename(Year = Time) %>%
  mutate(vast.est = Estimate * 1000, vast.lwr = lwr * 1000, vast.upr = upr * 1000)
vast.index.nbs.iid.dg.50kn <- read.csv(paste0(here::here(), "/NSRKC/VAST output/NSRKC/NSRKC_MALE_GE64_50kts_V2025_Epsilon2_Omega1_Omega2_disabled/Index.csv")) %>%
  mutate(lwr = Estimate - 1.96*Std..Error.for.Estimate, upr = Estimate + 1.96*Std..Error.for.Estimate) %>%
  rename(Year = Time) %>%
  mutate(vast.est = Estimate * 1000, vast.lwr = lwr * 1000, vast.upr = upr * 1000)
vast.index.nbs.iid.dg.100kn <- read.csv(paste0(here::here(), "/NSRKC/VAST output/NSRKC/NSRKC_MALE_GE64_100kts_V2025_Epsilon2_Omega1_Omega2_disabled/Index.csv")) %>%
  mutate(lwr = Estimate - 1.96*Std..Error.for.Estimate, upr = Estimate + 1.96*Std..Error.for.Estimate) %>%
  rename(Year = Time) %>%
  mutate(vast.est = Estimate * 1000, vast.lwr = lwr * 1000, vast.upr = upr * 1000)


ggplot(vast.index.nbs.iid.dg.30kn, aes(Year, vast.est)) + geom_line() +
  geom_ribbon(aes(ymin = vast.lwr, ymax = vast.upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)


vast.index.nbs.iid.dg <- rbind(
  vast.index.nbs.iid.dg.30kn %>% mutate(model = "IID, delta gamma, 30 kn"),
  vast.index.nbs.iid.dg.50kn %>% mutate(model = "IID, delta gamma, 50 kn"),
  vast.index.nbs.iid.dg.100kn %>% mutate(model = "IID, delta gamma, 100 kn")
)

vast.obs.pred.nbs <- abund.crab.nbs %>%
  left_join(vast.index.nbs.iid.dg) 


vast.nbs.fit.iid.dg <- ggplot(vast.obs.pred.nbs)+
  geom_line(aes(x = Year, y = vast.est, color = model))+
  geom_ribbon(aes(x = Year, y = vast.est, ymin = vast.lwr, ymax = vast.upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = Year, y = total.ab), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs.lwr95, ymax = obs.upr95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Red king crab abundance estimate (numbers)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  ggtitle("VAST")

ggsave(file.path(plotdir.ns, "vast_nsrkc_iid_dg_fit.png"), plot = vast.nbs.fit.iid.dg, height = 5, width = 7, units = "in")


# ************************************************************************************************
# compare sdmTMB predicted index to VAST index ----
# ************************************************************************************************

vast.sdmTMB.index.nbs <- rbind(
  vast.index.nbs.iid.dg.30kn %>% 
    mutate(model = "VAST, IID, delta gamma, 30 kn") %>%
    select(-c(lwr, upr)) %>%
    mutate(est = vast.est, lwr = vast.lwr, upr = vast.upr) %>%
    select(c(Year, est, lwr, upr, model)),
  index.nbs.rw.tw.30kn %>% 
    mutate(model = "sdmTMB, RW, Tweedie, 30 kn") %>%
    select(c(Year, est, lwr, upr, model))
)

vast.sdmTMB.obs.pred.nbs <- abund.crab.nbs %>%
  left_join(vast.sdmTMB.index.nbs) 


vast.sdmTMB.nbs.fit <- ggplot(vast.sdmTMB.obs.pred.nbs)+
  geom_line(aes(x = Year, y = est, color = model))+
  geom_ribbon(aes(x = Year, y = est, ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  geom_point(aes(x = Year, y = total.ab), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs.lwr95, ymax = obs.upr95), width = 0, color = "grey20")+
  xlab('Year') + 
  ylab('Red king crab abundance estimate (numbers)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom")

ggsave(file.path(plotdir.ns, "vast_sdmTMB_nsrkc_fit.png"), plot = vast.sdmTMB.nbs.fit, height = 5, width = 7, units = "in")


