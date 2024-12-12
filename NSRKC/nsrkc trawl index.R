# ************************************************************************************************************
# Generating a spatiotemporal model-based index for Norton Sound red king crab in trawl surveys
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

# **************************************************************************************************************
# read in and process survey data ----
# **************************************************************************************************************

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

# number of unique stations
length(unique(nsrkc_utm$ADFG_Station))

# check number of observations per year (this should be larger than the number of vertices in the mesh)
nsrkc_utm %>%
  group_by(Year) %>%
  count() %>%
  print(n = 50)

# plot survey data
ggplot(nsrkc_utm) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

ggplot(nsrkc_utm) + 
  geom_point(aes(x = Longitude, y = Latitude)) +
  theme_bw()

# **************************************************************************************************************
# make SPDE mesh ----
# **************************************************************************************************************

NS_spde_100kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 100, type = "kmeans")

plot(NS_spde_100kn)
NS_spde_100kn$mesh$n

NS_spde_70kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 70, type = "kmeans")

plot(NS_spde_70kn)
NS_spde_70kn$mesh$n

NS_spde_60kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 60, type = "kmeans")

plot(NS_spde_60kn)
NS_spde_60kn$mesh$n

NS_spde_50kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 50, type = "kmeans")

plot(NS_spde_50kn)
NS_spde_50kn$mesh$n

NS_spde_40kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 40, type = "kmeans")

plot(NS_spde_40kn)
NS_spde_40kn$mesh$n


NS_spde_30kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 30, type = "kmeans")

plot(NS_spde_30kn)
NS_spde_30kn$mesh$n

# plot mesh with data points

mesh.plot.100kn.ns <- ggplot(nsrkc_utm) + 
  inlabru::gg(NS_spde_100kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

mesh.plot.50kn.ns <- ggplot(nsrkc_utm) + 
  inlabru::gg(NS_spde_50kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

mesh.plot.30kn.ns <- ggplot(nsrkc_utm) + 
  inlabru::gg(NS_spde_30kn$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

# export mesh plots

ggsave(file.path(plotdir.ns, "mesh_100kn_ns.png"), plot = mesh.plot.100kn.ns, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "mesh_50kn_ns.png"), plot = mesh.plot.50kn.ns, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "mesh_30kn_ns.png"), plot = mesh.plot.30kn.ns, height = 4.2, width = 7, units = "in")

# **************************************************************************************************************
# fit and check models ----
# **************************************************************************************************************

# fit a GLMM with spatiotemporal = IID and 100 knots
m.nsrkc.iid.100kn <- sdmTMB(
  data = nsrkc_utm, 
  formula = crab.km2 ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "Year", 
  mesh = NS_spde_100kn, 
  spatiotemporal = "iid",
  #extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 50 knots
m.nsrkc.iid.50kn <- sdmTMB(
  data = nsrkc_utm, 
  formula = crab.km2 ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "Year", 
  mesh = NS_spde_50kn, 
  spatiotemporal = "iid",
  #extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))

# fit a GLMM with spatiotemporal = IID and 30 knots
m.nsrkc.iid.30kn <- sdmTMB(
  data = nsrkc_utm, 
  formula = crab.km2 ~ 0 + year_f, #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "Year", 
  mesh = NS_spde_30kn, 
  spatiotemporal = "iid",
  #extra_time = c(2020),
  silent = FALSE,
  anisotropy = TRUE,
  family = tweedie(link = "log"))


# save the fitted models
saveRDS(m.nsrkc.iid.100kn, file = paste0(modeldir.ns, "/m_nsrkc_iid_100kn.RDS"))
saveRDS(m.nsrkc.iid.50kn, file = paste0(modeldir.ns, "/m_nsrkc_iid_50kn.RDS"))
saveRDS(m.nsrkc.iid.30kn, file = paste0(modeldir.ns, "/m_nsrkc_iid_30kn.RDS"))

# print the model fit
m.nsrkc.iid.100kn
m.nsrkc.iid.100kn$sd_report

m.nsrkc.iid.50kn
m.nsrkc.iid.50kn$sd_report

m.nsrkc.iid.30kn
m.nsrkc.iid.30kn$sd_report

# view parameters
tidy.m.nsrkc.iid.50kn <- tidy(m.nsrkc.iid.50kn, conf.int = TRUE)
tidy.m.nsrkc.iid.50kn <- tidy(m.nsrkc.iid.50kn, effects = "ran_pars", conf.int = TRUE)

# run sanity checks
sanity(m.nsrkc.iid.100kn)
sanity(m.nsrkc.iid.50kn)
sanity(m.nsrkc.iid.30kn)



# **************************************************************************************************************
# model residuals ----
# **************************************************************************************************************

# DHARMa residuals
m.nsrkc.iid.100kn.resid <- simulate(update(m.nsrkc.iid.100kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.iid.100kn, return_DHARMa = TRUE)
plot(m.nsrkc.iid.100kn.resid)

m.nsrkc.iid.50kn.resid <- simulate(update(m.nsrkc.iid.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.iid.50kn, return_DHARMa = TRUE)
plot(m.nsrkc.iid.50kn.resid)

m.nsrkc.iid.30kn.resid <- simulate(update(m.nsrkc.iid.30kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.iid.30kn, return_DHARMa = TRUE)
plot(m.nsrkc.iid.30kn.resid)

DHARMa::testOutliers(m.nsrkc.iid.30kn.resid)
DHARMa::testOutliers(m.nsrkc.iid.50kn.resid)
DHARMa::testOutliers(m.nsrkc.iid.100kn.resid)

DHARMa::testQuantiles(m.nsrkc.iid.30kn.resid)
DHARMa::testQuantiles(m.nsrkc.iid.50kn.resid)
DHARMa::testQuantiles(m.nsrkc.iid.100kn.resid)

DHARMa::testDispersion(m.nsrkc.iid.30kn.resid)
DHARMa::testDispersion(m.nsrkc.iid.50kn.resid)
DHARMa::testDispersion(m.nsrkc.iid.100kn.resid)

DHARMa::testResiduals(m.nsrkc.iid.30kn.resid)
DHARMa::testResiduals(m.nsrkc.iid.50kn.resid)
DHARMa::testResiduals(m.nsrkc.iid.100kn.resid)

DHARMa::testSpatialAutocorrelation(m.nsrkc.iid.90kn.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)
DHARMa::testSpatialAutocorrelation(m.nsrkc.iid.120kn.dgs.resid, x = bkc_kgkm_utm$X, y = bkc_kgkm_utm$Y)

DHARMa::testZeroInflation(m.nsrkc.iid.30kn.resid)
DHARMa::testZeroInflation(m.nsrkc.iid.50kn.resid)
DHARMa::testZeroInflation(m.nsrkc.iid.100kn.resid)


# **************************************************************************************************************
# generate predictions ----
# **************************************************************************************************************

# replicate grid across all years
pred_grid.ns <- replicate_df(predgrid_utm.ns, "Year", unique(nsrkc_utm$Year))
pred_grid.ns$year_f <- factor(pred_grid.ns$Year)
dplyr::glimpse(pred_grid.ns)

# predictions on new data
predictions.iid.100kn.ns <- predict(m.nsrkc.iid.100kn, newdata = pred_grid.ns, return_tmb_object = T)
predictions.iid.50kn.ns <- predict(m.nsrkc.iid.50kn, newdata = pred_grid.ns, return_tmb_object = T)
predictions.iid.30kn.ns <- predict(m.nsrkc.iid.30kn, newdata = pred_grid.ns, return_tmb_object = T)

# save the predictions
write.csv(predictions.iid.100kn.ns$data, paste0(outdir, "/m_nsrkc_predictions_iid_100kn.csv"))
write.csv(predictions.iid.50kn.ns$data, paste0(outdir, "/m_nsrkc_predictions_iid_50kn.csv"))
write.csv(predictions.iid.30kn.ns$data, paste0(outdir, "/m_nsrkc_predictions_iid_30kn.csv"))

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

# heat map of predictions
pred.heat.iid.100kn.ns <- ggplot(predictions.iid.100kn.ns$data, aes(y = Lat, x = Lon, color = exp(est))) +
  #geom_tile(aes(y = Lat, x = Lon, color = exp(est1))) + 
  geom_point(size = 1, alpha = 0.5) +
  #scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Est RKC per sq.km") +
  #theme_gray() + 
  facet_wrap(~year_f)+
  ggtitle("NSRKC predicted abundance, 100 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom", axis.text.x = element_text(size = 8, angle = 90)) +
  scale_color_viridis_c(na.value = "white", limits = c(0, quantile(exp(predictions.iid.100kn.ns$data$est), 0.995)))

pred.heat.iid.50kn.ns <- ggplot(predictions.iid.50kn.ns$data, aes(y = Lat, x = Lon, color = exp(est))) +
  #geom_tile(aes(y = Lat, x = Lon, color = exp(est1))) + 
  geom_point(size = 1, alpha = 0.5) +
  #scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Est RKC per sq.km") +
  #theme_gray() + 
  facet_wrap(~year_f)+
  ggtitle("NSRKC predicted abundance, 50 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom", axis.text.x = element_text(size = 8, angle = 90)) +
  scale_color_viridis_c(na.value = "white", limits = c(0, quantile(exp(predictions.iid.50kn.ns$data$est), 0.995)))


pred.heat.iid.30kn.ns <- ggplot(predictions.iid.30kn.ns$data, aes(y = Lat, x = Lon, color = exp(est))) +
  #geom_tile(aes(y = Lat, x = Lon, color = exp(est1))) + 
  geom_point(size = 1, alpha = 0.5) +
  #scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       color = "Est RKC per sq.km") +
  #theme_gray() + 
  facet_wrap(~year_f)+
  ggtitle("NSRKC predicted abundance, 30 knots")+
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom", axis.text.x = element_text(size = 8, angle = 90)) +
  scale_color_viridis_c(na.value = "white", limits = c(0, quantile(exp(predictions.iid.30kn.ns$data$est), 0.995)))


ggsave(file.path(plotdir.ns, "pred_heat_iid_100kn.png"), plot = pred.heat.iid.100kn.ns, height = 7, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "pred_heat_iid_50kn.png"), plot = pred.heat.iid.50kn.ns, height = 7, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "pred_heat_iid_30kn.png"), plot = pred.heat.iid.30kn.ns, height = 7, width = 7, units = "in")

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

index.iid.100kn.ns <- get_index(predictions.iid.100kn.ns, area = predictions.iid.100kn.ns$data$Area_km2, bias_correct = TRUE)
index.iid.50kn.ns <- get_index(predictions.iid.50kn.ns, area = predictions.iid.50kn.ns$data$Area_km2, bias_correct = TRUE)
index.iid.30kn.ns <- get_index(predictions.iid.30kn.ns, area = predictions.iid.30kn.ns$data$Area_km2, bias_correct = TRUE)

# save index values
write.csv(index.iid.100kn.ns, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_iid_100kn.csv"))
write.csv(index.iid.50kn.ns, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_iid_50kn.csv"))
write.csv(index.iid.30kn.ns, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_iid_30kn.csv"))

# combine indices to compare models
mod.compare.iid.ns <- rbind(index.iid.100kn.ns %>% 
                           mutate(model = "100 knots"), 
                         index.iid.50kn.ns 
                         %>% mutate(model = "50 knots")) %>%
  rbind(index.iid.30kn.ns %>% mutate(model = "30 knots"))

# plot index
ggplot(mod.compare.iid.ns, aes(Year, est, group = model)) + 
  geom_line(aes(Year, est, color = model)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  xlab('Year') + ylab('Abundance estimate (t)') +
  scale_y_continuous(label=scales::comma)

# abundance estimates
mutate(index.iid.50kn.ns, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))


# **************************************************************************************************************
# compare predicted index to observations ----
# **************************************************************************************************************

# Abundance and CV of the NMFS 1976-1991 surveys were provided by NOAA (Jon Richer NOAA personal communication).  
# The abundance was estimated by averaging catch CPUE (#/nm2) of all stations (including survey stations out of Norton Sound)
# that was multiplied by standard Norton Sound Area (7600 nm2) (i.e., N = 7600*mean CPUE). The ADF&G survey abundance is 
# calculated at each station (i.e., n= CPUE*100 nm2) and summed across all surveyed stations (i.e., N = sum of 100*CPUEs) 
# (Bell and Hamazaki 2019).  Extent of the ADF&G survey coverage differed among years due to survey conditions, and survey 
# abundance has not been standardized.  NOAA NBS survey abundance is estimated by the author in similar manner as ADF&G 
# survey with the data limited to the Norton Sound survey area that overlaps the ADF&G survey area (5841 nm2) 

# abundance for 1976-1991 NMFS surveys
abund.crab.noaa <- nsrkc_utm %>%
  filter(Agent == "NOAA") %>%
  group_by(Agent, Year) %>%
  summarise(mean.crab.km2 = mean(crab.km2), 
            sd.crab.km2 = sd(crab.km2, na.rm = TRUE),
            n.station = length(ADFG_Station)) %>%
  mutate(se.crab.km2 = sd.crab.km2 / sqrt(n.station)) %>%
  mutate(total.ab = 7600 * 3.4299 * mean.crab.km2, 
         se.total.ab = 7600 * 3.4299 * se.crab.km2) %>%
  mutate(obs.lwr95 = total.ab - 1.96*se.total.ab, 
         obs.upr95 = total.ab + 1.96*se.total.ab)

# abundance for ADF&G survey
abund.crab.adfg <- nsrkc_utm %>%
  filter(Agent == "ADFG") %>%
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

# abundance for NBS survey
abund.crab.nbs <- nsrkc_utm %>%
  filter(Agent == "NBS") %>%
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

# total abundance
abund.crab <- rbind(abund.crab.noaa, abund.crab.adfg, abund.crab.nbs)

survey.plot.ns <- ggplot(abund.crab %>% rename("Survey" = "Agent"))+
  geom_point(aes(x = Year, y = total.ab, color = Survey))+
  geom_errorbar(aes(x = Year, ymin = obs.lwr95, ymax = obs.upr95, color = Survey), width = 0)+
  xlab('Year') + 
  ylab('Abundance estimate (number of crab)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.text = element_text(size = 12)) 

ggsave(file.path(plotdir.ns, "nsrkc_survey_abundance.png"), plot = survey.plot.ns, height = 5, width = 7, units = "in")

write.csv(abund.crab, paste0(here::here(), "/output/NSRKC_trawl_survey_total_abundance.csv"))

# **************************************************************************************************************
# compare predicted index to observations ----
# **************************************************************************************************************

obs.pred.ns <- left_join(mod.compare.iid.ns, abund.crab, by = "Year")

obs.pred.plot.ns <- ggplot(obs.pred.ns)+
  geom_line(aes(x = Year, y = est, color = model))+
  geom_ribbon(aes(x = Year, y = est, ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  #geom_point(aes(x = Year, y = total.ab, color = Agent))+
  #geom_errorbar(aes(x = Year, ymin = obs.lwr95, ymax = obs.upr95, color = Agent), width = 0)+
  geom_point(aes(x = Year, y = total.ab), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs.lwr95, ymax = obs.upr95), color = "grey20", width = 0)+
  xlab('Year') + 
  ylab('Abundance estimate (number of crab)') +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, NA))+
  scale_color_manual(values = cbpalette) +
  scale_fill_manual(values = cbpalette) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.text = element_text(size = 12)) 

ggsave(file.path(plotdir.ns, "nsrkc_model_fit.png"), plot = obs.pred.plot.ns, height = 5, width = 7, units = "in")

# ******************************************************************************