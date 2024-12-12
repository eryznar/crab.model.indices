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


# **************************************************************************************************************
# read in and process survey data ----
# **************************************************************************************************************

nsrkc.survey <- read.csv(paste0(here::here(), "/NSRKC/data/NSRKC_trawl_survey_abundance.csv"))

nsrkc.dt <- nsrkc.survey %>%
  select(Year, ADFG_Station, Latitude, Longitude, crab.km2) %>%
  arrange(Year, ADFG_Station) %>%
  mutate("year_f" = factor(Year))

nsrkc_utm_zone <- nsrkc.dt %>%  
  # get UTM zone
  mutate(zone = (floor((Longitude + 180)/6) %% 60) + 1) 

unique(nsrkc_utm_zone$zone)

nsrkc_utm_z2 <- nsrkc_utm_zone %>% filter(zone == 2)
nsrkc_utm_z3 <- nsrkc_utm_zone %>% filter(zone == 3)
nsrkc_utm_z4 <- nsrkc_utm_zone %>% filter(zone == 4)
nsrkc_utm_z13 <- nsrkc_utm_zone %>% filter(zone == 13)

d2 <- nsrkc_utm_z2
get_crs(d2, c("Longitude", "Latitude"))
d2u <- add_utm_columns(d2, c("Longitude", "Latitude"))

d3 <- nsrkc_utm_z3
get_crs(d3, c("Longitude", "Latitude"))
d3u <- add_utm_columns(d3, c("Longitude", "Latitude"))

d4 <- nsrkc_utm_z4
get_crs(d4, c("Longitude", "Latitude"))
d4u <- add_utm_columns(d4, c("Longitude", "Latitude"))

d13 <- nsrkc_utm_z13
get_crs(d13, c("Longitude", "Latitude"))
d13u <- add_utm_columns(d13, c("Longitude", "Latitude"))

nsrkc_utm <- rbind(d2u, d3u, d4u) %>%
  mutate(X = X / 1000, Y = Y / 1000) # convert UTM coordinates from meter to kilometers to reduce computing needs

# number of unique stations
length(unique(nsrkc_utm$ADFG_Station))

# check number of observations per year (this should be larger than the number of vertices in the mesh)
nsrkc_utm %>%
  group_by(Year) %>%
  count() %>%
  print(n = 50)




# **************************************************************************************************************
# make SPDE mesh ----
# **************************************************************************************************************

NS_spde_50kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 50, type = "kmeans")

plot(NS_spde_50kn)
NS_spde_50kn$mesh$n

NS_spde_40kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 40, type = "kmeans")

plot(NS_spde_40kn)
NS_spde_40kn$mesh$n


NS_spde_30kn <- make_mesh(nsrkc_utm, xy_cols = c("X","Y"), n_knots = 30, type = "kmeans")

plot(NS_spde_30kn)
NS_spde_30kn$mesh$n

# **************************************************************************************************************
# fit and check models ----
# **************************************************************************************************************

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



# save the fitted models
saveRDS(m.nsrkc.iid.50kn, file = paste0(modeldir.ns, "/m_nsrkc_iid_2024_12_02.RDS"))


# print the model fit
m.nsrkc.iid.50kn
m.nsrkc.iid.50kn$sd_report

# view parameters
tidy.nsrkc.ar1 <- tidy(m.nsrkc.ar1, conf.int = TRUE)
tidy.nsrkc.ran.ar1 <- tidy(m.nsrkc.ar1 , effects = "ran_pars", conf.int = TRUE)

# run sanity checks
sanity(m.nsrkc.iid.50kn)



# **************************************************************************************************************
# model residuals ----
# **************************************************************************************************************

# DHARMa residuals
m.nsrkc.iid.50kn.resid <- simulate(update(m.nsrkc.iid.50kn, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nsrkc.iid.50kn, return_DHARMa = TRUE)
plot(m.nsrkc.iid.resid.50kn)





# **************************************************************************************************************
# generate predictions ----
# **************************************************************************************************************

# replicate grid across all years
pred_grid.ns <- replicate_df(predgrid_utm.ns, "Year", unique(nsrkc_utm$Year))
pred_grid.ns$year_f <- factor(pred_grid.ns$Year)
dplyr::glimpse(pred_grid.ns)

# predictions on new data
predictions.iid.50kn.ns <- predict(m.nsrkc.iid.50kn, newdata = pred_grid.ns, return_tmb_object = T)

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

index.iid.50kn.ns <- get_index(predictions.iid.50kn.ns, area = predictions.iid.50kn.ns$data$Area_km2, bias_correct = TRUE)

# index values in metric tons rather than kg
index.iid.50kn.ns.t <- index.iid.50kn.ns %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# save index values
write.csv(index.ar1.ns.t, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_ar1_2024_12_02.csv"))
write.csv(index.50kn.ns.t, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_50kn_2024_12_02.csv"))
write.csv(index.iid.ns.t, paste0(here::here(), "/NSRKC/output/nsrkc_nmfs_trawl_index_iid_2024_12_02.csv"))

# plot index
ggplot(index.iid.50kn.ns.t, aes(Year, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (crab / km2)') +
  scale_y_continuous(label=scales::comma)

# abundance estimates
mutate(index.iid.50kn.ns.t, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))


# **************************************************************************************************************
# compare predicted index to observations ----
# **************************************************************************************************************


# ******************************************************************************