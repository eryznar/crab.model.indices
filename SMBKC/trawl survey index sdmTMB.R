# ************************************************************************************************************
# Generating a spatiotemporal model-based index for St Matthew Island blue king crab in the NMFS trawl survey
# September 2024
# Caitlin Stern
# ************************************************************************************************************

# following the sdmTMB index standardization vignette: 
# https://pbs-assess.github.io/sdmTMB/articles/web_only/index-standardization.html


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
#cbpalette <- c(cbpalette, "tomato3", "turquoise4", "orangered4")

# **************************************************************************************************************
# read in and process data
# **************************************************************************************************************

haul_bkc <- read.csv(paste0(here::here(), '/SMBKC/data/trawl_survey/EBSCrab_Haul_', cur_yr, '.csv'), skip = 5)

# error check
haul_bkc %>% filter(SAMPLING_FACTOR >= 1 & is.na(LENGTH) == TRUE)

# read in grid from Jon Richar to use for predictions and convert to UTM
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

predgrid_utm <- rbind(u1u, u2u)

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
  rename("SURVEY_YEAR" = "AKFIN_SURVEY_YEAR", "LATITUDE" = "MID_LATITUDE", "LONGITUDE" = "MID_LONGITUDE")
  
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

bkc_kgkm_utm <- rbind(d1u, d2u, d3u)

# **************************************************************************************************************
# make SPDE mesh ----
# **************************************************************************************************************

#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), cutoff="10") 
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 100)
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 50)
BK_spde <- make_mesh(bkc_kgkm_utm, c("X","Y"), n_knots = 120)
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 120)
plot(BK_spde)
BK_spde$mesh$n

#BK_spde.2 <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 10)

# look at spatial range of the model (Matern range - want 2 knots per range distance)


# **************************************************************************************************************
# fit and check model ----
# **************************************************************************************************************

# fit a GLMM
m.smbkc <- sdmTMB(
  data = bkc_kgkm, 
  formula = kg.km ~ 0 + as.factor(SURVEY_YEAR), #the 0 is there so there is a factor predictor for each time slice
  time = "SURVEY_YEAR", mesh = BK_spde, family = tweedie(link = "log"))

m.smbkc.utm <- sdmTMB(
  data = bkc_kgkm_utm, 
  formula = kg.km ~ 0 + as.factor(SURVEY_YEAR), #the 0 is there so there is a factor predictor for each time slice
  spatial = "on",
  time = "SURVEY_YEAR", 
  mesh = BK_spde, 
  spatiotemporal = "IID",
  family = tweedie(link = "log"))

# barrier effect? To take into account the island
# size comps - could just fit 3 different models

# print the model fit
m.smbkc
m.smbkc$sd_report

m.smbkc.utm
m.smbkc.utm$sd_report

# view parameters
tidy.smbkc <- tidy(m.smbkc, conf.int = TRUE)
tidy.smbkc.ran <- tidy(m.smbkc, effects = "ran_pars", conf.int = TRUE)

tidy.smbkc.utm <- tidy(m.smbkc.utm, conf.int = TRUE)
tidy.smbkc.ran.utm <- tidy(m.smbkc.utm , effects = "ran_pars", conf.int = TRUE)

# run sanity checks
sanity(m.smbkc)
sanity(m.smbkc.utm)

# inspect randomized quantile residuals
bkc_kgkm$resids <- residuals(m.smbkc) # randomized quantile residuals
bkc_kgkm_utm$resids <- residuals(m.smbkc.utm) # randomized quantile residuals
bkc_kgkm_utm$dharma_resids <- dharma_residuals(m.smbkc.utm)


fit_dl <- update(m.smbkc.utm, family = delta_lognormal())
simulate(fit_dl, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.utm)

ret <- simulate(m.smbkc.utm.dl, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.smbkc.utm, return_DHARMa = TRUE)
plot(ret)

#> Note what used to be the default sdmTMB residuals (before version 0.4.3.9005)
#> are now `type = 'mle-eb'`. We recommend using the current default `'mle-mvn'`,
#> which takes one sample from the approximate posterior of the random effects or
#> `dharma_residuals()` using a similar approach.
#> 
#> # Also see residuals(..., type = "mle-mcmc") which are better but slower

hist(bkc_kgkm$resids)
qqnorm(bkc_kgkm$resids)
abline(a=0, b=1)


ggplot(bkc_kgkm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()


ggplot(bkc_kgkm_utm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()

# **************************************************************************************************************
# generate predictions
# **************************************************************************************************************

# now predict a fine scale grid over entire domain
grid_new <- expand.grid(LONGITUDE = seq(from= min(bkc_kgkm$LONGITUDE), to= max(bkc_kgkm$LONGITUDE), by=0.05),
                       LATITUDE = seq(from= min(bkc_kgkm$LATITUDE), to= max(bkc_kgkm$LATITUDE), by=0.05)
) 

grid_new_utm <- expand.grid(X = seq(from= min(bkc_kgkm_utm$X), to= max(bkc_kgkm_utm$X), by=2),
                        Y = seq(from= min(bkc_kgkm_utm$Y), to= max(bkc_kgkm_utm$Y), by=2)
) 

grid_noland_4km <- expand.grid(LONGITUDE = grid4km$Lon, LATITUDE = grid4km$Lat)

grid_noland_4km_utm <- expand.grid(X = predgrid_utm$X, Y = predgrid_utm$Y)

# look at predicted grid
plot(grid_new) 
plot(grid_noland_4km)
plot(grid_new_utm) 

# replicate grid across all years
grid_yrs <- replicate_df(grid_new, "SURVEY_YEAR", unique(bkc_kgkm$SURVEY_YEAR))
dplyr::glimpse(grid_yrs)

grid_yrs_noland4km_utm <- replicate_df(predgrid_utm, "SURVEY_YEAR", unique(bkc_kgkm$SURVEY_YEAR))
dplyr::glimpse(grid_yrs_noland4km_utm)

grid_yrs_utm <- replicate_df(grid_new_utm, "SURVEY_YEAR", unique(bkc_kgkm_utm$SURVEY_YEAR))
dplyr::glimpse(grid_yrs_utm)

# predictions on new data
predictions <- predict(m.smbkc, newdata = grid_yrs, return_tmb_object = T)

predictions_noland4km_utm <- predict(m.smbkc.utm, newdata = grid_yrs_noland4km_utm, return_tmb_object = T)

predictions_utm <- predict(m.smbkc.utm, newdata = grid_yrs_utm, return_tmb_object = T)

# make a map function
plot_map <- function(dat, column) {
  ggplot(dat, aes(LONGITUDE, LATITUDE, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~SURVEY_YEAR) +
    coord_fixed()
}

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
# generate index
# **************************************************************************************************************

#use get_index() function to extract the total biomass calculations and standard errors
# not sure what argument to use for area. Area should be equal to area of grid cells that you're predicting over.
# check grid set up after converting to UTM
# when CPUE is catch/km^2
# can change aggregation function - if response variable is density, makes sense to sum. But could make more sense
# to take average over the area. 
# UTM is in meters

index <- get_index(predictions, area = 1, bias_correct = TRUE)

index_noland4km <- get_index(predictions_noland4km_utm, area = 16, bias_correct = TRUE)

index_utm <- get_index(predictions_utm, area = 10, bias_correct = TRUE)

index_utm4 <- get_index(predictions_utm, area = 4, bias_correct = TRUE)

index_utm4.t <- index_utm4 %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# index values in metric tons rather than kg
index.t <- index_noland4km %>% 
  mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000) %>%
  mutate(est.t2 = est / 100, lwr.t2 = lwr/100, upr.t2 = upr/100)

# plot index
ggplot(index.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(index_utm4.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

ggplot(index_noland4km, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

# biomass estimates
mutate(index_noland4km, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))

# plot predicted vs. observed index

obs.biom <- read.csv(paste0(here::here(), '/SMBKC/data/trawl_survey/survey_biomass_mt2.csv'))

obs.pred <- left_join(index.t, obs.biom, by = "SURVEY_YEAR") %>% filter(SURVEY_YEAR >= 1978) %>%
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

ggplot(index.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

########################################################
