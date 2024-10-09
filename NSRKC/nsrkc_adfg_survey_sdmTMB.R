# ************************************************************************************************************
# Generating a spatiotemporal model-based index for Norton Sound red king crab in the NBS trawl survey
# September 2024
# Caitlin Stern
# ************************************************************************************************************

# following the sdmTMB index standardization vignette: 
# https://pbs-assess.github.io/sdmTMB/articles/web_only/index-standardization.html

# load libraries
library(tidyverse)
library(sdmTMB)
library(ggplot2)

plotdir <- paste0(here::here(), "/NSRKC/plots")

# read in data
haul_rkc <- read.csv(paste0(here::here(), '/NSRKC/data/nmfs_survey/NBSCRAB - Haul Report.csv'))

# basic QAQC
haul_rkc %>% filter(SAMPLING_FACTOR >= 1 & is.na(WIDTH_MM) == TRUE)

haul_rkc %>% filter(SAMPLING_FACTOR > 1)

# prepare data
rkc_kgkm <- haul_rkc %>%
  dplyr::select(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, AREA_SWEPT, SPECIES_NAME, SEX, WIDTH_MM, SAMPLING_FACTOR, CALCULATED_WEIGHT) %>% 
  filter(SEX == 1 & WIDTH_MM >= 64) %>% 
  mutate(wt.kg = CALCULATED_WEIGHT / 1000) %>%
  group_by(AKFIN_SURVEY_YEAR, GIS_STATION) %>% 
  mutate(total.wt.kg = sum(wt.kg)) %>% 
  mutate(AREA_SWEPT_km2 = AREA_SWEPT/0.29155335) %>% # convert from square nautical miles to square km
  mutate(kg.km = total.wt.kg / AREA_SWEPT_km2) %>%
  select(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, AREA_SWEPT, SPECIES_NAME, total.wt.kg, kg.km) %>%
  distinct() %>%
  arrange(AKFIN_SURVEY_YEAR, GIS_STATION) %>%
  rename("SURVEY_YEAR" = "AKFIN_SURVEY_YEAR", "LATITUDE" = "MID_LATITUDE", "LONGITUDE" = "MID_LONGITUDE")


# make SPDE mesh
RK_spde <- make_mesh(rkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 80) 
plot(RK_spde)
RK_spde$mesh$n

# fit a GLMM
m.nsrkc <- sdmTMB(
  data = rkc_kgkm, 
  formula = kg.km ~ 0 + as.factor(SURVEY_YEAR), #the 0 is there so there is a factor predictor for each time slice
  time = "SURVEY_YEAR", mesh = RK_spde, family = tweedie(link = "log"))
#Warning message:
#The model may not have converged. Maximum final gradient: 0.0363312230861308. 

# print the model fit
m.nsrkc

# view parameters
tidynsrkc <- tidy(m.nsrkc, conf.int = TRUE)
tidy.nsrkc.ran <- tidy(m.nsrkc, effects = "ran_pars", conf.int = TRUE)

# run sanity checks
sanity(m.nsrkc)

# inspect randomized quantile residuals
rkc_kgkm$resids <- residuals(m.nsrkc) # randomized quantile residuals

#> Note what used to be the default sdmTMB residuals (before version 0.4.3.9005)
#> are now `type = 'mle-eb'`. We recommend using the current default `'mle-mvn'`,
#> which takes one sample from the approximate posterior of the random effects or
#> `dharma_residuals()` using a similar approach.
#> 
#> # Also see residuals(..., type = "mle-mcmc") which are better but slower

hist(rkc_kgkm$resids)
qqnorm(rkc_kgkm$resids)
abline(a=0, b=1)


ggplot(rkc_kgkm, aes(LONGITUDE, LATITUDE, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~SURVEY_YEAR) + coord_fixed()

# now predict a fine scale grid over entire domain
grid_new <- expand.grid(LONGITUDE = seq(from= min(rkc_kgkm$LONGITUDE), to= max(rkc_kgkm$LONGITUDE), by=0.05),
                        LATITUDE = seq(from= min(rkc_kgkm$LATITUDE), to= max(rkc_kgkm$LATITUDE), by=0.05)
) 

# look at predicted grid
plot(grid_new) 

# replicate grid across all years
grid_yrs <- replicate_df(grid_new, "SURVEY_YEAR", unique(rkc_kgkm$SURVEY_YEAR))
dplyr::glimpse(grid_yrs)

# predictions on new data
predictions <- predict(m.nsrkc, newdata = grid_yrs, return_tmb_object = T)

# make a map function
plot_map <- function(dat, column) {
  ggplot(dat, aes(LONGITUDE, LATITUDE, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~SURVEY_YEAR) +
    coord_fixed()
}

# types of predictions
## predictions that incorporate all fixef and ranef
map.fe.re <- plot_map(predictions$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

map.fe.re.yrs <- plot_map(predictions$data %>% filter(SURVEY_YEAR > 2003), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") + #sqrt so we can see changes better?
  ggtitle("Prediction (fixed effects + all random effects)")

ggsave(file.path(plotdir, "map_fixedeff_randomeff.png"), plot = map.fe.re, height = 4.2, width = 7, units = "in")

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

#use get_index() function to extract the total biomass calculations and standard errors
# not sure what argument to use for area
index <- get_index(predictions, area = 1, bias_correct = TRUE)

# index values in metric tons rather than kg
index.t <- index %>% mutate(est.t = est / 1000, lwr.t = lwr/1000, upr.t = upr/1000)

# plot index
ggplot(index.t, aes(SURVEY_YEAR, est.t)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.t, ymax = upr.t), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (t)') +
  scale_y_continuous(label=scales::comma)

# biomass estimates
mutate(index, cv = sqrt(exp(se^2) - 1)) %>% 
  select(-log_est, -se) %>%
  knitr::kable(format = "pandoc", digits = c(0, 0, 0, 0, 2))

########################################################
