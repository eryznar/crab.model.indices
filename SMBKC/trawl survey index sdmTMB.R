# ************************************************************************************************************
# Generating a spatiotemporal model-based index for St Matthew Island blue king crab in the NMFS trawl survey
# September 2024
# Caitlin Stern
# ************************************************************************************************************

# following the sdmTMB index standardization vignette: 
# https://pbs-assess.github.io/sdmTMB/articles/web_only/index-standardization.html

# load libraries
library(tidyverse)
library(sdmTMB)
library(ggplot2)

plotdir <- paste0(here::here(), "/SMBKC/plots")

# read in data
#SMBKC <- read.csv("adfg_smbkc_survey_cpuenum_VAST.csv")
haul_bkc <- read.csv(paste0(here::here(), '/SMBKC/data/trawl_survey/EBSCrab_Haul_', cur_yr, '.csv'), skip = 5)

haul_bkc %>% filter(SAMPLING_FACTOR >= 1 & is.na(LENGTH) == TRUE)

# prepare data
bkc_kgkm <- haul_bkc %>% filter(MID_LATITUDE > 58.5) %>% 
  dplyr::select(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, AREA_SWEPT, SPECIES_NAME, SEX, LENGTH, SAMPLING_FACTOR) %>% 
  filter(SEX == 1 & LENGTH >= 90) %>% 
  # calculate weight in kg. Mean weight by stage: 0.7 kg for Stage-1, 1.2 kg for Stage-2, and 1.9 kg for Stage-3
  mutate(wt.kg1 = case_when(
    LENGTH >= 90 & LENGTH < 104 ~ 0.7,
    LENGTH >= 104 & LENGTH < 120 ~ 1.2,
    LENGTH >= 120 ~ 1.9
  )) %>%
  mutate(wt.kg = wt.kg1 * SAMPLING_FACTOR) %>%
  group_by(AKFIN_SURVEY_YEAR, GIS_STATION) %>% 
  mutate(total.wt.kg = sum(wt.kg)) %>% 
  mutate(AREA_SWEPT_km2 = AREA_SWEPT/0.29155335) %>% # convert from square nautical miles to square km
  mutate(kg.km = total.wt.kg / AREA_SWEPT_km2) %>%
  select(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, AREA_SWEPT, SPECIES_NAME, total.wt.kg, kg.km) %>%
  distinct() %>%
  arrange(AKFIN_SURVEY_YEAR, GIS_STATION) %>%
  rename("SURVEY_YEAR" = "AKFIN_SURVEY_YEAR", "LATITUDE" = "MID_LATITUDE", "LONGITUDE" = "MID_LONGITUDE")

# make SPDE mesh
#BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), cutoff="10") 
BK_spde <- make_mesh(bkc_kgkm, c("LONGITUDE","LATITUDE"), n_knots = 100) 
plot(BK_spde)
BK_spde$mesh$n

# fit a GLMM
m <- sdmTMB(
  data = bkc_kgkm, 
  formula = kg.km ~ 0 + as.factor(SURVEY_YEAR), #the 0 is there so there is a factor predictor for each time slice
  time = "SURVEY_YEAR", mesh = BK_spde, family = tweedie(link = "log"))

# inspect randomized quantile residuals
bkc_kgkm$resids <- residuals(m) # randomized quantile residuals

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

# now predict a fine scale grid over entire domain
grid_new <- expand.grid(LONGITUDE = seq(from= min(bkc_kgkm$LONGITUDE), to= max(bkc_kgkm$LONGITUDE), by=0.05),
                       LATITUDE = seq(from= min(bkc_kgkm$LATITUDE), to= max(bkc_kgkm$LATITUDE), by=0.05)
) 

# look at predicted grid
plot(grid_new) 

# replicate grid across all years
grid_yrs <- replicate_df(grid_new, "SURVEY_YEAR", unique(bkc_kgkm$SURVEY_YEAR))
dplyr::glimpse(grid_yrs)

# predictions on new data
predictions <- predict(m, newdata = grid_yrs, return_tmb_object = T)

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
