# ************************************************************************************************
# Generating a spatiotemporal model-based index for Norton Sound red king crab 
# in each of the three surveys 
# March 2025
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
library(scales)
library(fmesher)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(terra)
library(ggpubr)
library(grid)
library(remotes)
library(sdmTMBextra)
library(units)
library(spdep)
library(stars)

# install.packages("remotes")
#remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)

cur_yr.ns <- 2024

plotdir.ns <- paste0(here::here(), "/NSRKC/plots")
modeldir.ns <- paste0(here::here(), "/NSRKC/models")
outdir.ns <- paste0(here::here(), "/NSRKC/output")
datdir.ns <- paste0(here::here(), "/NSRKC/data")

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
# read in and process survey data ----
# *************************************************************************************************

nsrkc.survey <- read.csv(paste0(here::here(), "/NSRKC/data/NSRKC_trawl_survey_abundance.csv"))

nsrkc.dt <- nsrkc.survey %>%
  #select(Year, ADFG_Station, Latitude, Longitude, crab.km2) %>%
  arrange(Year, ADFG_Station) %>%
  mutate("year_f" = factor(Year)) %>%
  dplyr::select(-X) 

# check that all Lat/Lon combinations are valid. Norton Sound includes 
# 160 < longitude < 168 (source: https://www.arlis.org/docs/vol1/RIR/162606772.pdf)
# 61.49 < latitude < 66 (source: https://www.adfg.alaska.gov/static/applications/dcfnewsrelease/750733004.pdf)
filter(nsrkc.dt, Latitude > 66 | Latitude < 61.49)
filter(nsrkc.dt, Longitude > -160 | Longitude < -170)

nsrkc.dt.invalid <- nsrkc.dt %>%
  filter((Year == 2018 & ADFG_Station == 99 & Haul == 28) | (Year == 2018 & ADFG_Station == 127 & Haul == 26)) 

nsrkc.dt2 <- nsrkc.dt %>%
  #filter(Latitude <= 65.9) %>%
  #filter(Latitude >= 61.49) %>%
  #filter(Longitude <= -160) %>%
  #filter(Longitude >= -168) %>%
  mutate(Longitude = case_when(
    Year == 1999 & Agent == "ADFG" & Haul == 24 & ADFG_Station == 184 ~ -165.3512, # fixing typo
    Longitude > 0 ~ Longitude * -1, # fixing longitudes that should be negative 
    .default = Longitude
  )) %>%
  filter(Latitude <= 67) %>%
  filter(Latitude >= 61) %>%
  filter(Longitude <= -159) %>%
  filter(Longitude >= -169) %>%
  anti_join(nsrkc.dt.invalid)

filter(nsrkc.dt2, Latitude > 66 | Latitude < 61.49)
filter(nsrkc.dt2, Longitude > -160 | Longitude < -170)

# project the survey lat/lons to UTM

nsrkc_utm1 <- nsrkc.dt2 %>%
  dplyr::select(c(Longitude, Latitude)) %>%
  as.matrix() %>%
  terra::vect(crs="+proj=longlat +datum=WGS84")

nsrkc_utm2 <- project(nsrkc_utm1, "+proj=utm +zone=3 +datum=WGS84  +units=m")

nsrkc_utm3 <- geom(nsrkc_utm2)[, c("x", "y")]
head(nsrkc_utm3, 3)
nsrkc_utm4 <- nsrkc_utm3 %>%
  data.frame() %>%
  rename("X" = x, "Y" = y)

nsrkc_utm <- cbind(nsrkc.dt2, nsrkc_utm4) %>%
  # center and scale depth by its standard deviation
  mutate(depth_scaled = scale(Depth_m, center = TRUE, scale = TRUE)) %>%
  mutate(depth_scaled2 = depth_scaled ^ 2)

# plot all survey data
ggplot(nsrkc_utm) + 
  geom_point(aes(x = Longitude, y = Latitude, color = Agent)) +
  theme_bw() 

ggplot(nsrkc_utm) + 
  geom_point(aes(x = X, y = Y, color = Agent)) +
  theme_bw() 

# include only data from the ADF&G survey
nsrkc_utm_adfg <- nsrkc_utm %>%
  filter(Agent == "ADFG") %>%
  # convert UTM from meters to kilometers to reduce computing needs
  mutate(X = X/1000, Y = Y/1000)

# ADF&G survey data with depth information
nsrkc_utm_adfg_d <- nsrkc_utm_adfg %>% filter(is.na(Depth_m) == FALSE)

# include only data from the NOAA NS survey
nsrkc_utm_noaa <- nsrkc_utm %>%
  filter(Agent == "NOAA") %>%
  # convert UTM from meters to kilometers to reduce computing needs
  mutate(X = X/1000, Y = Y/1000)

# NOAA NS survey data with depth information
nsrkc_utm_noaa_d <- nsrkc_utm_noaa %>% filter(is.na(Depth_m) == FALSE)

# include only data from the NOAA NBS survey
nsrkc_utm_nbs <- nsrkc_utm %>%
  filter(Agent == "NBS") %>%
  # convert UTM from meters to kilometers to reduce computing needs
  mutate(X = X/1000, Y = Y/1000)

# NOAA NS survey data with depth information
nsrkc_utm_nbs_d <- nsrkc_utm_nbs %>% filter(is.na(Depth_m) == FALSE)

# mean depth per year
nsrkc_utm_adfg_d %>% group_by(Year) %>% summarise(mean.depth = mean(Depth_m))
nsrkc_utm_noaa_d %>% group_by(Year) %>% summarise(mean.depth = mean(Depth_m))
nsrkc_utm_nbs_d %>% group_by(Year) %>% summarise(mean.depth = mean(Depth_m))

# number of unique stations
length(unique(nsrkc_utm_adfg$ADFG_Station))

# check number of observations per year (this should be larger than the number of vertices in the mesh)
n_yr_adfg <- nsrkc_utm_adfg %>%
  group_by(Year) %>%
  count() 
median(n_yr_adfg$n)
max(n_yr_adfg$n)

n_yr_noaans <- nsrkc_utm_noaa %>%
  group_by(Year) %>%
  count() 
median(n_yr_noaans$n)
max(n_yr_noaans$n)

n_yr_nbs <- nsrkc_utm_nbs %>%
  group_by(Year) %>%
  count() 
median(n_yr_nbs$n)
max(n_yr_nbs$n)

# plot survey data - histograms
ggplot(nsrkc_utm_adfg) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "X", y = "Y")

hist.adfg <- ggplot(nsrkc_utm_adfg, aes(x = crab.km2)) +
  geom_histogram(color=FALSE, fill=cbpalette[1], binwidth = 70) +
  xlab(bquote("Estimated crab per" ~~ km^2)) +
  ylab("Number of estimates") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 450)) +
  geom_text(x=5500, y=420, label="ADF&G trawl survey", size = 4) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

hist.noaa <- ggplot(nsrkc_utm_noaa, aes(x = crab.km2)) +
  geom_histogram(color=FALSE, fill=cbpalette[2], binwidth = 30) +
  xlab(bquote("Estimated crab per" ~~ km^2)) +
  ylab("Number of estimates") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 250)) +
  geom_text(x=2000, y=220, label="NOAA Norton Sound \ntrawl survey", size = 4) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

hist.nbs <- ggplot(nsrkc_utm_nbs, aes(x = crab.km2)) +
  geom_histogram(color=FALSE, fill=cbpalette[3], binwidth = 10) +
  xlab(bquote("Estimated crab per" ~~ km^2)) +
  ylab("Number of estimates") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 125)) +
  geom_text(x=450, y=110, label="NOAA Northern Bering \nSea trawl survey", size = 4) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

hist.surveys <- ggarrange(hist.adfg, hist.noaa, hist.nbs, ncol = 2, nrow = 2)
ggsave(file.path(plotdir.ns, "hist_surveys.png"), plot = hist.surveys, height = 7, width = 7, units = "in")

# abundance for ADF&G survey
abund.crab.adfg <- nsrkc.dt %>%
  filter(Agent == "ADFG") %>%
  group_by(Year) %>%
  mutate(crab.ab = crab.km2 * area.km2) %>%
  mutate(obs_index = sum(crab.ab), sd.crab.ab = sd(crab.ab), mean.crab.ab = mean(crab.ab), n.station = length(ADFG_Station)) %>%
  mutate(se.crab.ab = sd.crab.ab / sqrt(n.station)) %>%
  mutate(cv.crab.ab = sd.crab.ab / mean.crab.ab) %>%
  mutate(obs.mean.lwr95 = mean.crab.ab - 1.96*se.crab.ab,
         obs.mean.upr95 = mean.crab.ab + 1.96*se.crab.ab,
         obs_l95 = obs.mean.lwr95 * n.station,
         obs_u95 = obs.mean.upr95 * n.station) %>%
  dplyr::select(Year, obs_index, mean.crab.ab, n.station, se.crab.ab, cv.crab.ab, obs.mean.lwr95, obs.mean.upr95, obs_l95, obs_u95) %>%
  ungroup() %>%
  distinct()

# abundance for surveys from assessment model data input file
abund.mod <- read.csv("C:/Users/castern/OneDrive - State of Alaska/Documents/BSAI_crab_assessments/NSRKC/nsrkc_2025_05/tables/index_obs_pred.csv") %>%
  dplyr::select(-pred_index)
  
abund.adfg.mod <- abund.mod %>%
  filter(fleet == "ADFG_Trawl")

abund.noaa.mod <- abund.mod %>%
  filter(fleet == "NMFS_Trawl")

abund.nbs.mod <- abund.mod %>%
  filter(fleet == "NBS_Trawl")

abund.adfg.mod %>% 
  ggplot() +
  geom_point(aes(x = year, y = obs_index), color = "grey20")+
  geom_errorbar(aes(x = year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20") +
  theme_bw() +
  scale_y_continuous(limits = c(0, max(abund.adfg.mod$obs_u95)), labels = label_comma())

abund.noaa.mod %>% 
  ggplot() +
  geom_point(aes(x = year, y = obs_index), color = "grey20")+
  geom_errorbar(aes(x = year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20") +
  theme_bw() +
  scale_y_continuous(limits = c(0, max(abund.adfg.mod$obs_u95)), labels = label_comma())

abund.nbs.mod %>% 
  ggplot() +
  geom_point(aes(x = year, y = obs_index), color = "grey20")+
  geom_errorbar(aes(x = year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20") +
  theme_bw() +
  scale_y_continuous(limits = c(0, max(abund.adfg.mod$obs_u95)), labels = label_comma())

abund.crab.adfg %>% 
  ggplot() +
  geom_point(aes(x = Year, y = obs_index), color = "grey20")+
  geom_errorbar(aes(x = Year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20") +
  theme_bw() +
  #scale_y_continuous(limits = c(0, max(abund.crab.adfg$obs_u95)), labels = label_comma())
  scale_y_continuous(labels = label_comma())

# plot showing distribution of observations
ggplot(nsrkc_utm_adfg, aes(x = ADFG_Station, y = crab.km2)) +
  geom_point()

# proportion of zeros in each data set
count.zero.adfg <- filter(nsrkc_utm_adfg, crab.km2 == 0) %>% count() 
count.obs.adfg <- length(nsrkc_utm_adfg$crab.km2)
prop.zero.adfg <- count.zero.adfg/count.obs.adfg
count.zero.noaa <- filter(nsrkc_utm_noaa, crab.km2 == 0) %>% count() 
count.obs.noaa <- length(nsrkc_utm_noaa$crab.km2)
prop.zero.noaa <- count.zero.noaa/count.obs.noaa
count.zero.nbs <- filter(nsrkc_utm_nbs, crab.km2 == 0) %>% count() 
count.obs.nbs <- length(nsrkc_utm_nbs$crab.km2)
prop.zero.nbs <- count.zero.nbs/count.obs.nbs

# *************************************************************************************************
# read in bathymetry data ----
# *************************************************************************************************

#depth1 <- rast(paste0(here::here(), "/NSRKC/data/so_ak_crm_v1_tif/so_ak_crm_v1.tif"))
#depth1 <- raster::raster(paste0(here::here(), "/NSRKC/data/so_ak_crm_v1_tif/so_ak_crm_v1.tif"))

# this raster has longitude in East longitude. To convert to West longitude,
# for E < 180, W = 360 - E
# for E > 180, W = E - 360
# Norton Sound bounds are -170 to -160 W = 190 to 200 E

# crop raster to Norton Sound area
#depth1.crop <- raster::crop(depth1, extent(190, 200, 61.49, 66.5))

#depth1.crop.df <- as.data.frame(depth1.crop, xy = TRUE) %>%
  #mutate(Longitude = -1*(360 - x), Latitude = y)

# join to pred_grid_ns_km
#depth.grid <- left_join(pred_grid_ns_km %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)), 
                        #depth1.crop.df %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)),
                        #by = c("Lon_join", "Lat_join")) %>%
  #distinct(X, Y, .keep_all = TRUE) %>%
  #rename(depth = so_ak_crm_v1)

# ggOceanMaps

library(ggOceanMaps)

extent.dat <- data.frame(lon = c(-170, -170, -163, -163), lat = c(62, 65, 62, 65))

basemap(data = extent.dat, bathymetry = TRUE) + 
  geom_polygon(data = transform_coord(extent.dat), aes(x = lon, y = lat), 
               color = "red", fill = NA)


# *************************************************************************************************
# read in grid used in VAST to use for predictions and convert to UTM ----
# *************************************************************************************************

grid5km.ns <- readRDS(paste0(here::here(), "/NSRKC/data/NSRKC_5km_grid_V2024/NSRKC_5km_grid_V2024.rds")) %>%
  # get UTM zone
  mutate(zone = (floor((Lon + 180)/6) %% 60) + 1) 

unique(grid5km.ns$zone)

predgrid_utm_z3.ns <- grid5km.ns %>% filter(zone == 3)
predgrid_utm_z4.ns <- grid5km.ns %>% filter(zone == 4)

u3.ns <- predgrid_utm_z3.ns
get_crs(u3.ns, c("Lon", "Lat"))
u3u.ns <- add_utm_columns(u3.ns, c("Lon", "Lat"))

u4.ns <- predgrid_utm_z4.ns
get_crs(u4.ns, c("Lon", "Lat"))
u4u.ns <- add_utm_columns(u4.ns, c("Lon", "Lat"))

predgrid_utm.ns <- rbind(u3u.ns, u4u.ns) %>%
  mutate(X = X / 1000, Y = Y / 1000) %>% # convert UTM coordinates from meter to kilometers to reduce computing needs
  rename("Longitude" = Lon, "Latitude" = Lat)

# plot prediction grid
predgrid5kmplot <- ggplot(predgrid_utm.ns, aes(x = Longitude, y = Latitude, color = Area_km2)) +
  geom_point() + 
  theme(legend.position="none") 

ggsave(file.path(plotdir.ns, "prediction_grid_5km_V2024.png"), plot = predgrid5kmplot, height = 4.2, width = 7, units = "in")

# plot all survey data with prediction grid
vast_grid_survey <- ggplot() +
  geom_point(data = predgrid_utm.ns, aes(x = Longitude, y = Latitude), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette)

ggsave(file.path(plotdir.ns, "vast_grid_survey.png"), plot = vast_grid_survey, height = 4.2, width = 7, units = "in")


# *************************************************************************************************
# create prediction grid including all survey sampling locations ----
# *************************************************************************************************

coastmap <- ne_countries(country = "United States of America", scale = 'medium', returnclass = 'sf') 

# create sf object
coastmap2 <- st_as_sf(coastmap)
coastmap1 <- st_as_sf(coastmap, coords = c("lon", "lat"), crs=32603)

# check object units
sf::st_length(coastmap2)
st_crs(coastmap2, parameters = TRUE)$units_gdal

# check projection
st_crs(coastmap1)$proj4string

# crop map to Norton Sound area
nsound <- st_crop(coastmap2, xmin=-170, xmax=-160, ymin=61.49, ymax=66.5)
plot(nsound$geometry)
st_crs(nsound, parameters = TRUE)$units_gdal
ggplot(nsound) +
  geom_sf()

# project the data to get a grid in units of meters
nsound2 <- st_transform(nsound, crs=32603)
st_crs(nsound2)$proj4string
plot(nsound2$geometry)
ggplot(nsound2) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# crop again to get area for prediction grid
nsound3 <- st_crop(nsound2, xmin=280000, xmax=730000, ymin=6970000, ymax=7300000)
ggplot(nsound3) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# now can set cellsize in meters; 5 km^2 = cellsize of 5000 x 5000
# ADFG survey and NOAA Norton Sound survey have grid size of 10 nm; 
# NBS survey has grid size of 20 nm. 1 nm = 1.852 km so grid size of 10 nm
# corresponds to grid size of 18.52 km and grid cell area of 342.99 km.
# Grid size of 20 nm = grid size of 37.04 km and grid cell area of 1371.962 km.
grid.sea <- nsound3 %>% 
  st_make_grid(cellsize = c(18520,18520), what = "centers") %>% # grid of points
  st_difference(nsound3) %>%
  st_crop(nsound2, xmin=280000, xmax=730000, ymin=6970000, ymax=7259000)

ggplot(grid.sea) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

grid.sea.nbs <- nsound3 %>% 
  st_make_grid(cellsize = c(37040,37040), what = "centers") %>% # grid of points
  st_difference(nsound3) %>%
  st_crop(nsound2, xmin=280000, xmax=730000, ymin=6970000, ymax=7259000)

ggplot(grid.sea.nbs) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# plot survey data on the grid to check overlap
ggplot() +
  geom_sf(data = grid.sea) +
  coord_sf(datum=st_crs(32603)) +
  geom_point(data = nsrkc_utm, aes(x = X, y = Y, color = Agent))

ggplot() +
  geom_sf(data = grid.sea.nbs) +
  coord_sf(datum=st_crs(32603)) +
  geom_point(data = nsrkc_utm, aes(x = X, y = Y, color = Agent))

# extract UTM coordinates from the prediction grid and plot with survey data
grid.sea.co <- st_coordinates(grid.sea)
grid.sea.nbs.co <- st_coordinates(grid.sea.nbs)

ns.pred.grid <- ggplot() +
  geom_point(data = grid.sea.co, aes(x = X/1000, y = Y/1000), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm %>% filter(Agent %in% c("NOAA","ADFG")), aes(x = X/1000, y = Y/1000, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  scale_x_continuous(breaks=seq(300, 700, by = 100), limits = c(280, 700)) +
  scale_y_continuous(breaks=seq(7000, 7260, by = 100), limits = c(6970, 7270))

nbs.pred.grid <- ggplot() +
  geom_point(data = grid.sea.nbs.co, aes(x = X/1000, y = Y/1000), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm %>% filter(Agent == "NBS"), aes(x = X/1000, y = Y/1000, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  scale_x_continuous(breaks=seq(300, 700, by = 100), limits = c(280, 700)) +
  scale_y_continuous(breaks=seq(7000, 7260, by = 100), limits = c(6970, 7270))

# export plot of prediction grid
ggsave(file.path(plotdir.ns, "ns_pred_grid.png"), plot = ns.pred.grid, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nbs_pred_grid.png"), plot = nbs.pred.grid, height = 4.2, width = 7, units = "in")

# convert the utm coordinates of the prediction grid back to lat/lon
v <- terra::vect(grid.sea.co, crs="+proj=utm +zone=3 +datum=WGS84  +units=m")
y <- project(v, "+proj=longlat +datum=WGS84")
lonlat <- geom(y)[, c("x", "y")]
head(lonlat, 3)
lonlat1 <- lonlat %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y)

v.nbs <- terra::vect(grid.sea.nbs.co, crs="+proj=utm +zone=3 +datum=WGS84  +units=m")
y.nbs <- project(v.nbs, "+proj=longlat +datum=WGS84")
lonlat.nbs <- geom(y.nbs)[, c("x", "y")]
head(lonlat.nbs, 3)
lonlat1.nbs <- lonlat.nbs %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y)

# plot survey data on the grid to check overlap
ggplot() +
  geom_point(data = lonlat1, aes(x = Longitude, y = Latitude), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette)

ggplot() +
  geom_point(data = lonlat1.nbs, aes(x = Longitude, y = Latitude), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette)


# complete prediction grid data frame
pred_grid_ns <- cbind(grid.sea.co, lonlat) %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y) %>%
  mutate(Area.km2 = 18.52^2)

pred_grid_nbs <- cbind(grid.sea.nbs.co, lonlat.nbs) %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y) %>%
  mutate(Area.km2 = (18.52*2)^2)

# prediction grid to use
# convert UTM coordinate units from m to km to reduce computing needs
pred_grid_ns_km <- pred_grid_ns %>%
  mutate(X = X/1000, Y = Y/1000)

pred_grid_nbs_km <- pred_grid_nbs %>%
  mutate(X = X/1000, Y = Y/1000)

# plot only ADFG survey data to check for errors
ggplot() +
  geom_point(data = lonlat1, aes(x = Longitude, y = Latitude)) +
  geom_point(data = nsrkc_utm %>% filter(Agent == "ADFG"), aes(x = Longitude, y = Latitude, color = factor(Year)))

#ggsave(file.path(plotdir.ns, "ns_adfg_mapped.png"), plot = ns.coastmap.adfg)

# map of area in lat/lon space
nsound3_latlon <- st_transform(nsound3, crs = "+proj=longlat +datum=WGS84")

ggplot() +
  geom_sf(data = nsound3_latlon) +
  geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent)) 

# map of survey data available by year ----

g <- purrr::map(levels(nsrkc_utm$year_f),
                function(x) {
                  ggplot() +
                    geom_sf(data = nsound3_latlon) +
                    geom_point(data = filter(nsrkc_utm, year_f == x), aes(x = Longitude, y = Latitude, color = factor(Agent, levels = c("ADFG", "NOAA", "NBS"))), show.legend = TRUE, size = 0.75) +
                    scale_color_manual(name = "Survey", 
                                       values = c("ADFG" = "#0072B2", "NOAA" = "#009E73", "NBS" = "#E69F00"),
                                       drop = FALSE) +
                    scale_x_continuous(breaks=seq(round(min(nsrkc_utm$Longitude),0), round(max(nsrkc_utm$Longitude),0), by = 6)) +
                    scale_y_continuous(breaks=seq(round(min(nsrkc_utm$Latitude),0), round(max(nsrkc_utm$Latitude),0), by = 1)) +
                    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
                })

survey_by_year1 <- ggarrange(plotlist = g, labels = levels(nsrkc_utm$year_f), common.legend = TRUE, legend="right")

survey_by_year <- annotate_figure(survey_by_year1, left = textGrob("Latitude", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Longitude", gp = gpar(cex = 1.3)))

ggsave(file.path(plotdir.ns, "survey_by_year.png"), plot = survey_by_year, height = 7.5, width = 10, units = "in", bg = "white")

# plot survey data on grid to check overlap
#ggplot() +
  #geom_sf(data = grid.sea) +
  #geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent))

# get coordinates of grid
#grid.sea.c <- st_coordinates(grid.sea) %>%
  #st_drop_geometry() %>%
 # data.frame() %>%
  #mutate(X = X/1000000, Y = Y/1000000, Area_km2 = 25, crs = 32603)

#ggplot() +
 # geom_point(data = grid.sea.c, aes(x = X, y = Y)) +
  #geom_point(data = nsrkc_utm, aes(x = X, y = Y, color = Agent))

# convert UTM to lat/lon
#grid.sea.c.ll <- grid.sea.c %>%
 # vect(geom=c("X", "Y"), crs="+proj=utm +zone=3") %>%
  #project("+proj=longlat") %>%
  #terra::as.data.frame(geom = "XY")

#grid.sea %>%
 # st_transform(crs = "+proj=longlat +datum=WGS84") %>%
  #st_coordinates() %>%
  #data.frame()

#grid.sea.c.ll2 <- grid.sea %>%
#  st_transform(crs = "+proj=longlat +datum=WGS84") %>%
#  st_coordinates() %>%
#  data.frame() %>%
#  rename("Longitude" = X, "Latitude" = Y)

#grid.sea.c.ll2 <- grid.sea.c %>% 
 # st_as_sf(coords = c("X", "Y"), crs = "+proj=utm +zone=3") %>%
  #st_transform(crs = "+proj=longlat +datum=WGS84") %>%
  #st_coordinates(st_cast(.$geometry,"MULTIPOINT")) %>%
  #data.frame() %>%
  #rename("Longitude" = X, "Latitude" = Y) %>%
  #cbind(grid.sea.c)
  
 
  
# not correct. Maybe try using lat/lons?
#start_crs <- "+proj=longlat +datum=WGS84 +no_defs"
#grid.sea.ll <- st_transform(grid.sea, crs = start_crs)

#grid.sea.ll.c <- st_coordinates(grid.sea.ll) %>%
  #data.frame() %>%
  #rename("Longitude" = X, "Latitude" = Y) %>%
  #mutate(Area_km2 = 25)

#ggplot() +
 # geom_point(data = grid.sea.c.ll, aes(x = Longitude, y = Latitude)) +
  #geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent))

# adjustments to prediction grid

#grid.sea.ll.use <- grid.sea.ll.c %>%
 # filter(Longitude > -168.1) %>%
  #filter(Latitude > 62.75)

#ggplot() +
 # geom_point(data = grid.sea.ll.use, aes(x = Longitude, y = Latitude)) +
  #geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent))

# plots of survey abundance ----

abund_adfg <- ggplot() + 
  geom_point(data = nsrkc_utm_adfg %>% filter(crab.km2 != 0), aes(x = Longitude, y = Latitude, size = crab.km2), color = cbpalette[1]) +
  geom_point(data = nsrkc_utm_adfg, aes(x = Longitude, y = Latitude), color = "black", size = 0.25) +
  geom_sf(data = nsound3_latlon) +
  scale_size_continuous(range = c(1, 4), breaks = seq(25, 8000, 2000)) +
  theme_bw() +
  scale_x_continuous(breaks=seq(round(min(nsrkc_utm$Longitude),0), round(max(nsrkc_utm$Longitude),0), by = 6)) +
  scale_y_continuous(breaks=seq(round(min(nsrkc_utm$Latitude),0), round(max(nsrkc_utm$Latitude),0), by = 1)) +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year) +
  ggtitle(ADFG~trawl~survey~estimated~crab~per~km^2)

abund_noaa <- ggplot() + 
  geom_point(data = nsrkc_utm_noaa %>% filter(crab.km2 != 0), aes(x = Longitude, y = Latitude, size = crab.km2), color = cbpalette[1]) +
  geom_point(data = nsrkc_utm_noaa, aes(x = Longitude, y = Latitude), color = "black", size = 0.25) +
  geom_sf(data = nsound3_latlon) +
  scale_size_continuous(range = c(1, 4), breaks = seq(20, 3000, 500)) +
  theme_bw() +
  scale_x_continuous(breaks=seq(round(min(nsrkc_utm$Longitude),0), round(max(nsrkc_utm$Longitude),0), by = 6)) +
  scale_y_continuous(breaks=seq(round(min(nsrkc_utm$Latitude),0), round(max(nsrkc_utm$Latitude),0), by = 1)) +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year) +
  ggtitle(NOAA~Norton~Sound~trawl~survey~estimated~crab~per~km^2)
       
abund_nbs <- ggplot() + 
  geom_point(data = nsrkc_utm_nbs %>% filter(crab.km2 != 0), aes(x = Longitude, y = Latitude, size = crab.km2), color = cbpalette[1]) +
  geom_point(data = nsrkc_utm_nbs, aes(x = Longitude, y = Latitude), color = "black", size = 0.25) +
  geom_sf(data = nsound3_latlon) +
  scale_size_continuous(range = c(1, 4), breaks = seq(25, 720, 100)) +
  scale_x_continuous(breaks=seq(round(min(nsrkc_utm$Longitude),0), round(max(nsrkc_utm$Longitude),0), by = 6)) +
  scale_y_continuous(breaks=seq(round(min(nsrkc_utm$Latitude),0), round(max(nsrkc_utm$Latitude),0), by = 1)) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year) +
  ggtitle(NOAA~Northern~Bering~Sea~trawl~survey~estimated~crab~per~km^2)
              
ggsave(file.path(plotdir.ns, "survey_abund_adfg.png"), plot = abund_adfg, height = 7, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "survey_abund_noaa.png"), plot = abund_noaa, height = 4.5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "survey_abund_nbs.png"), plot = abund_nbs, height = 4.5, width = 7, units = "in")

# *************************************************************************************************
# create prediction grid including only ADFG survey sampling locations since 2010 ----
# *************************************************************************************************

# crop map to Norton Sound area
nsound.s <- st_crop(coastmap2, xmin=-170, xmax=-160, ymin=61.49, ymax=66.5)
plot(nsound$geometry)
st_crs(nsound.s, parameters = TRUE)$units_gdal
ggplot(nsound.s) +
  geom_sf()

# project the data to get a grid in units of meters
nsound2.s <- st_transform(nsound.s, crs=32603)
st_crs(nsound2.s)$proj4string
plot(nsound2.s$geometry)
ggplot(nsound2.s) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# crop again to get area for prediction grid
nsound3.s <- st_crop(nsound2.s, xmin=280000, xmax=730000, ymin=6970000, ymax=7300000)
ggplot(nsound3.s) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# now can set cellsize in meters; 5 km^2 = cellsize of 5000 x 5000
# ADFG survey and NOAA Norton Sound survey have grid size of 10 nm; 
# NBS survey has grid size of 20 nm. 1 nm = 1.852 km so grid size of 10 nm
# corresponds to grid size of 18.52 km and grid cell area of 342.99 km.
# Grid size of 20 nm = grid size of 37.04 km and grid cell area of 1371.962 km.
grid.sea.s <- nsound3.s %>% 
  st_make_grid(cellsize = c(18520,18520), what = "centers") %>% # grid of points
  st_difference(nsound3.s) %>%
  st_crop(nsound2, xmin=380000, xmax=730000, ymin=7050000, ymax=7200000)

ggplot(grid.sea.s) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

grid.sea.nbs.s <- nsound3.s %>% 
  st_make_grid(cellsize = c(37040,37040), what = "centers") %>% # grid of points
  st_difference(nsound3.s) %>%
  st_crop(nsound2.s, xmin=380000, xmax=730000, ymin=7050000, ymax=7200000)

ggplot(grid.sea.nbs.s) +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# plot survey data on the grid to check overlap
nsrkc_utm_adfg10 <- nsrkc_utm %>%
  filter(Agent == "ADFG" & Year > 2010)

ggplot() +
  geom_sf(data = grid.sea.s) +
  coord_sf(datum=st_crs(32603)) +
  geom_point(data = nsrkc_utm_adfg10, aes(x = X, y = Y, color = Agent))

ggplot() +
  geom_sf(data = grid.sea.nbs.s) +
  coord_sf(datum=st_crs(32603)) +
  geom_point(data = nsrkc_utm_adfg10, aes(x = X, y = Y, color = Agent))

# extract UTM coordinates from the prediction grid and plot with survey data
grid.sea.co.s <- st_coordinates(grid.sea.s)
grid.sea.nbs.co.s <- st_coordinates(grid.sea.nbs.s)

ns.pred.grid.s <- ggplot() +
  geom_point(data = grid.sea.co.s, aes(x = X/1000, y = Y/1000), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm %>% filter(Agent %in% c("NOAA","ADFG")), aes(x = X/1000, y = Y/1000, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  scale_x_continuous(breaks=seq(300, 700, by = 100), limits = c(280, 700)) +
  scale_y_continuous(breaks=seq(7000, 7260, by = 100), limits = c(6970, 7270))

nbs.pred.grid.s <- ggplot() +
  geom_point(data = grid.sea.nbs.co.s, aes(x = X/1000, y = Y/1000), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm %>% filter(Agent == "NBS"), aes(x = X/1000, y = Y/1000, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  scale_x_continuous(breaks=seq(300, 700, by = 100), limits = c(280, 700)) +
  scale_y_continuous(breaks=seq(7000, 7260, by = 100), limits = c(6970, 7270))

# export plot of prediction grid
ggsave(file.path(plotdir.ns, "ns_pred_grid_s.png"), plot = ns.pred.grid.s, height = 4.2, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nbs_pred_grid_s.png"), plot = nbs.pred.grid.s, height = 4.2, width = 7, units = "in")

# convert the utm coordinates of the prediction grid back to lat/lon
v.s <- terra::vect(grid.sea.co.s, crs="+proj=utm +zone=3 +datum=WGS84  +units=m")
y.s <- project(v.s, "+proj=longlat +datum=WGS84")
lonlat.s <- geom(y.s)[, c("x", "y")]
head(lonlat.s, 3)
lonlat1.s <- lonlat.s %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y)

v.nbs.s <- terra::vect(grid.sea.nbs.co.s, crs="+proj=utm +zone=3 +datum=WGS84  +units=m")
y.nbs.s <- project(v.nbs.s, "+proj=longlat +datum=WGS84")
lonlat.nbs.s <- geom(y.nbs.s)[, c("x", "y")]
head(lonlat.nbs.s, 3)
lonlat1.nbs.s <- lonlat.nbs.s %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y)

# plot survey data on the grid to check overlap
ggplot() +
  geom_point(data = lonlat1.s, aes(x = Longitude, y = Latitude), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette)

ggplot() +
  geom_point(data = lonlat1.nbs.s, aes(x = Longitude, y = Latitude), shape = 4, size = 1) +
  geom_point(data = nsrkc_utm, aes(x = Longitude, y = Latitude, color = Agent)) + 
  labs(color = "Survey") +
  scale_color_discrete(type = cbpalette)

# complete prediction grid data frame
pred_grid_ns.s <- cbind(grid.sea.co.s, lonlat.s) %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y) %>%
  mutate(Area.km2 = 18.52^2)

pred_grid_nbs.s <- cbind(grid.sea.nbs.co.s, lonlat.nbs.s) %>%
  data.frame() %>%
  rename("Longitude" = x, "Latitude" = y) %>%
  mutate(Area.km2 = (18.52*2)^2)

# prediction grid to use
# convert UTM coordinate units from m to km to reduce computing needs
pred_grid_ns_km.s <- pred_grid_ns.s %>%
  mutate(X = X/1000, Y = Y/1000)

pred_grid_nbs_km.s <- pred_grid_nbs.s %>%
  mutate(X = X/1000, Y = Y/1000)


# calculate grid areas for comparison ----

# first transform grid back into polygon object, then calculate area

grid.sea.poly <- st_as_sf(grid.sea, coords = c("x", "y"), crs = 32603) %>% 
  summarise() %>%
  st_concave_hull(ratio = 0.1)

ggplot() +
  geom_sf(data= grid.sea.poly) +
  coord_sf(datum=st_crs(32603))

grid.sea.nbs.poly <- st_as_sf(grid.sea.nbs, coords = c("x", "y"), crs = 32603) %>% 
  summarise() %>%
  st_concave_hull(ratio = 0.1)

ggplot() +
  geom_sf(data= grid.sea.nbs.poly) +
  coord_sf(datum=st_crs(32603))

grid.sea.s.poly <- st_as_sf(grid.sea.s, coords = c("x", "y"), crs = 32603) %>% 
  summarise() %>%
  st_concave_hull(ratio = 0.1)

grid.sea.s.poly <- st_as_sf(grid.sea.s, crs = 32603) %>% 
  summarise() %>%
  st_concave_hull(ratio = 0.1)

ggplot() +
  geom_sf(data= grid.sea.s.poly) +
  coord_sf(datum=st_crs(32603))

grid.sea.nbs.s.poly <- st_as_sf(grid.sea.nbs.s, coords = c("x", "y"), crs = 32603) %>% 
  summarise() %>%
  st_concave_hull(ratio = 0.1)

ggplot() +
  geom_sf(data= grid.sea.nbs.s.poly) +
  coord_sf(datum=st_crs(32603))

st_area(grid.sea.poly)
st_area(grid.sea.nbs.poly)
st_area(grid.sea.s.poly)
st_area(grid.sea.nbs.s.poly)
# use the values for the non-NBS grids because when transforming grid to polygon,
# the lower spatial resolution results in lower calculated area even though the
# area used to generate the grid was equivalent.

# area currently used for the ADF&G and NBS survey abundance estimate is 5641 nm^2
# area currently used for the NOAA survey abundance estimate is 7600 nm^2
# convert to km^2: multiply by 3.4299
7600*3.4299
5641*3.4299

# *************************************************************************************************
# add depth information to all prediction grids ----
# *************************************************************************************************

# depth data source: https://www.ngdc.noaa.gov/mgg/coastal/s_alaska.html

depth1 <- raster::raster(paste0(here::here(), "/NSRKC/data/so_ak_crm_v1_tif/so_ak_crm_v1.tif"))

# this raster has longitude in East longitude. To convert to West longitude,
# for E < 180, W = 360 - E
# for E > 180, W = E - 360
# Norton Sound bounds are -170 to -160 W = 190 to 200 E

# crop raster to Norton Sound area
depth1.crop <- raster::crop(depth1, extent(190, 200, 61.49, 66.5))

# convert to data frame and convert longitudes to West
depth1.crop.df <- as.data.frame(depth1.crop, xy = TRUE) %>%
  mutate(Longitude = -1*(360 - x), Latitude = y) %>%
  rename(depth = so_ak_crm_v1) %>%
  mutate(depth_scaled = scale(depth, center = TRUE, scale = TRUE)) 

# join to prediction grids. Have to round lat/lons to 2 decimal places to get matches
pred_grid_ns_km_dep <- left_join(pred_grid_ns_km %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)), 
                        depth1.crop.df %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)),
                        by = c("Lon_join", "Lat_join")) %>%
  distinct(X, Y, .keep_all = TRUE)

pred_grid_nbs_km_dep <- left_join(pred_grid_nbs_km %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)), 
                                 depth1.crop.df %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)),
                                 by = c("Lon_join", "Lat_join")) %>%
  distinct(X, Y, .keep_all = TRUE)

pred_grid_ns_km.s_dep <- left_join(pred_grid_ns_km.s %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)), 
                                 depth1.crop.df %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)),
                                 by = c("Lon_join", "Lat_join")) %>%
  distinct(X, Y, .keep_all = TRUE)

pred_grid_nbs_km.s_dep <- left_join(pred_grid_nbs_km.s %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)), 
                                 depth1.crop.df %>% mutate(Lon_join = round(Longitude,2), Lat_join = round(Latitude,2)),
                                 by = c("Lon_join", "Lat_join")) %>%
  distinct(X, Y, .keep_all = TRUE)
  

# *************************************************************************************************
# spatial commercial harvest plots ----
# *************************************************************************************************

ns.ft <- read.csv(paste0(here::here(), "/NSRKC/data/NSRKC_fish_tickets.csv"))

# summarize harvest by year and stat area
ns.ft.sa.yr <- ns.ft %>%
  filter(Species.Code == 921 & King.Mgt.Area == "Q" & King.Mgt.District == "NORTSD") %>%
  group_by(Batch.Year, Stat.Area) %>%
  summarize(total.lb = sum(Whole.Weight..sum., na.rm = TRUE), total.crab = sum(Number.Of.Animals, na.rm = TRUE), dist = n_distinct(Vessel.CFEC.ID)) %>%
  mutate(mill.lb = total.lb/1000000) %>%
  rename(year = Batch.Year)

# summarize harvest by stat area
ns.ft.sa <- ns.ft %>%
  filter(Species.Code == 921 & King.Mgt.Area == "Q" & King.Mgt.District == "NORTSD") %>%
  group_by(Stat.Area) %>%
  summarize(total.lb = sum(Whole.Weight..sum., na.rm = TRUE), total.crab = sum(Number.Of.Animals, na.rm = TRUE), dist = n_distinct(Vessel.CFEC.ID)) %>%
  mutate(mill.lb = total.lb/1000000) 

# bring in stat areas shape file
# from https://soa-adfg.opendata.arcgis.com/datasets/groundfish-statistical-areas-2001
st_layers(paste0(datdir.ns, "/Groundfish_Statistical_Areas_2001/PVG_Statewide_2001_Present_GCS_WGS1984.shp"))

st.area <- st_read(paste0(datdir.ns, "/Groundfish_Statistical_Areas_2001/PVG_Statewide_2001_Present_GCS_WGS1984.shp"))

# define Northern Bering Sea stat areas from this map: https://www.adfg.alaska.gov/static/fishing/PDFs/commercial/chart04_nbs.pdf
nbs.st.areas <- c(
  # < 62 degrees latitude
  696130, 686130, 676130, 666131, 666132, 
  # > 62 and < 62.5 degrees latitude
  696200, 686200, 676200, 666200, 656201, 656202,
  # > 62.5 and < 63 degrees latitude
  696232, 696231, 686230, 676230, 666230, 656231, 646231, 656232, 646232,
  # > 63 and < 63.5 degrees latitude
  696301, 696302, 696303, 696304, 686301, 686302, 676300, 666300, 656300, 646301, 646302, 636301, 636302, 626301, 626302,
  # > 63.5 and < 64 degrees latitude
  696330, 686330, 676330, 666330, 656330, 646330, 636330, 626331, 626332, 616331, 616332, 616333,
  # > 64 and < 64.5 degrees latitude
  696400, 686400, 676400, 666401, 666402, 656401, 656402, 656403, 646401, 646402, 646403, 636401, 636402, 636403, 
  626401, 626402, 626403, 616401, 616402, 616403, 
  # > 64.5 and < 65 degrees latitude
  696430, 686431, 686432, 676430, 666431, 666432, 626433, 616432,
  # > 65 degrees latitude
  696500, 686500, 676501, 676502, 666501, 666502
)

st.area.nbs <- st.area %>%
  filter(STAT_AREA %in% nbs.st.areas) %>%
  dplyr::select(STAT_AREA, geometry) %>%
  rename("Stat.Area" = STAT_AREA) %>%
  # project to UTM
  st_transform(nsound.s, crs=32603)

st.area.nbs %>% 
  ggplot() +
  geom_sf() +
  coord_sf(datum=st_crs(32603))

# plot with all crab caught in all years
st.area.nbs.ft <- st.area.nbs %>%
  left_join(ns.ft.sa, by = "Stat.Area") %>%
  # harvest present/absent field needed due to confidentiality requirements
  mutate(harvest.present = case_when(
    is.na(total.crab) == TRUE ~ "No",
    is.na(total.crab) == FALSE ~ "Yes"
  ))

ggplot(data = st.area.nbs.ft) +
  geom_sf(data = st.area.nbs.ft,  aes(fill = total.crab)) +
  coord_sf(datum=st_crs(32603)) 

# plot with crab caught in each year
st.area.nbs.ft.yr <-
  # replicate all stat areas for all years 1985 - present
  map(seq(1985, 2024), ~st.area.nbs) %>%
  bind_rows(.id = "id") %>%
  mutate(id = as.numeric(id)) %>%
  mutate(year = id + 1984) %>%
  # join to harvest information
  full_join(ns.ft.sa.yr, by = c("year", "Stat.Area")) %>%
  # harvest present/absent field needed due to confidentiality requirements
  mutate(harvest.present = case_when(
    is.na(total.crab) == TRUE ~ "No",
    is.na(total.crab) == FALSE ~ "Yes"
    ))

ggplot(data = st.area.nbs.ft.yr) +
  geom_sf(data = st.area.nbs.ft.yr,  aes(fill = total.crab)) +
  facet_wrap(~year) +
  coord_sf(datum=st_crs(32603))

# plot of full ADFG and NOAA prediction grid with harvest info - all years combined
harvest.full.all <- ggplot() +
  geom_sf(data = st.area.nbs.ft,  aes(fill = harvest.present)) +
  geom_point(data = grid.sea.co, aes(x = X, y = Y), color = "white") +
  #scale_fill_viridis_c(na.value = "gray") +
  scale_fill_manual(values = c("gray", cbpalette[2]), na.value = "gray") +
  scale_x_continuous(labels=function(x)x/1000, breaks=seq(300000, 700000, by = 100000), limits = c(280000, 700000)) +
  scale_y_continuous(labels=function(x)x/1000, breaks=seq(7000000, 7260000, by = 100000), limits = c(6970000, 7270000)) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  coord_sf(datum=st_crs(32603)) + 
  guides(fill=guide_legend(title="Crab \nharvested", reverse = T))

# plot of full ADFG and NOAA prediction grid with harvest info - year by year
harvest.full.yrs <- ggplot() +
  geom_sf(data = st.area.nbs.ft.yr %>% filter(year > 2016),  aes(fill = harvest.present)) +
  geom_point(data = grid.sea.co, aes(x = X, y = Y), color = "white", size = 0.5) +
  #scale_fill_viridis_c(na.value = "gray") +
  scale_fill_manual(values = c("gray", cbpalette[2]), na.value = "gray") +
  scale_x_continuous(labels=function(x)x/1000, breaks=seq(300000, 700000, by = 100000), limits = c(280000, 700000)) +
  scale_y_continuous(labels=function(x)x/1000, breaks=seq(7000000, 7260000, by = 100000), limits = c(6970000, 7270000)) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  facet_wrap(~year) +
  coord_sf(datum=st_crs(32603)) + 
  guides(fill=guide_legend(title="Crab \nharvested", reverse = T))

# plot of reduced ADFG and NOAA prediction grid with harvest info - all years combined
harvest.reduced.all <- ggplot(data = st.area.nbs.ft) +
  #geom_sf(data = st.area.nbs.ft,  aes(fill = total.crab)) +
  geom_sf(data = st.area.nbs.ft,  aes(fill = harvest.present)) +
  coord_sf(datum=st_crs(32603)) +
  geom_point(data = grid.sea.co.s, aes(x = X, y = Y), color = "white") +
  #scale_fill_viridis_c(na.value = "gray") +
  scale_fill_manual(values = c("gray", cbpalette[2]), na.value = "gray") +
  scale_x_continuous(labels=function(x)x/1000, breaks=seq(300000, 700000, by = 100000), limits = c(280000, 700000)) +
  scale_y_continuous(labels=function(x)x/1000, breaks=seq(7000000, 7260000, by = 100000), limits = c(6970000, 7270000)) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") + 
  guides(fill=guide_legend(title="Crab \nharvested", reverse = T))

# plot of reduced ADFG and NOAA prediction grid with harvest info - year by year
harvest.reduced.yrs <- ggplot() +
  #geom_sf(data = st.area.nbs.ft.yr %>% filter(year > 2016),  aes(fill = total.crab)) +
  geom_sf(data = st.area.nbs.ft.yr %>% filter(year > 2016),  aes(fill = harvest.present)) +
  geom_point(data = grid.sea.co.s, aes(x = X, y = Y), color = "white", size = 0.5) +
  #scale_fill_viridis_c(na.value = "gray") +
  scale_fill_manual(values = c("gray", cbpalette[2]), na.value = "gray") +
  scale_x_continuous(labels=function(x)x/1000, breaks=seq(300000, 700000, by = 100000), limits = c(280000, 700000)) +
  scale_y_continuous(labels=function(x)x/1000, breaks=seq(7000000, 7260000, by = 100000), limits = c(6970000, 7270000)) +
  xlab("Eastings (km)") +
  ylab("Northings (km)") + 
  guides(fill=guide_legend(title="Crab \nharvested", reverse = T)) +
  facet_wrap(~year) +
  coord_sf(datum=st_crs(32603))

ggsave(file.path(plotdir.ns, "harvest_full_yrs.png"), plot = harvest.full.yrs, height = 7, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "harvest_full_all.png"), plot = harvest.full.all, height = 7, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "harvest_reduced_yrs.png"), plot = harvest.reduced.yrs, height = 7, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "harvest_reduced_all.png"), plot = harvest.reduced.all, height = 7, width = 7, units = "in")

# *************************************************************************************************
# make SPDE mesh ----
# *************************************************************************************************

# ADFG mesh
# median stations = 54.5, max = 100

adfg_mesh <- make_mesh(nsrkc_utm_adfg, 
                            xy_cols = c("X","Y"), 
                            fmesher_func = fmesher::fm_mesh_2d_inla,
                            cutoff = 40)
                            #offset = c(10, 70))
plot(adfg_mesh)
adfg_mesh$mesh$n # 51

adfg_barriermesh <- sdmTMBextra::add_barrier_mesh(
  spde_obj = adfg_mesh,
  barrier_sf = nsound3,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
plot(adfg_barriermesh)

adfg_mesh_depth <- make_mesh(nsrkc_utm_adfg_d, 
                             xy_cols = c("X","Y"), 
                             fmesher_func = fmesher::fm_mesh_2d_inla,
                             cutoff = 40)
plot(adfg_mesh_depth)
adfg_mesh_depth$mesh$n

adfg_barriermesh_depth <- sdmTMBextra::add_barrier_mesh(
  spde_obj = adfg_mesh_depth,
  barrier_sf = nsound3,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
plot(adfg_barriermesh_depth)

# NOAA NS mesh
# median stations = 78.5, max = 104

noaa_mesh <- make_mesh(nsrkc_utm_noaa, 
                       xy_cols = c("X","Y"), 
                       fmesher_func = fmesher::fm_mesh_2d_inla,
                       cutoff = 25)
                       #offset = c(10, 70))
plot(noaa_mesh)
noaa_mesh$mesh$n # 77

noaa_barriermesh <- sdmTMBextra::add_barrier_mesh(
  spde_obj = noaa_mesh,
  barrier_sf = nsound3,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
plot(noaa_barriermesh)

noaa_mesh_depth <- make_mesh(nsrkc_utm_noaa_d, 
                       xy_cols = c("X","Y"), 
                       fmesher_func = fmesher::fm_mesh_2d_inla,
                       cutoff = 25)
#offset = c(10, 70))
plot(noaa_mesh_depth)
noaa_mesh_depth$mesh$n # 77

noaa_barriermesh_depth <- sdmTMBextra::add_barrier_mesh(
  spde_obj = noaa_mesh_depth,
  barrier_sf = nsound3,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
plot(noaa_barriermesh_depth)

# NBS mesh
# median = 35, range = 34 to 35

nbs_mesh <- make_mesh(nsrkc_utm_nbs, 
                        xy_cols = c("X","Y"), 
                        fmesher_func = fmesher::fm_mesh_2d_inla,
                        cutoff = 45)
                        #offset = c(10, 70))
plot(nbs_mesh)
nbs_mesh$mesh$n # 33

nbs_barriermesh <- sdmTMBextra::add_barrier_mesh(
  spde_obj = nbs_mesh,
  barrier_sf = nsound3,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
plot(nbs_barriermesh)

nbs_mesh_depth <- make_mesh(nsrkc_utm_nbs_d, 
                             xy_cols = c("X","Y"), 
                             fmesher_func = fmesher::fm_mesh_2d_inla,
                             cutoff = 25)
#offset = c(10, 70))
plot(nbs_mesh_depth)
nbs_mesh_depth$mesh$n # 77

nbs_barriermesh_depth <- sdmTMBextra::add_barrier_mesh(
  spde_obj = nbs_mesh_depth,
  barrier_sf = nsound3,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
plot(nbs_barriermesh_depth)

# ADFG mesh plots

mesh.plot.adfg <- ggplot(nsrkc_utm_adfg) + 
  inlabru::gg(adfg_barriermesh$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "Eastings (km)", y = "Northings (km)", title = "SPDE mesh for ADF&G trawl survey with 51 vertices")

ggsave(file.path(plotdir.ns, "mesh_adfg.png"), plot = mesh.plot.adfg, height = 4.2, width = 7, units = "in")

# plot barrier mesh vertices with land
ggplot() +
  geom_sf(data = nsound3) +
  geom_sf(data = adfg_barriermesh$mesh_sf) #+
  #coord_sf(datum=st_crs(32603))
   
ggplot() + 
  inlabru::gg(adfg_barriermesh$mesh)

ggplot() +
  #geom_point(data = adfg_mesh_100kn$loc_xy, aes(x = X, y = Y)) +
  geom_sf(data = nsound3) +
  coord_sf(datum=st_crs(32603))

#ggplot() +
 # inlabru::gg(adfg_mesh_100kn$mesh$graph$vt) + 
  #geom_point() +
  #geom_sf(data = nsound3) +
  #coord_sf(datum=st_crs(32603))

ggplot(nsrkc_utm_adfg) + 
  #inlabru::gg(adfg_mesh_100kn$mesh) + 
  geom_point(aes(x = X*1000, y = Y*1000), color = "blue") +
  theme_bw() +
  labs(x = "X", y = "Y") +
  geom_sf(data = nsound3) +
  geom_sf(data = adfg_barriermesh$mesh_sf)

ggplot() + 
  theme_bw() +
  geom_sf(data = nsound3) +
  geom_sf(data = adfg_barriermesh$mesh_sf)

# NOAA mesh plot

mesh.plot.noaa <- ggplot(nsrkc_utm_noaa) + 
  inlabru::gg(noaa_barriermesh$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "Eastings (km)", y = "Northings (km)", title = "SPDE mesh for NOAA Norton Sound trawl survey with 77 vertices")

ggsave(file.path(plotdir.ns, "mesh_noaa.png"), plot = mesh.plot.noaa, height = 4.2, width = 7, units = "in")

# NBS mesh plot

mesh.plot.nbs <- ggplot(nsrkc_utm_nbs) + 
  inlabru::gg(nbs_barriermesh$mesh) + 
  geom_point(aes(x = X, y = Y)) +
  theme_bw() +
  labs(x = "Eastings (km)", y = "Northings (km)", title = "SPDE mesh for NOAA Northern Bering Sea trawl survey with 33 vertices")

ggsave(file.path(plotdir.ns, "mesh_nbs.png"), plot = mesh.plot.nbs, height = 4.2, width = 7, units = "in")


# *************************************************************************************************
# fit and check models ----
# *************************************************************************************************

# ADFG survey models ----

# IID Tweedie model
m.adfg.iid.tw <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  #extra_time = c(1997,1998,2000,2001,2003,2004,2005,2007,2009,2010,2012,2013,2015,2016,2022),
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"))

# IID delta gamma model
m.adfg.iid.dg <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"))

m.adfg.iid.dgl <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(link1 = "logit", link2 = "log"))

m.adfg.iid.dgp <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "poisson-link"))

m.adfg.iid.dgm <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma_mix(link1 = "logit", link2 = "log"))

# IID delta lognormal model
m.adfg.iid.dl <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

m.adfg.dl.d <- sdmTMB(
  data = nsrkc_utm_adfg_d, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

# IID negative bionomial model

m.adfg.iid.nb <- sdmTMB(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg$Year), max(nsrkc_utm_adfg$Year)),
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = nbinom1(link = "log"))

# ADFG models with depth covariate ----

m.adfg.tw.d <- sdmTMB(
  data = nsrkc_utm_adfg_d,
  formula = crab.km2 ~ depth_scaled + year_f, 
  #formula = crab.km2 ~ depth_scaled2, 
  #formula = crab.km2 ~ Depth_m + year_f,
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg_d$Year), max(nsrkc_utm_adfg_d$Year)),
  mesh = adfg_barriermesh_depth,
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"))

m.adfg.dg.d <- sdmTMB(
  data = nsrkc_utm_adfg_d, 
  #formula = crab.km2 ~ Depth_m + year_f, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg_d$Year), max(nsrkc_utm_adfg_d$Year)),
  mesh = adfg_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"))

m.adfg.dl.d <- sdmTMB(
  data = nsrkc_utm_adfg_d, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg_d$Year), max(nsrkc_utm_adfg_d$Year)),
  mesh = adfg_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

# NOAA Norton Sound models ----

# IID Tweedie model
m.noaa.iid.tw <- sdmTMB(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"))

# IID delta gamma model
m.noaa.iid.dg <- sdmTMB(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"))

# IID delta lognormal model
m.noaa.iid.dl <- sdmTMB(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

# IID negative bionomial model
m.noaa.iid.nb <- sdmTMB(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = nbinom1(link = "log"))

# NOAA Norton Sound models with depth covariate ----

m.noaa.tw.d <- sdmTMB(
  data = nsrkc_utm_noaa_d,
  formula = crab.km2 ~ depth_scaled + year_f, 
  #formula = crab.km2 ~ depth_scaled2, 
  #formula = crab.km2 ~ Depth_m + year_f,
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa_d$Year), max(nsrkc_utm_noaa_d$Year)),
  mesh = noaa_barriermesh_depth,
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"))

m.noaa.dg.d <- sdmTMB(
  data = nsrkc_utm_noaa_d, 
  #formula = crab.km2 ~ Depth_m + year_f, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa_d$Year), max(nsrkc_utm_noaa_d$Year)),
  mesh = noaa_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"))

m.noaa.dl.d <- sdmTMB(
  data = nsrkc_utm_noaa_d, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa_d$Year), max(nsrkc_utm_noaa_d$Year)),
  mesh = noaa_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

# NOAA NBS models ----

# IID Tweedie model
m.nbs.iid.tw <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs$Year), max(nsrkc_utm_nbs$Year)),
  mesh = nbs_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"))

# IID delta gamma model
m.nbs.iid.dg <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs$Year), max(nsrkc_utm_nbs$Year)),
  mesh = nbs_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"))

# IID delta lognormal model
m.nbs.iid.dl <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs$Year), max(nsrkc_utm_nbs$Year)),
  mesh = nbs_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

# IID negative bionomial model
m.nbs.iid.nb <- sdmTMB(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs$Year), max(nsrkc_utm_nbs$Year)),
  mesh = nbs_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = nbinom1(link = "log"))

# NOAA Northern Bering Sea models with depth covariate ----

m.nbs.tw.d <- sdmTMB(
  data = nsrkc_utm_nbs_d,
  formula = crab.km2 ~ depth_scaled + year_f, 
  #formula = crab.km2 ~ depth_scaled2, 
  #formula = crab.km2 ~ Depth_m + year_f,
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs_d$Year), max(nsrkc_utm_nbs_d$Year)),
  mesh = nbs_barriermesh_depth,
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"))

m.nbs.dg.d <- sdmTMB(
  data = nsrkc_utm_nbs_d, 
  #formula = crab.km2 ~ Depth_m + year_f, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs_d$Year), max(nsrkc_utm_nbs_d$Year)),
  mesh = nbs_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"))

m.nbs.dl.d <- sdmTMB(
  data = nsrkc_utm_nbs_d, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs_d$Year), max(nsrkc_utm_nbs_d$Year)),
  mesh = nbs_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"))

# save fitted models ----

# save ADFG models
saveRDS(m.adfg.iid.tw, file = paste0(modeldir.ns, "/m_adfg_iid_tw.RDS"))
saveRDS(m.adfg.tw.d, file = paste0(modeldir.ns, "/m_adfg_tw_d.RDS"))
saveRDS(m.adfg.dg.d, file = paste0(modeldir.ns, "/m_adfg_dg_d.RDS"))

# save NOAA NS models
saveRDS(m.noaa.iid.tw, file = paste0(modeldir.ns, "/m_noaa_iid_tw.RDS"))
saveRDS(m.noaa.iid.dg, file = paste0(modeldir.ns, "/m_noaa_iid_dg.RDS"))
saveRDS(m.noaa.iid.dl, file = paste0(modeldir.ns, "/m_noaa_iid_dl.RDS"))
saveRDS(m.noaa.tw.d, file = paste0(modeldir.ns, "/m_noaa_tw_d.RDS"))
saveRDS(m.noaa.dg.d, file = paste0(modeldir.ns, "/m_noaa_dg_d.RDS"))
saveRDS(m.noaa.dl.d, file = paste0(modeldir.ns, "/m_noaa_dl_d.RDS"))

# save NOAA NBS models
saveRDS(m.nbs.iid.tw, file = paste0(modeldir.ns, "/m_nbs_iid_tw.RDS"))
saveRDS(m.nbs.tw.d, file = paste0(modeldir.ns, "/m_nbs_tw_d.RDS"))

# read in saved models
#m.adfg.iid.tw <- readRDS(paste0(modeldir.ns, "/m_adfg_iid_tw.RDS"))
#m.adfg.tw.d <- readRDS(paste0(modeldir.ns, "/m_adfg_tw_d.RDS"))
#m.adfg.dg.d <- readRDS(paste0(modeldir.ns, "/m_adfg_dg_d.RDS"))
#m.noaa.iid.tw <- readRDS(paste0(modeldir.ns, "/m_noaa_iid_tw.RDS"))
#m.noaa.iid.dg <- readRDS(paste0(modeldir.ns, "/m_noaa_iid_dg.RDS"))
#m.noaa.iid.dl <- readRDS(paste0(modeldir.ns, "/m_noaa_iid_dl.RDS"))
#m.noaa.tw.d <- readRDS(paste0(modeldir.ns, "/m_noaa_tw.d.RDS"))
#m.noaa.dg.d <- readRDS(paste0(modeldir.ns, "/m_noaa_dg.d.RDS"))
#m.noaa.dl.d <- readRDS(paste0(modeldir.ns, "/m_noaa_dl.d.RDS"))
#m.nbs.iid.tw <- readRDS(paste0(modeldir.ns, "/m_nbs_iid_tw.RDS"))
#m.nbs.tw.d <- readRDS(paste0(modeldir.ns, "/m_nbs_tw.d.RDS"))

# run sanity checks ----

# ADFG models
sanity(m.adfg.iid.tw) # pass
sanity(m.adfg.iid.dg) # fail
sanity(m.adfg.iid.dgp) # fail
sanity(m.adfg.iid.dgm) # fail
sanity(m.adfg.iid.dgl) # fail
sanity(m.adfg.iid.dl) # fail
sanity(m.adfg.iid.nb) # fail

# ADFG models - year + depth
sanity(m.adfg.tw.d) # pass
sanity(m.adfg.dg.d) # pass
sanity(m.adfg.dl.d) # fail

# NOAA Norton Sound models
sanity(m.noaa.iid.tw) # pass
sanity(m.noaa.iid.dg) # pass
sanity(m.noaa.iid.dl) # pass
sanity(m.noaa.iid.nb) # pass

# NOAA Norton Sound models - year + depth
sanity(m.noaa.tw.d) # pass
sanity(m.noaa.dg.d) # pass
sanity(m.noaa.dl.d) # pass

# NOAA NBS models
sanity(m.nbs.iid.tw) # pass
sanity(m.nbs.iid.dg) # fail
sanity(m.nbs.iid.dl) # fail
sanity(m.nbs.iid.nb) # pass

# NOAA NBS models - year + depth
sanity(m.nbs.tw.d) # pass
sanity(m.nbs.dg.d) # fail
sanity(m.nbs.dl.d) # fail


# ************************************************************************************************
# compare models using cross-validation ----
# ************************************************************************************************

# ADFG models - cross-validation ----

m.adfg.iid.tw.cv <- sdmTMB_cv(
  data = nsrkc_utm_adfg, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = adfg_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"),
  k_folds = 10)

m.adfg.tw.d.cv <- sdmTMB_cv(
  data = nsrkc_utm_adfg_d,
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg_d$Year), max(nsrkc_utm_adfg_d$Year)),
  mesh = adfg_barriermesh_depth,
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"),
  k_folds = 10)

m.adfg.dg.d.cv <- sdmTMB_cv(
  data = nsrkc_utm_adfg_d, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_adfg_d$Year), max(nsrkc_utm_adfg_d$Year)),
  mesh = adfg_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"),
  k_folds = 10)

# compare log-likelihood values
m.adfg.iid.tw.cv$sum_loglik
m.adfg.tw.d.cv$sum_loglik
m.adfg.dg.d.cv$sum_loglik

# compare RMSE (root mean square error) and MAE (mean absolute error)
rmse.adfg.tw.d <- sqrt(mean((m.adfg.tw.d.cv$data$crab.km2 - m.adfg.tw.d.cv$data$cv_predicted)^2))
mae.adfg.tw.d <- mean(abs(m.adfg.dg.d.cv$data$crab.km2 - m.adfg.dg.d.cv$data$cv_predicted))

rmse.adfg.dg.d <- sqrt(mean((m.adfg.dg.d.cv$data$crab.km2 - m.adfg.dg.d.cv$data$cv_predicted)^2))
mae.adfg.dg.d <- mean(abs(m.adfg.dg.d.cv$data$crab.km2 - m.adfg.dg.d.cv$data$cv_predicted))


# save cross validation models
saveRDS(m.adfg.iid.tw.cv, file = paste0(modeldir.ns, "/m_adfg_iid_tw_cv.RDS"))
saveRDS(m.adfg.tw.d.cv, file = paste0(modeldir.ns, "/m_adfg_tw_d_cv.RDS"))
saveRDS(m.adfg.dg.d.cv, file = paste0(modeldir.ns, "/m_adfg_dg_d_cv.RDS"))


# NOAA models - cross-validation ----

m.noaa.iid.tw.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"),
  k_folds = 10)

m.noaa.iid.dg.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"),
  k_folds = 10)

m.noaa.iid.dl.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"),
  k_folds = 10)

m.noaa.iid.nb.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa$Year), max(nsrkc_utm_noaa$Year)),
  mesh = noaa_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = nbinom1(link = "log"),
  k_folds = 10)

m.noaa.tw.d.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa_d,
  formula = crab.km2 ~ depth_scaled + year_f, 
  #formula = crab.km2 ~ depth_scaled2, 
  #formula = crab.km2 ~ Depth_m + year_f,
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa_d$Year), max(nsrkc_utm_noaa_d$Year)),
  mesh = noaa_barriermesh_depth,
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"),
  k_folds = 10)

m.noaa.dg.d.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa_d, 
  #formula = crab.km2 ~ Depth_m + year_f, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa_d$Year), max(nsrkc_utm_noaa_d$Year)),
  mesh = noaa_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_gamma(type = "standard"),
  k_folds = 10)

m.noaa.dl.d.cv <- sdmTMB_cv(
  data = nsrkc_utm_noaa_d, 
  formula = crab.km2 ~ depth_scaled + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_noaa_d$Year), max(nsrkc_utm_noaa_d$Year)),
  mesh = noaa_barriermesh_depth, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = delta_lognormal(type = "standard"),
  k_folds = 10)

# save cross validation models
saveRDS(m.noaa.iid.tw.cv, file = paste0(modeldir.ns, "/m_noaa_iid_tw_cv.RDS"))
saveRDS(m.noaa.iid.dg.cv, file = paste0(modeldir.ns, "/m_noaa_iid_dg_cv.RDS"))
saveRDS(m.noaa.iid.dl.cv, file = paste0(modeldir.ns, "/m_noaa_iid_dl_cv.RDS"))
saveRDS(m.noaa.tw.d.cv, file = paste0(modeldir.ns, "/m_noaa_tw_d_cv.RDS"))
saveRDS(m.noaa.dg.d.cv, file = paste0(modeldir.ns, "/m_noaa_dg_d_cv.RDS"))
saveRDS(m.noaa.dl.d.cv, file = paste0(modeldir.ns, "/m_noaa_dl_d_cv.RDS"))

# compare model predictive skill
m.noaa.iid.tw.cv$sum_loglik 
m.noaa.iid.dg.cv$sum_loglik 
m.noaa.iid.dl.cv$sum_loglik 
#m.noaa.iid.nb.cv$sum_loglik
m.noaa.tw.d.cv$sum_loglik
m.noaa.dg.d.cv$sum_loglik
m.noaa.dl.d.cv$sum_loglik

# compare RMSE (root mean square error) and MAE (mean absolute error)
rmse.noaa.tw.d <- sqrt(mean((m.noaa.tw.d.cv$data$crab.km2 - m.noaa.tw.d.cv$data$cv_predicted)^2))
mae.noaa.tw.d <- mean(abs(m.noaa.dg.d.cv$data$crab.km2 - m.noaa.dg.d.cv$data$cv_predicted))

rmse.noaa.dg.d <- sqrt(mean((m.noaa.dg.d.cv$data$crab.km2 - m.noaa.dg.d.cv$data$cv_predicted)^2))
mae.noaa.dg.d <- mean(abs(m.noaa.dg.d.cv$data$crab.km2 - m.noaa.dg.d.cv$data$cv_predicted))

rmse.noaa.dl.d <- sqrt(mean((m.noaa.dl.d.cv$data$crab.km2 - m.noaa.dl.d.cv$data$cv_predicted)^2))
mae.noaa.dl.d <- mean(abs(m.noaa.dl.d.cv$data$crab.km2 - m.noaa.dl.d.cv$data$cv_predicted))

# NBS models - cross-validation ----

m.nbs.iid.tw.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs, 
  formula = crab.km2 ~ 0 + year_f, 
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs$Year), max(nsrkc_utm_nbs$Year)),
  mesh = nbs_barriermesh, 
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"),
  k_folds = 10)

m.nbs.tw.d.cv <- sdmTMB_cv(
  data = nsrkc_utm_nbs_d,
  formula = crab.km2 ~ depth_scaled + year_f, 
  #formula = crab.km2 ~ depth_scaled2, 
  #formula = crab.km2 ~ Depth_m + year_f,
  spatial = "on",
  time = "Year", 
  extra_time = seq(min(nsrkc_utm_nbs_d$Year), max(nsrkc_utm_nbs_d$Year)),
  mesh = nbs_barriermesh_depth,
  spatiotemporal = "iid",
  silent = FALSE,
  anisotropy = FALSE,
  family = tweedie(link = "log"),
  k_folds = 10)

# compare model predictive skill
m.nbs.iid.tw.cv$sum_loglik 
m.nbs.tw.d.cv$sum_loglik

# compare RMSE (root mean square error) and MAE (mean absolute error)
sqrt(mean((m.nbs.tw.d.cv$data$density - m.nbs.tw.d.cv$data$cv_predicted)^2)) 
mean(abs(m.nbs.tw.d.cv$data$density - m.nbs.tw.d.cv$data$cv_predicted))

# save cross validation models
saveRDS(m.nbs.iid.tw.cv, file = paste0(modeldir.ns, "/m_nbs_iid_tw_cv.RDS"))
saveRDS(m.nbs.tw.d.cv, file = paste0(modeldir.ns, "/m_nbs_tw_d_cv.RDS"))

# ************************************************************************************************
# model residuals ----
# *************************************************************************************************

# ADFG model residuals ----
# TW + year
m.adfg.iid.tw.resid <- simulate(update(m.adfg.iid.tw, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.adfg.iid.tw, return_DHARMa = TRUE)
plot(m.adfg.iid.tw.resid)
# TW + year + depth
m.adfg.tw.d.resid <- simulate(update(m.adfg.tw.d, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.adfg.tw.d, return_DHARMa = TRUE)
# DG + year + depth
m.adfg.dg.d.resid <- simulate(update(m.adfg.dg.d, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.adfg.dg.d, return_DHARMa = TRUE)

saveRDS(m.adfg.iid.tw.resid, file = paste0(modeldir.ns, "/m_adfg_resid_iid_tw.RDS"))
saveRDS(m.adfg.tw.d.resid, file = paste0(modeldir.ns, "/m_adfg_resid_iid_tw.RDS"))
saveRDS(m.adfg.dg.d.resid, file = paste0(modeldir.ns, "/m_adfg_resid_iid_tw.RDS"))
#m.adfg.iid.tw.resid <- readRDS(paste0(modeldir.ns, "/m_adfg_resid_iid_tw.RDS"))

# residuals tests
# TW model year only
adfg.iid.tw.resid.outl <- DHARMa::testOutliers(m.adfg.iid.tw.resid) 
adfg.iid.tw.resid.quan <- DHARMa::testQuantiles(m.adfg.iid.tw.resid) 
adfg.iid.tw.resid.disp <- DHARMa::testDispersion(m.adfg.iid.tw.resid) 
adfg.iid.tw.resid.qq <- DHARMa::testResiduals(m.adfg.iid.tw.resid)
adfg.iid.tw.resid.zero <- DHARMa::testZeroInflation(m.adfg.iid.tw.resid) 
# TW model + depth
adfg.tw.d.resid.outl <- DHARMa::testOutliers(m.adfg.tw.d.resid) 
adfg.tw.d.resid.quan <- DHARMa::testQuantiles(m.adfg.tw.d.resid) 
adfg.tw.d.resid.disp <- DHARMa::testDispersion(m.adfg.tw.d.resid) 
adfg.tw.d.resid.qq <- DHARMa::testResiduals(m.adfg.tw.d.resid)
adfg.tw.d.resid.zero <- DHARMa::testZeroInflation(m.adfg.tw.d.resid) 
# DG model + depth
adfg.dg.d.resid.outl <- DHARMa::testOutliers(m.adfg.dg.d.resid) 
adfg.dg.d.resid.quan <- DHARMa::testQuantiles(m.adfg.dg.d.resid) 
adfg.dg.d.resid.disp <- DHARMa::testDispersion(m.adfg.dg.d.resid) 
adfg.dg.d.resid.qq <- DHARMa::testResiduals(m.adfg.dg.d.resid)
adfg.dg.d.resid.zero <- DHARMa::testZeroInflation(m.adfg.dg.d.resid) 

# plots of residuals tests
DHARMa::plotQQunif(m.adfg.iid.tw.resid) 
DHARMa::plotResiduals(m.adfg.iid.tw.resid) 
hist(m.adfg.iid.tw.resid)
plot(m.adfg.iid.tw.resid)

png(paste0(plotdir.ns, "/qq_adfg_iid_tw.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.adfg.iid.tw.resid) 
DHARMa::plotResiduals(m.adfg.iid.tw.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_adfg_tw_d.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.adfg.tw.d.resid) 
DHARMa::plotResiduals(m.adfg.tw.d.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_adfg_dg_d.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.adfg.dg.d.resid) 
DHARMa::plotResiduals(m.adfg.dg.d.resid) 
par(mfrow = c(1, 1))
dev.off()

# NOAA model residuals ----

m.noaa.iid.tw.resid <- simulate(update(m.noaa.iid.tw, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.noaa.iid.tw, return_DHARMa = TRUE)
plot(m.noaa.iid.tw.resid)

m.noaa.iid.dg.resid <- simulate(update(m.noaa.iid.dg, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.noaa.iid.dg, return_DHARMa = TRUE)
plot(m.noaa.iid.dg.resid)

m.noaa.iid.dl.resid <- simulate(update(m.noaa.iid.dl, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.noaa.iid.dl, return_DHARMa = TRUE)
plot(m.noaa.iid.dl.resid)

m.noaa.tw.d.resid <- simulate(update(m.noaa.tw.d, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.noaa.tw.d, return_DHARMa = TRUE)
plot(m.noaa.tw.d.resid)

m.noaa.dg.d.resid <- simulate(update(m.noaa.dg.d, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.noaa.dg.d, return_DHARMa = TRUE)
plot(m.noaa.dg.d.resid)

m.noaa.dl.d.resid <- simulate(update(m.noaa.dl.d, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.noaa.dl.d, return_DHARMa = TRUE)
plot(m.noaa.dl.d.resid)

saveRDS(m.noaa.iid.tw.resid, file = paste0(modeldir.ns, "/m_noaa_resid_iid_tw.RDS"))
saveRDS(m.noaa.iid.dg.resid, file = paste0(modeldir.ns, "/m_noaa_resid_iid_dg.RDS"))
saveRDS(m.noaa.iid.dl.resid, file = paste0(modeldir.ns, "/m_noaa_resid_iid_dl.RDS"))
saveRDS(m.noaa.tw.d.resid, file = paste0(modeldir.ns, "/m_noaa_resid_tw_d.RDS"))
saveRDS(m.noaa.dg.d.resid, file = paste0(modeldir.ns, "/m_noaa_resid_dg_d.RDS"))
saveRDS(m.noaa.dl.d.resid, file = paste0(modeldir.ns, "/m_noaa_resid_dl_d.RDS"))
#saveRDS(m.noaa.iid.nb.resid, file = paste0(modeldir.ns, "/m_noaa_resid_iid_nb.RDS"))
#m.noaa.iid.tw.resid <- readRDS(paste0(modeldir.ns, "/m_noaa_resid_iid_tw.RDS"))
#m.noaa.iid.dg.resid <- readRDS(paste0(modeldir.ns, "/m_noaa_resid_iid_dg.RDS"))
#m.noaa.iid.dl.resid <- readRDS(paste0(modeldir.ns, "/m_noaa_resid_iid_dl.RDS"))
#m.noaa.iid.nb.resid <- readRDS(paste0(modeldir.ns, "/m_noaa_resid_iid_nb.RDS"))

# residuals tests
noaa.iid.tw.resid.outl <- DHARMa::testOutliers(m.noaa.iid.tw.resid) 
noaa.iid.tw.resid.quan <- DHARMa::testQuantiles(m.noaa.iid.tw.resid) 
noaa.iid.tw.resid.disp <- DHARMa::testDispersion(m.noaa.iid.tw.resid) 
noaa.iid.tw.resid.qq <- DHARMa::testResiduals(m.noaa.iid.tw.resid)
noaa.iid.tw.resid.zero <- DHARMa::testZeroInflation(m.noaa.iid.tw.resid) 
noaa.iid.dg.resid.outl <- DHARMa::testOutliers(m.noaa.iid.dg.resid) 
noaa.iid.dg.resid.quan <- DHARMa::testQuantiles(m.noaa.iid.dg.resid) 
noaa.iid.dg.resid.disp <- DHARMa::testDispersion(m.noaa.iid.dg.resid) 
noaa.iid.dg.resid.qq <- DHARMa::testResiduals(m.noaa.iid.dg.resid)
noaa.iid.dg.resid.zero <- DHARMa::testZeroInflation(m.noaa.iid.dg.resid) 
noaa.iid.dl.resid.outl <- DHARMa::testOutliers(m.noaa.iid.dl.resid) 
noaa.iid.dl.resid.quan <- DHARMa::testQuantiles(m.noaa.iid.dl.resid) 
noaa.iid.dl.resid.disp <- DHARMa::testDispersion(m.noaa.iid.dl.resid) 
noaa.iid.dl.resid.qq <- DHARMa::testResiduals(m.noaa.iid.dl.resid)
noaa.iid.dl.resid.zero <- DHARMa::testZeroInflation(m.noaa.iid.dl.resid) 
noaa.tw.d.resid.outl <- DHARMa::testOutliers(m.noaa.tw.d.resid) 
noaa.tw.d.resid.quan <- DHARMa::testQuantiles(m.noaa.tw.d.resid) 
noaa.tw.d.resid.disp <- DHARMa::testDispersion(m.noaa.tw.d.resid) 
noaa.tw.d.resid.qq <- DHARMa::testResiduals(m.noaa.tw.d.resid)
noaa.tw.d.resid.zero <- DHARMa::testZeroInflation(m.noaa.tw.d.resid) 
noaa.dg.d.resid.outl <- DHARMa::testOutliers(m.noaa.dg.d.resid) 
noaa.dg.d.resid.quan <- DHARMa::testQuantiles(m.noaa.dg.d.resid) 
noaa.dg.d.resid.disp <- DHARMa::testDispersion(m.noaa.dg.d.resid) 
noaa.dg.d.resid.qq <- DHARMa::testResiduals(m.noaa.dg.d.resid)
noaa.dg.d.resid.zero <- DHARMa::testZeroInflation(m.noaa.dg.d.resid) 
noaa.dl.d.resid.outl <- DHARMa::testOutliers(m.noaa.dl.d.resid) 
noaa.dl.d.resid.quan <- DHARMa::testQuantiles(m.noaa.dl.d.resid) 
noaa.dl.d.resid.disp <- DHARMa::testDispersion(m.noaa.dl.d.resid) 
noaa.dl.d.resid.qq <- DHARMa::testResiduals(m.noaa.dl.d.resid)
noaa.dl.d.resid.zero <- DHARMa::testZeroInflation(m.noaa.dl.d.resid) 

# plots of residuals tests
png(paste0(plotdir.ns, "/qq_noaa_iid_tw.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.noaa.iid.tw.resid) 
DHARMa::plotResiduals(m.noaa.iid.tw.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_noaa_iid_dg.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.noaa.iid.dg.resid) 
DHARMa::plotResiduals(m.noaa.iid.dg.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_noaa_iid_dl.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.noaa.iid.dl.resid) 
DHARMa::plotResiduals(m.noaa.iid.dl.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_noaa_tw_d.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.noaa.tw.d.resid) 
DHARMa::plotResiduals(m.noaa.tw.d.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_noaa_dg_d.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.noaa.dg.d.resid) 
DHARMa::plotResiduals(m.noaa.dg.d.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_noaa_dl_d.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.noaa.dl.d.resid) 
DHARMa::plotResiduals(m.noaa.dl.d.resid) 
par(mfrow = c(1, 1))
dev.off()

# NBS model residuals ----

m.nbs.iid.tw.resid <- simulate(update(m.nbs.iid.tw, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nbs.iid.tw, return_DHARMa = TRUE)
plot(m.nbs.iid.tw.resid)

m.nbs.tw.d.resid <- simulate(update(m.nbs.tw.d, family = delta_lognormal()), nsim = 100, type = "mle-mvn") |>
  dharma_residuals(m.nbs.tw.d, return_DHARMa = TRUE)
plot(m.nbs.tw.d.resid)

saveRDS(m.nbs.iid.tw.resid, file = paste0(modeldir.ns, "/m_nbs_resid_iid_tw.RDS"))
saveRDS(m.nbs.tw.d.resid, file = paste0(modeldir.ns, "/m_nbs_resid_tw_d.RDS"))
#m.nbs.iid.tw.resid <- readRDS(paste0(modeldir.ns, "/m_nbs_resid_iid_tw.RDS"))
#m.nbs.iid.nb.resid <- readRDS(paste0(modeldir.ns, "/m_nbs_resid_iid_nb.RDS"))

# residuals tests
nbs.iid.tw.resid.outl <- DHARMa::testOutliers(m.nbs.iid.tw.resid) 
nbs.iid.tw.resid.quan <- DHARMa::testQuantiles(m.nbs.iid.tw.resid) 
nbs.iid.tw.resid.disp <- DHARMa::testDispersion(m.nbs.iid.tw.resid) 
nbs.iid.tw.resid.qq <- DHARMa::testResiduals(m.nbs.iid.tw.resid)
nbs.iid.tw.resid.zero <- DHARMa::testZeroInflation(m.nbs.iid.tw.resid) 
nbs.tw.d.resid.outl <- DHARMa::testOutliers(m.nbs.tw.d.resid) 
nbs.tw.d.resid.quan <- DHARMa::testQuantiles(m.nbs.tw.d.resid) 
nbs.tw.d.resid.disp <- DHARMa::testDispersion(m.nbs.tw.d.resid) 
nbs.tw.d.resid.qq <- DHARMa::testResiduals(m.nbs.tw.d.resid)
nbs.tw.d.resid.zero <- DHARMa::testZeroInflation(m.nbs.tw.d.resid) 

# plots of residuals tests
png(paste0(plotdir.ns, "/qq_nbs_iid_tw.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.nbs.iid.tw.resid) 
DHARMa::plotResiduals(m.nbs.iid.tw.resid) 
par(mfrow = c(1, 1))
dev.off()

png(paste0(plotdir.ns, "/qq_nbs_tw_d.png"), width = 7, height = 4, units = "in", res = 720)
par(mfrow = c(1, 2)) 
DHARMa::plotQQunif(m.nbs.tw.d.resid) 
DHARMa::plotResiduals(m.nbs.tw.d.resid) 
par(mfrow = c(1, 1))
dev.off()

# make table with loglikelihood and residuals tests results ----
compare.models <- data.frame(Survey=character(),
                                  model=character(),
                                  Effects=character(),
                                  loglik=numeric(),
                                  ks = numeric(), quantile = numeric(), 
                                  dispersion = numeric(), outliers = numeric(), 
                                  zero_inf = numeric(), stringsAsFactors=FALSE) %>%
  # ADFG survey model
  add_row(Survey = "ADF&G", model = "Tweedie", Effects="year", loglik = m.adfg.iid.tw.cv$sum_loglik, ks = adfg.iid.tw.resid.qq$uniformity$p.value, quantile = adfg.iid.tw.resid.quan$p.value, dispersion = adfg.iid.tw.resid.disp$p.value, outliers = adfg.iid.tw.resid.outl$p.value, zero_inf = adfg.iid.tw.resid.zero$p.value) %>%
  add_row(Survey = "ADF&G", model = "Tweedie", Effects="year, depth", loglik = m.adfg.tw.d.cv$sum_loglik, ks = adfg.tw.d.resid.qq$uniformity$p.value, quantile = adfg.tw.d.resid.quan$p.value, dispersion = adfg.tw.d.resid.disp$p.value, outliers = adfg.tw.d.resid.outl$p.value, zero_inf = adfg.tw.d.resid.zero$p.value) %>%
  add_row(Survey = "ADF&G", model = "DG", Effects="year, depth", loglik = m.adfg.dg.d.cv$sum_loglik, ks = adfg.dg.d.resid.qq$uniformity$p.value, quantile = adfg.dg.d.resid.quan$p.value, dispersion = adfg.dg.d.resid.disp$p.value, outliers = adfg.dg.d.resid.outl$p.value, zero_inf = adfg.dg.d.resid.zero$p.value) %>%
  # NOAA survey models
  add_row(Survey = "NOAA NS", model = "Tweedie", Effects="year", loglik = m.noaa.iid.tw.cv$sum_loglik, ks = noaa.iid.tw.resid.qq$uniformity$p.value, quantile = noaa.iid.tw.resid.quan$p.value, dispersion = noaa.iid.tw.resid.disp$p.value, outliers = noaa.iid.tw.resid.outl$p.value, zero_inf = noaa.iid.tw.resid.zero$p.value) %>%
  add_row(Survey = "NOAA NS", model = "DG", Effects="year", loglik = m.noaa.iid.dg.cv$sum_loglik, ks = noaa.iid.dg.resid.qq$uniformity$p.value, quantile = noaa.iid.dg.resid.quan$p.value, dispersion = noaa.iid.dg.resid.disp$p.value, outliers = noaa.iid.dg.resid.outl$p.value, zero_inf = noaa.iid.dg.resid.zero$p.value) %>%
  add_row(Survey = "NOAA NS", model = "DL", Effects="year", loglik = m.noaa.iid.dl.cv$sum_loglik, ks = noaa.iid.dl.resid.qq$uniformity$p.value, quantile = noaa.iid.dl.resid.quan$p.value, dispersion = noaa.iid.dl.resid.disp$p.value, outliers = noaa.iid.dl.resid.outl$p.value, zero_inf = noaa.iid.dl.resid.zero$p.value) %>%
  add_row(Survey = "NOAA NS", model = "Tweedie", Effects="year, depth", loglik = m.noaa.tw.d.cv$sum_loglik, ks = noaa.tw.d.resid.qq$uniformity$p.value, quantile = noaa.tw.d.resid.quan$p.value, dispersion = noaa.tw.d.resid.disp$p.value, outliers = noaa.tw.d.resid.outl$p.value, zero_inf = noaa.tw.d.resid.zero$p.value) %>%
  add_row(Survey = "NOAA NS", model = "DG", Effects="year, depth", loglik = m.noaa.dg.d.cv$sum_loglik, ks = noaa.dg.d.resid.qq$uniformity$p.value, quantile = noaa.dg.d.resid.quan$p.value, dispersion = noaa.dg.d.resid.disp$p.value, outliers = noaa.dg.d.resid.outl$p.value, zero_inf = noaa.dg.d.resid.zero$p.value) %>%
  add_row(Survey = "NOAA NS", model = "DL", Effects="year, depth", loglik = m.noaa.dl.d.cv$sum_loglik, ks = noaa.dl.d.resid.qq$uniformity$p.value, quantile = noaa.dl.d.resid.quan$p.value, dispersion = noaa.dl.d.resid.disp$p.value, outliers = noaa.dl.d.resid.outl$p.value, zero_inf = noaa.dl.d.resid.zero$p.value) %>%
  #add_row(Survey = "NOAA NS", model = "delta lognormal", loglik = m.noaa.iid.nb.cv$sum_loglik, ks = noaa.iid.nb.resid.qq$uniformity$p.value, quantile = noaa.iid.nb.resid.quan$p.value, dispersion = noaa.iid.nb.resid.disp$p.value, outliers = noaa.iid.nb.resid.outl$p.value, zero_inf = noaa.iid.nb.resid.zero$p.value) %>%
  # NBS survey models
  add_row(Survey = "NOAA NBS", model = "Tweedie", Effects="year", loglik = m.nbs.iid.tw.cv$sum_loglik, ks = nbs.iid.tw.resid.qq$uniformity$p.value, quantile = nbs.iid.tw.resid.quan$p.value, dispersion = nbs.iid.tw.resid.disp$p.value, outliers = nbs.iid.tw.resid.outl$p.value, zero_inf = nbs.iid.tw.resid.zero$p.value) %>%
  add_row(Survey = "NOAA NBS", model = "Tweedie", Effects="year, depth", loglik = m.nbs.tw.d.cv$sum_loglik, ks = nbs.tw.d.resid.qq$uniformity$p.value, quantile = nbs.tw.d.resid.quan$p.value, dispersion = nbs.tw.d.resid.disp$p.value, outliers = nbs.tw.d.resid.outl$p.value, zero_inf = nbs.tw.d.resid.zero$p.value) %>%
  #add_row(Survey = "NOAA NBS", model = "NB2", loglik = m.nbs.iid.nb.cv$sum_loglik, ks = nbs.iid.nb.resid.qq$uniformity$p.value, quantile = nbs.iid.nb.resid.quan$p.value, dispersion = nbs.iid.nb.resid.disp$p.value, outliers = nbs.iid.nb.resid.outl$p.value, zero_inf = nbs.iid.nb.resid.zero$p.value) %>%
  # format display of numbers
  mutate(across(loglik, ~ format(., big.mark = ",", scientific = F, digits = 1))) %>%
  mutate(across(ks:zero_inf, ~ format(., scientific = F, digits = 2))) %>%
  mutate(quantile = "<0.01") %>%
  mutate(dispersion = case_when(
    dispersion == "0.00" ~ "<0.01",
    TRUE ~ dispersion
  ))
  
# save table with predictive skill values and residuals tests results
write.csv(compare.models, paste0(here::here(), "/NSRKC/output/ns_compare_models.csv"))


# spatial residuals plots ----
# add the scaled residuals to data frame and then make plots

# spatial residuals plotting function
spat_resid_plot <- function(dat, survey, model) {
  ggplot(dat) + 
    geom_point(aes(y = Latitude, x = Longitude, color = scaledResiduals), size = 1) +
    scale_color_gradient2(midpoint = 0) + 
    labs(y = "Latitude",
         x = "Longitude") +
    theme_gray() + 
    ggtitle(paste0(survey, " trawl survey, " , model))+
    facet_wrap(~Year)+
    theme(axis.title = element_text(size = 10),
          legend.position = "bottom") +
    labs(color = "scaled residuals") +
    scale_x_continuous(limits = c(-169, -161), breaks=seq(-169, -161, by = 6)) +
    scale_y_continuous(limits = c(62.8, 66), breaks=seq(62.8, 66, by = 1))
}

# add ADF&G model spatial residuals to data 
spres_adfg_tw <- cbind(nsrkc_utm_adfg, 
                               data.frame(m.adfg.iid.tw.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.adfg.iid.tw.resid.scaledResiduals)
spres_adfg_tw_d <- cbind(nsrkc_utm_adfg_d, 
                       data.frame(m.adfg.tw.d.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.adfg.tw.d.resid.scaledResiduals)
spres_adfg_dg_d <- cbind(nsrkc_utm_adfg_d, 
                       data.frame(m.adfg.dg.d.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.adfg.dg.d.resid.scaledResiduals)

# ADF&G spatial residuals plots
spres.adfg.tw.plot <- spat_resid_plot(spres_adfg_tw, "ADF&G", "Tweedie, year")
spres.adfg.tw.d.plot <- spat_resid_plot(spres_adfg_tw_d, "ADF&G", "Tweedie, year + depth")
spres.adfg.dg.d.plot <- spat_resid_plot(spres_adfg_dg_d, "ADF&G", "DG, year + depth")

# save ADF&G spatial residuals plots
ggsave(paste0(plotdir.ns, "/spat_resid_adfg_iid_tw.png"), spres.adfg.tw.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_adfg_tw_d.png"), spres.adfg.tw.d.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_adfg_dg_d.png"), spres.adfg.dg.d.plot, height = 5, width = 7, units = "in")

# NOAA model spatial residuals
spres_noaa_tw <- cbind(nsrkc_utm_noaa, 
                      data.frame(m.noaa.iid.tw.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.noaa.iid.tw.resid.scaledResiduals)
spres_noaa_dg <- cbind(nsrkc_utm_noaa, 
                       data.frame(m.noaa.iid.dg.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.noaa.iid.dg.resid.scaledResiduals)
spres_noaa_dl <- cbind(nsrkc_utm_noaa, 
                       data.frame(m.noaa.iid.dl.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.noaa.iid.dl.resid.scaledResiduals)
spres_noaa_tw_d <- cbind(nsrkc_utm_noaa_d, 
                       data.frame(m.noaa.tw.d.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.noaa.tw.d.resid.scaledResiduals)
spres_noaa_dg_d <- cbind(nsrkc_utm_noaa_d, 
                       data.frame(m.noaa.dg.d.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.noaa.dg.d.resid.scaledResiduals)
spres_noaa_dl_d <- cbind(nsrkc_utm_noaa_d, 
                       data.frame(m.noaa.dl.d.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.noaa.dl.d.resid.scaledResiduals)
  
# NOAA model spatial residuals plots
spres.noaa.tw.plot <- spat_resid_plot(spres_noaa_tw, "NOAA NS", "Tweedie, year")
spres.noaa.dg.plot <- spat_resid_plot(spres_noaa_dg, "NOAA NS", "DG, year")
spres.noaa.dl.plot <- spat_resid_plot(spres_noaa_dl, "NOAA NS", "DL, year")
spres.noaa.tw.d.plot <- spat_resid_plot(spres_noaa_tw_d, "NOAA NS", "Tweedie, year + depth")
spres.noaa.dg.d.plot <- spat_resid_plot(spres_noaa_dg_d, "NOAA NS", "DG, year + depth")
spres.noaa.dl.d.plot <- spat_resid_plot(spres_noaa_dl_d, "NOAA NS", "DL, year + depth")
  
# save NOAA spatial residuals plots  
ggsave(paste0(plotdir.ns, "/spat_resid_noaa_iid_tw.png"), spres.noaa.tw.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_noaa_iid_dg.png"), spres.noaa.dg.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_noaa_iid_dl.png"), spres.noaa.dl.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_noaa_tw_d.png"), spres.noaa.tw.d.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_noaa_dg_d.png"), spres.noaa.dg.d.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_noaa_dl_d.png"), spres.noaa.dl.d.plot, height = 5, width = 7, units = "in")


# NBS model spatial residuals
spres_nbs_tw <- cbind(nsrkc_utm_nbs, 
                       data.frame(m.nbs.iid.tw.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.nbs.iid.tw.resid.scaledResiduals)
spres_nbs_tw_d <- cbind(nsrkc_utm_nbs_d, 
                      data.frame(m.nbs.tw.d.resid$scaledResiduals)) %>%
  rename("scaledResiduals" = m.nbs.tw.d.resid.scaledResiduals)

spres.nbs.tw.plot <- spat_resid_plot(spres_nbs_tw, "NOAA NBS", "Tweedie, year")
spres.nbs.tw.d.plot <- spat_resid_plot(spres_nbs_tw_d, "NOAA NBS", "Tweedie, year + depth")

ggsave(paste0(plotdir.ns, "/spat_resid_nbs_iid_tw.png"), spres.nbs.tw.plot, height = 5, width = 7, units = "in")
ggsave(paste0(plotdir.ns, "/spat_resid_nbs_tw_d.png"), spres.nbs.tw.d.plot, height = 5, width = 7, units = "in")


# Moran's I analysis ----
# this tests for spatial autocorrelation across the spatial domain in a given year
# p-value < 0.05 provides evidence that we should reject the null hypothesis
# the null hypothesis is that the scaled residual values are randomly distributed in space
# significance in Moran's I indicates spatial autocorrelation when the score is positive
# and negative spatial autocorrelation when the score is negative

# check for skew and outliers in distribution of scaled residuals, since Moran's I statistic is not robust to those
hist(spres_adfg_tw$scaledResiduals)
hist(spres_adfg_tw_d$scaledResiduals)
hist(spres_noaa_tw$scaledResiduals)
hist(spres_noaa_dg$scaledResiduals)
hist(spres_noaa_dl$scaledResiduals)
hist(spres_noaa_tw_d$scaledResiduals)
hist(spres_noaa_dg_d$scaledResiduals)
hist(spres_noaa_dl_d$scaledResiduals)
hist(spres_nbs_tw$scaledResiduals)
hist(spres_nbs_tw_d$scaledResiduals)

# use the 10 nearest neighbors to calculate weights
# extract longitude and latitude from the data
# within a year
spres_adfg_2024 <- spres_adfg_tw %>% filter(Year == 2024)
longlats <- cbind(long = spres_adfg_2024$Longitude, lat = spres_adfg_2024$Latitude) %>% as.data.frame()
nb_list <- knn2nb(knearneigh(longlats, k=10, longlat = TRUE, use_kd_tree=FALSE))
nb_weights <- nb2listw(nb_list, style="W", zero.policy=TRUE)
spdep::moran.test(spres_adfg_2024$scaledResiduals, nb_weights, alternative = "greater")
moran.mc.adfg <- spdep::moran.mc(spres_adfg_2024$scaledResiduals, nb_weights, alternative = "greater", nsim = 1000)
plot(moran.mc.adfg)
moran.plot(spres_adfg_2024$scaledResiduals, nb_weights)

# function to calculate Moran's I statistic for each year
moran.i.yr <- function(dat, year, survey, model){
  dat.yr <- dat %>% filter(Year == year)
  longlats <- cbind(long = dat.yr$Longitude, lat = dat.yr$Latitude) %>% as.data.frame()
  nb_list <- knn2nb(knearneigh(longlats, k=10, longlat = TRUE, use_kd_tree=FALSE))
  nb_weights <- nb2listw(nb_list, style="W", zero.policy=TRUE)
  test.results <- spdep::moran.mc(dat.yr$scaledResiduals, nb_weights, alternative = "greater", nsim = 1000)
  output <- data.frame(survey, model, year, test.results$statistic, test.results$p.value)
  rownames(output) <- NULL
  colnames(output) <- c("survey", "model", "year", "statistic", "p.value")
  output
}

moran.adfg.2024 <- moran.i.yr(spres_adfg_tw, 2024, "ADF&G", "Tweedie")

# years of the survey
adfg.years <- as.vector(unique(spres_adfg_tw$Year))
noaa.years <- as.vector(unique(spres_noaa_tw$Year))
nbs.years <- as.vector(unique(spres_nbs_tw$Year))

# apply moran.i.yr function across all years of the survey for each model
# ADFG model
moran.i.adfg.tw <- purrr::map_dfr(adfg.years, moran.i.yr, dat = spres_adfg_tw, survey = "ADF&G", model = "Tweedie, year")
moran.i.adfg.tw.d <- purrr::map_dfr(adfg.years, moran.i.yr, dat = spres_adfg_tw_d, survey = "ADF&G", model = "Tweedie, year + depth")
moran.i.adfg.dg.d <- purrr::map_dfr(adfg.years, moran.i.yr, dat = spres_adfg_dg_d, survey = "ADF&G", model = "DG, year + depth")
moran.i.noaa.tw <- purrr::map_dfr(noaa.years, moran.i.yr, dat = spres_noaa_tw, survey = "NOAA NS", model = "Tweedie, year")
moran.i.noaa.dg <- purrr::map_dfr(noaa.years, moran.i.yr, dat = spres_noaa_dg, survey = "NOAA NS", model = "DG, year")
moran.i.noaa.dl <- purrr::map_dfr(noaa.years, moran.i.yr, dat = spres_noaa_dl, survey = "NOAA NS", model = "DL, year")
moran.i.noaa.tw.d <- purrr::map_dfr(noaa.years, moran.i.yr, dat = spres_noaa_tw_d, survey = "NOAA NS", model = "Tweedie, year + depth")
moran.i.noaa.dg.d <- purrr::map_dfr(noaa.years, moran.i.yr, dat = spres_noaa_dg_d, survey = "NOAA NS", model = "DG, year + depth")
moran.i.noaa.dl.d <- purrr::map_dfr(noaa.years, moran.i.yr, dat = spres_noaa_dl_d, survey = "NOAA NS", model = "DL, year + depth")
moran.i.nbs.tw <- purrr::map_dfr(nbs.years, moran.i.yr, dat = spres_nbs_tw, survey = "NOAA NBS", model = "Tweedie, year")
moran.i.nbs.tw.d <- purrr::map_dfr(nbs.years, moran.i.yr, dat = spres_nbs_tw_d, survey = "NOAA NBS", model = "Tweedie, year + depth")

moran.i.all <- rbind(moran.i.adfg.tw, moran.i.adfg.tw.d, moran.i.adfg.dg.d, moran.i.noaa.tw, moran.i.noaa.dg, moran.i.noaa.dl, moran.i.noaa.tw.d, moran.i.noaa.dg.d, moran.i.noaa.dl.d, moran.i.nbs.tw, moran.i.nbs.tw.d)

write.csv(moran.i.all, paste0(here::here(), "/NSRKC/output/moran_i_all.csv"))

# ************************************************************************************************
# generate predictions ----
# *************************************************************************************************

# replicate grid across all years
# full survey area grids
predgrid.utm.adfg <- replicate_df(pred_grid_ns_km, "Year", unique(nsrkc_utm_adfg$Year))
predgrid.utm.adfg$year_f <- factor(predgrid.utm.adfg$Year)
dplyr::glimpse(predgrid.utm.adfg)

predgrid.utm.noaa <- replicate_df(pred_grid_ns_km, "Year", unique(nsrkc_utm_noaa$Year))
predgrid.utm.noaa$year_f <- factor(predgrid.utm.noaa$Year)
dplyr::glimpse(predgrid.utm.noaa)

predgrid.utm.nbs <- replicate_df(pred_grid_nbs_km, "Year", unique(nsrkc_utm_nbs$Year))
predgrid.utm.nbs$year_f <- factor(predgrid.utm.nbs$Year)
dplyr::glimpse(predgrid.utm.nbs)

# reduced area grids
predgrid.utm.adfg.s <- replicate_df(pred_grid_ns_km.s, "Year", unique(nsrkc_utm_adfg$Year))
predgrid.utm.adfg.s$year_f <- factor(predgrid.utm.adfg.s$Year)
dplyr::glimpse(predgrid.utm.adfg.s)

predgrid.utm.noaa.s <- replicate_df(pred_grid_ns_km.s, "Year", unique(nsrkc_utm_noaa$Year))
predgrid.utm.noaa.s$year_f <- factor(predgrid.utm.noaa.s$Year)
dplyr::glimpse(predgrid.utm.noaa.s)

predgrid.utm.nbs.s <- replicate_df(pred_grid_nbs_km.s, "Year", unique(nsrkc_utm_nbs$Year))
predgrid.utm.nbs.s$year_f <- factor(predgrid.utm.nbs.s$Year)
dplyr::glimpse(predgrid.utm.nbs.s)

# full survey area grids with depth
predgrid.utm.adfg.d <- replicate_df(pred_grid_ns_km_dep, "Year", unique(nsrkc_utm_adfg_d$Year))
predgrid.utm.adfg.d$year_f <- factor(predgrid.utm.adfg.d$Year)

predgrid.utm.noaa.d <- replicate_df(pred_grid_ns_km_dep, "Year", unique(nsrkc_utm_noaa_d$Year))
predgrid.utm.noaa.d$year_f <- factor(predgrid.utm.noaa.d$Year)

predgrid.utm.nbs.d <- replicate_df(pred_grid_nbs_km_dep, "Year", unique(nsrkc_utm_nbs_d$Year))
predgrid.utm.nbs.d$year_f <- factor(predgrid.utm.nbs.d$Year)

# reduced area grids with depth
predgrid.utm.adfg.s.d <- replicate_df(pred_grid_ns_km.s_dep, "Year", unique(nsrkc_utm_adfg_d$Year))
predgrid.utm.adfg.s.d$year_f <- factor(predgrid.utm.adfg.s.d$Year)
dplyr::glimpse(predgrid.utm.adfg.s.d)

predgrid.utm.noaa.s.d <- replicate_df(pred_grid_ns_km.s_dep, "Year", unique(nsrkc_utm_noaa_d$Year))
predgrid.utm.noaa.s.d$year_f <- factor(predgrid.utm.noaa.s.d$Year)
dplyr::glimpse(predgrid.utm.noaa.s.d)

predgrid.utm.nbs.s.d <- replicate_df(pred_grid_nbs_km.s_dep, "Year", unique(nsrkc_utm_nbs_d$Year))
predgrid.utm.nbs.s.d$year_f <- factor(predgrid.utm.nbs.s.d$Year)
dplyr::glimpse(predgrid.utm.nbs.s.d)

# predictions on new data ----
# predict over both full area grid and reduced area grid for each model
# ADFG model predictions
pred.adfg.iid.tw <- predict(m.adfg.iid.tw, newdata = predgrid.utm.adfg, return_tmb_object = T)
pred.adfg.iid.tw.s <- predict(m.adfg.iid.tw, newdata = predgrid.utm.adfg.s, return_tmb_object = T)
pred.adfg.tw.d <- predict(m.adfg.tw.d, newdata = predgrid.utm.adfg.d, return_tmb_object = T)
pred.adfg.tw.s.d <- predict(m.adfg.tw.d, newdata = predgrid.utm.adfg.s.d, return_tmb_object = T)
pred.adfg.dg.d <- predict(m.adfg.dg.d, newdata = predgrid.utm.adfg.d, return_tmb_object = T)
pred.adfg.dg.s.d <- predict(m.adfg.dg.d, newdata = predgrid.utm.adfg.s.d, return_tmb_object = T)
# NOAA model predictions
pred.noaa.iid.tw <- predict(m.noaa.iid.tw, newdata = predgrid.utm.noaa, return_tmb_object = T)
pred.noaa.iid.tw.s <- predict(m.noaa.iid.tw, newdata = predgrid.utm.noaa.s, return_tmb_object = T)
pred.noaa.iid.dg <- predict(m.noaa.iid.dg, newdata = predgrid.utm.noaa, return_tmb_object = T)
pred.noaa.iid.dg.s <- predict(m.noaa.iid.dg, newdata = predgrid.utm.noaa.s, return_tmb_object = T)
pred.noaa.iid.dl <- predict(m.noaa.iid.dl, newdata = predgrid.utm.noaa, return_tmb_object = T)
pred.noaa.iid.dl.s <- predict(m.noaa.iid.dl, newdata = predgrid.utm.noaa.s, return_tmb_object = T)
pred.noaa.iid.tw.d <- predict(m.noaa.tw.d, newdata = predgrid.utm.noaa.d, return_tmb_object = T)
pred.noaa.iid.tw.s.d <- predict(m.noaa.tw.d, newdata = predgrid.utm.noaa.s.d, return_tmb_object = T)
pred.noaa.iid.dg.d <- predict(m.noaa.dg.d, newdata = predgrid.utm.noaa.d, return_tmb_object = T)
pred.noaa.iid.dg.s.d <- predict(m.noaa.dg.d, newdata = predgrid.utm.noaa.s.d, return_tmb_object = T)
pred.noaa.iid.dl.d <- predict(m.noaa.dl.d, newdata = predgrid.utm.noaa.d, return_tmb_object = T)
pred.noaa.iid.dl.s.d <- predict(m.noaa.dl.d, newdata = predgrid.utm.noaa.s.d, return_tmb_object = T)
#pred.noaa.iid.nb <- predict(m.noaa.iid.nb, newdata = predgrid.utm.noaa, return_tmb_object = T)
#pred.noaa.iid.nb.s <- predict(m.noaa.iid.nb, newdata = predgrid.utm.noaa.s, return_tmb_object = T)
# NBS model predictions
pred.nbs.iid.tw <- predict(m.nbs.iid.tw, newdata = predgrid.utm.nbs, return_tmb_object = T)
pred.nbs.iid.tw.s <- predict(m.nbs.iid.tw, newdata = predgrid.utm.nbs.s, return_tmb_object = T)
pred.nbs.iid.tw.d <- predict(m.nbs.tw.d, newdata = predgrid.utm.nbs.d, return_tmb_object = T)
pred.nbs.iid.tw.s.d <- predict(m.nbs.tw.d, newdata = predgrid.utm.nbs.s.d, return_tmb_object = T)
#pred.nbs.iid.nb <- predict(m.nbs.iid.nb, newdata = predgrid.utm.nbs, return_tmb_object = T)
#pred.nbs.iid.nb.s <- predict(m.nbs.iid.nb, newdata = predgrid.utm.nbs.s, return_tmb_object = T)

# save predictions
saveRDS(pred.adfg.iid.tw, file = paste0(modeldir.ns, "/m_adfg_pred_iid_tw.RDS"))
saveRDS(pred.adfg.iid.tw.s, file = paste0(modeldir.ns, "/m_adfg_pred_iid_tw_s.RDS"))
saveRDS(pred.noaa.iid.tw, file = paste0(modeldir.ns, "/m_noaa_pred_iid_tw.RDS"))
saveRDS(pred.noaa.iid.tw.s, file = paste0(modeldir.ns, "/m_noaa_pred_iid_tw_s.RDS"))
saveRDS(pred.noaa.iid.dg, file = paste0(modeldir.ns, "/m_noaa_pred_iid_dg.RDS"))
saveRDS(pred.noaa.iid.dg.s, file = paste0(modeldir.ns, "/m_noaa_pred_iid_dg_s.RDS"))
saveRDS(pred.noaa.iid.dl, file = paste0(modeldir.ns, "/m_noaa_pred_iid_dl.RDS"))
saveRDS(pred.noaa.iid.dl.s, file = paste0(modeldir.ns, "/m_noaa_pred_iid_dl_s.RDS"))
#saveRDS(pred.noaa.iid.nb, file = paste0(modeldir.ns, "/m_noaa_pred_iid_nb.RDS"))
#saveRDS(pred.noaa.iid.nb.s, file = paste0(modeldir.ns, "/m_noaa_pred_iid_nb_s.RDS"))
saveRDS(pred.nbs.iid.tw, file = paste0(modeldir.ns, "/m_nbs_pred_iid_tw.RDS"))
saveRDS(pred.nbs.iid.tw.s, file = paste0(modeldir.ns, "/m_nbs_pred_iid_tw_s.RDS"))
#saveRDS(pred.nbs.iid.nb, file = paste0(modeldir.ns, "/m_nbs_pred_iid_nb.RDS"))
#saveRDS(pred.nbs.iid.nb.s, file = paste0(modeldir.ns, "/m_nbs_pred_iid_nb_s.RDS"))

# bring in saved predictions
#pred.adfg.iid.tw <- readRDS(paste0(modeldir.ns, "/m_adfg_pred_iid_tw.RDS"))
#pred.adfg.iid.tw.s <- readRDS(paste0(modeldir.ns, "/m_adfg_pred_iid_tw_s.RDS"))

# estimate the uncertainty in spatiotemporal density predictions using simulations ----
# from the joint precision matrix 

# ADFG sims
pred.adfg.iid.tw.sim <- predict(m.adfg.iid.tw, newdata = predgrid.utm.adfg, nsim = 100)
pred.adfg.iid.tw.s.sim <- predict(m.adfg.iid.tw, newdata = predgrid.utm.adfg.s, nsim = 100)
pred.adfg.iid.tw.d.sim <- predict(m.adfg.tw.d, newdata = predgrid.utm.adfg.d, nsim = 100)
pred.adfg.iid.tw.s.d.sim <- predict(m.adfg.tw.d, newdata = predgrid.utm.adfg.s.d, nsim = 100)
pred.adfg.iid.dg.d.sim <- predict(m.adfg.dg.d, newdata = predgrid.utm.adfg.d, nsim = 100)
pred.adfg.iid.dg.s.d.sim <- predict(m.adfg.dg.d, newdata = predgrid.utm.adfg.s.d, nsim = 100)
# NOAA sims
pred.noaa.iid.tw.sim <- predict(m.noaa.iid.tw, newdata = predgrid.utm.noaa, nsim = 100)
pred.noaa.iid.tw.s.sim <- predict(m.noaa.iid.tw, newdata = predgrid.utm.noaa.s, nsim = 100)
pred.noaa.iid.dg.sim <- predict(m.noaa.iid.dg, newdata = predgrid.utm.noaa, nsim = 100)
pred.noaa.iid.dg.s.sim <- predict(m.noaa.iid.dg, newdata = predgrid.utm.noaa.s, nsim = 100)
pred.noaa.iid.dl.sim <- predict(m.noaa.iid.dl, newdata = predgrid.utm.noaa, nsim = 100)
pred.noaa.iid.dl.s.sim <- predict(m.noaa.iid.dl, newdata = predgrid.utm.noaa.s, nsim = 100)
pred.noaa.iid.tw.d.sim <- predict(m.noaa.tw.d, newdata = predgrid.utm.noaa.d, nsim = 100)
pred.noaa.iid.tw.s.d.sim <- predict(m.noaa.tw.d, newdata = predgrid.utm.noaa.s.d, nsim = 100)
pred.noaa.iid.dg.d.sim <- predict(m.noaa.dg.d, newdata = predgrid.utm.noaa.d, nsim = 100)
pred.noaa.iid.dg.s.d.sim <- predict(m.noaa.dg.d, newdata = predgrid.utm.noaa.s.d, nsim = 100)
pred.noaa.iid.dl.d.sim <- predict(m.noaa.dl.d, newdata = predgrid.utm.noaa.d, nsim = 100)
pred.noaa.iid.dl.s.d.sim <- predict(m.noaa.dl.d, newdata = predgrid.utm.noaa.s.d, nsim = 100)
#pred.noaa.iid.nb.sim <- predict(m.noaa.iid.nb, newdata = predgrid.utm.noaa, nsim = 100)
#pred.noaa.iid.nb.s.sim <- predict(m.noaa.iid.nb, newdata = predgrid.utm.noaa.s, nsim = 100)
# NBS sims
pred.nbs.iid.tw.sim <- predict(m.nbs.iid.tw, newdata = predgrid.utm.nbs, nsim = 100)
pred.nbs.iid.tw.s.sim <- predict(m.nbs.iid.tw, newdata = predgrid.utm.nbs.s, nsim = 100)
pred.nbs.iid.tw.d.sim <- predict(m.nbs.tw.d, newdata = predgrid.utm.nbs.d, nsim = 100)
pred.nbs.iid.tw.s.d.sim <- predict(m.nbs.tw.d, newdata = predgrid.utm.nbs.s.d, nsim = 100)
#pred.nbs.iid.nb.sim <- predict(m.nbs.iid.nb, newdata = predgrid.utm.nbs, nsim = 100)
#pred.nbs.iid.nb.s.sim <- predict(m.nbs.iid.nb, newdata = predgrid.utm.nbs.s, nsim = 100)

# from vignette
sim_last <- pred.adfg.iid.tw.sim[predgrid.utm.adfg$Year == max(predgrid.utm.adfg$Year), ] # just plot last year
pred_last <- pred.adfg.iid.tw$data[pred.adfg.iid.tw$data$Year == max(predgrid.utm.adfg$Year), ]
pred_last$lwr <- apply(exp(sim_last), 1, quantile, probs = 0.025)
pred_last$upr <- apply(exp(sim_last), 1, quantile, probs = 0.975)
pred_last$sd <- round(apply(exp(sim_last), 1, function(x) sd(x)), 2)
pred_last$cv <- round(apply(exp(sim_last), 1, function(x) sd(x) / mean(x)), 2)

# for testing
sim_year <- subset(pred.adfg.iid.tw.sim, rownames(pred.adfg.iid.tw.sim) %in% c(2024))
pred_year <- filter(pred.adfg.iid.tw$data, Year == 2024)
pred_year$lwr <- apply(exp(sim_year), 1, quantile, probs = 0.025)
pred_year$upr <- apply(exp(sim_year), 1, quantile, probs = 0.975)
pred_year$sd <- round(apply(exp(sim_year), 1, function(x) sd(x)), 2)
pred_year$cv <- round(apply(exp(sim_year), 1, function(x) sd(x) / mean(x)), 2)
pred_year

# function to use
cv_year <- function(year_cv, pred.dat, predsim.dat){
  sim_year <- subset(predsim.dat, rownames(predsim.dat) %in% c(year_cv))
  pred_year <- filter(pred.dat$data, Year == year_cv)
  pred_year$lwr <- apply(exp(sim_year), 1, quantile, probs = 0.025)
  pred_year$upr <- apply(exp(sim_year), 1, quantile, probs = 0.975)
  pred_year$sd <- round(apply(exp(sim_year), 1, function(x) sd(x)), 2)
  pred_year$cv <- round(apply(exp(sim_year), 1, function(x) sd(x) / mean(x)), 2)
  pred_year
}

# years of the survey
adfg.years <- unique(rownames(pred.adfg.iid.tw.sim))
noaa.years <- unique(rownames(pred.noaa.iid.tw.sim))
nbs.years <- unique(rownames(pred.nbs.iid.tw.sim))

# apply cv_year function across all years of the survey for each model
# ADFG models
pred.cv.adfg.tw <- purrr::map_dfr(adfg.years, cv_year, pred.dat = pred.adfg.iid.tw, predsim.dat = pred.adfg.iid.tw.sim)
pred.cv.adfg.tw.s <- purrr::map_dfr(adfg.years, cv_year, pred.dat = pred.adfg.iid.tw.s, predsim.dat = pred.adfg.iid.tw.s.sim)
pred.cv.adfg.tw.d <- purrr::map_dfr(adfg.years, cv_year, pred.dat = pred.adfg.tw.d, predsim.dat = pred.adfg.iid.tw.d.sim)
pred.cv.adfg.tw.s.d <- purrr::map_dfr(adfg.years, cv_year, pred.dat = pred.adfg.tw.s.d, predsim.dat = pred.adfg.iid.tw.s.d.sim)
pred.cv.adfg.dg.d <- purrr::map_dfr(adfg.years, cv_year, pred.dat = pred.adfg.dg.d, predsim.dat = pred.adfg.iid.dg.d.sim)
pred.cv.adfg.dg.s.d <- purrr::map_dfr(adfg.years, cv_year, pred.dat = pred.adfg.dg.s.d, predsim.dat = pred.adfg.iid.dg.s.d.sim)
# NOAA models
pred.cv.noaa.tw <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.tw, predsim.dat = pred.noaa.iid.tw.sim)
pred.cv.noaa.tw.s <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.tw.s, predsim.dat = pred.noaa.iid.tw.s.sim)
pred.cv.noaa.dg <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dg, predsim.dat = pred.noaa.iid.dg.sim)
pred.cv.noaa.dg.s <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dg.s, predsim.dat = pred.noaa.iid.dg.s.sim)
pred.cv.noaa.dl <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dl, predsim.dat = pred.noaa.iid.dl.sim)
pred.cv.noaa.dl.s <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dl.s, predsim.dat = pred.noaa.iid.dl.s.sim)
pred.cv.noaa.tw.d <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.tw.d, predsim.dat = pred.noaa.iid.tw.d.sim)
pred.cv.noaa.tw.s.d <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.tw.s.d, predsim.dat = pred.noaa.iid.tw.s.d.sim)
pred.cv.noaa.dg.d <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dg.d, predsim.dat = pred.noaa.iid.dg.d.sim)
pred.cv.noaa.dg.s.d <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dg.s.d, predsim.dat = pred.noaa.iid.dg.s.d.sim)
pred.cv.noaa.dl.d <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dl.d, predsim.dat = pred.noaa.iid.dl.d.sim)
pred.cv.noaa.dl.s.d <- purrr::map_dfr(noaa.years, cv_year, pred.dat = pred.noaa.iid.dl.s.d, predsim.dat = pred.noaa.iid.dl.s.d.sim)
# NBS models
pred.cv.nbs.tw <- purrr::map_dfr(nbs.years, cv_year, pred.dat = pred.nbs.iid.tw, predsim.dat = pred.nbs.iid.tw.sim)
pred.cv.nbs.tw.s <- purrr::map_dfr(nbs.years, cv_year, pred.dat = pred.nbs.iid.tw.s, predsim.dat = pred.nbs.iid.tw.s.sim)
pred.cv.nbs.tw.d <- purrr::map_dfr(nbs.years, cv_year, pred.dat = pred.nbs.iid.tw.d, predsim.dat = pred.nbs.iid.tw.d.sim)
pred.cv.nbs.tw.s.d <- purrr::map_dfr(nbs.years, cv_year, pred.dat = pred.nbs.iid.tw.s.d, predsim.dat = pred.nbs.iid.tw.s.d.sim)

# plot the CV on the estimates

# for testing
ggplot(pred_last, aes(X, Y, fill = cv)) +
  geom_raster() +
  scale_fill_viridis_c()

cvpred.adfg.tw <- ggplot(pred.cv.adfg.tw, aes(X, Y, fill = cv)) +
  geom_raster() +
  scale_fill_viridis_c() +
  xlab("Eastings (km)") +
  ylab("Northings (km)") +
  scale_x_continuous(breaks=seq(300, 600, by = 300)) +
  facet_wrap(~Year) +
  labs(title = "Coefficient of variation for predicted male abundance from ADF&G trawl survey",
       subtitle = "Tweedie, full prediction area") +
  guides(fill = guide_colourbar(title = "CV"))

# plotting function
cvpred_plot <- function(dat, survey, family, area) {
  ggplot(dat, aes(X, Y, fill = cv)) +
    geom_raster() +
    scale_fill_viridis_c(limits = c(0,10), na.value = "grey") +
    xlab("Eastings (km)") +
    ylab("Northings (km)") +
    scale_x_continuous(breaks=seq(350, 650, by = 300), limits = c(275, 725)) +
    scale_y_continuous(breaks=seq(7000, 7200, by = 100), limits = c(6950, 7275)) +
    facet_wrap(~Year) +
    labs(title = paste0("Coefficient of variation for predicted male abundance"),
         subtitle = paste0(survey, " trawl survey, ", family, ", ", area, " prediction area")) +
    guides(fill = guide_colourbar(title = "CV"))
}

# ADFG plots
pred.cv.plot.adfg.tw <- cvpred_plot(pred.cv.adfg.tw, "ADF&G", "Tweedie, year", "full")
pred.cv.plot.adfg.tw.s <- cvpred_plot(pred.cv.adfg.tw.s, "ADF&G", "Tweedie, year", "reduced")
pred.cv.plot.adfg.tw.d <- cvpred_plot(pred.cv.adfg.tw.d, "ADF&G", "Tweedie, year + depth", "full")
pred.cv.plot.adfg.tw.s.d <- cvpred_plot(pred.cv.adfg.tw.s.d, "ADF&G", "Tweedie, year + depth", "reduced")
pred.cv.plot.adfg.dg.d <- cvpred_plot(pred.cv.adfg.dg.d, "ADF&G", "DG, year + depth", "full")
pred.cv.plot.adfg.dg.s.d <- cvpred_plot(pred.cv.adfg.dg.s.d, "ADF&G", "DG, year + depth", "reduced")
# NOAA plots
pred.cv.plot.noaa.tw <- cvpred_plot(pred.cv.noaa.tw, "NOAA Norton Sound", "Tweedie, year", "full")
pred.cv.plot.noaa.tw.s <- cvpred_plot(pred.cv.noaa.tw.s, "NOAA Norton Sound", "Tweedie, year", "reduced")
pred.cv.plot.noaa.dg <- cvpred_plot(pred.cv.noaa.dg, "NOAA Norton Sound", "DG, year", "full")
pred.cv.plot.noaa.dg.s <- cvpred_plot(pred.cv.noaa.dg.s, "NOAA Norton Sound", "DG, year", "reduced")
pred.cv.plot.noaa.dl <- cvpred_plot(pred.cv.noaa.dl, "NOAA Norton Sound", "DL, year", "full")
pred.cv.plot.noaa.dl.s <- cvpred_plot(pred.cv.noaa.dl.s, "NOAA Norton Sound", "DL, year", "reduced")
pred.cv.plot.noaa.tw.d <- cvpred_plot(pred.cv.noaa.tw.d, "NOAA Norton Sound", "Tweedie, year + depth", "full")
pred.cv.plot.noaa.tw.s.d <- cvpred_plot(pred.cv.noaa.tw.s.d, "NOAA Norton Sound", "Tweedie, year + depth", "reduced")
pred.cv.plot.noaa.dg.d <- cvpred_plot(pred.cv.noaa.dg.d, "NOAA Norton Sound", "DG, year + depth", "full")
pred.cv.plot.noaa.dg.s.d <- cvpred_plot(pred.cv.noaa.dg.s.d, "NOAA Norton Sound", "DG, year + depth", "reduced")
pred.cv.plot.noaa.dl.d <- cvpred_plot(pred.cv.noaa.dl.d, "NOAA Norton Sound", "DL, year + depth", "full")
pred.cv.plot.noaa.dl.s.d <- cvpred_plot(pred.cv.noaa.dl.s.d, "NOAA Norton Sound", "DL, year + depth", "reduced")
# NBS plots
pred.cv.plot.nbs.tw <- cvpred_plot(pred.cv.nbs.tw, "NOAA NBS", "Tweedie, year", "full")
pred.cv.plot.nbs.tw.s <- cvpred_plot(pred.cv.nbs.tw.s, "NOAA NBS", "Tweedie, year", "reduced")
pred.cv.plot.nbs.tw.d <- cvpred_plot(pred.cv.nbs.tw.d, "NOAA NBS", "Tweedie, year + depth", "full")
pred.cv.plot.nbs.tw.s.d <- cvpred_plot(pred.cv.nbs.tw.s.d, "NOAA NBS", "Tweedie, year + depth", "reduced")

# save plots
# ADFG plots
ggsave(file.path(plotdir.ns, "predcv_adfg_tw.png"), plot = pred.cv.plot.adfg.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_adfg_tw_s.png"), plot = pred.cv.plot.adfg.tw.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_adfg_tw_d.png"), plot = pred.cv.plot.adfg.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_adfg_tw_s_d.png"), plot = pred.cv.plot.adfg.tw.s.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_adfg_dg_d.png"), plot = pred.cv.plot.adfg.dg.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_adfg_dg_s_d.png"), plot = pred.cv.plot.adfg.dg.s.d, height = 5, width = 7, units = "in")
# NOAA plots
ggsave(file.path(plotdir.ns, "predcv_noaa_tw.png"), plot = pred.cv.plot.noaa.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_tw_s.png"), plot = pred.cv.plot.noaa.tw.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dg.png"), plot = pred.cv.plot.noaa.dg, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dg_s.png"), plot = pred.cv.plot.noaa.dg.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dl.png"), plot = pred.cv.plot.noaa.dl, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dl_s.png"), plot = pred.cv.plot.noaa.dl.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_tw_d.png"), plot = pred.cv.plot.noaa.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_tw_s_d.png"), plot = pred.cv.plot.noaa.tw.s.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dg_d.png"), plot = pred.cv.plot.noaa.dg.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dg_s_d.png"), plot = pred.cv.plot.noaa.dg.s.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dl_d.png"), plot = pred.cv.plot.noaa.dl.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_noaa_dl_s_d.png"), plot = pred.cv.plot.noaa.dl.s.d, height = 5, width = 7, units = "in")
# NBS plots
ggsave(file.path(plotdir.ns, "predcv_nbs_tw.png"), plot = pred.cv.plot.nbs.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_nbs_tw_s.png"), plot = pred.cv.plot.nbs.tw.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_nbs_tw_d.png"), plot = pred.cv.plot.nbs.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "predcv_nbs_tw_s_d.png"), plot = pred.cv.plot.nbs.tw.s.d, height = 5, width = 7, units = "in")

# spatial biomass prediction plots ----

# plotting function
plot_map <- function(dat, column, estimate, survey, dist, area) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c(
      trans = "sqrt",
      # trim extreme high values to make spatial variation more visible
      na.value = "yellow", limits = c(0, quantile(exp(estimate), 0.995))
    ) +
    xlab("Eastings (km)") +
    ylab("Northings (km)") +
    #scale_x_continuous(breaks=seq(round(min(dat$X),0), round(max(dat$X),0), by = 200)) +
    scale_x_continuous(breaks=seq(350, 650, by = 300), limits = c(275, 725)) +
    scale_y_continuous(breaks=seq(7000, 7200, by = 100), limits = c(6950, 7275)) +
    facet_wrap(~Year) +
    labs(title = paste0("Predicted male abundance from ", survey, " trawl survey"),
         subtitle = paste0(dist, ", ", area, " prediction grid")) +
    guides(fill = guide_colourbar(title = expression(atop(textstyle("Crab per km"^2)))))
}

# create plots
# ADFG plots
sppred.adfg.iid.tw <- plot_map(dat=pred.adfg.iid.tw$data, column=exp(est), estimate = pred.adfg.iid.tw$data$est, survey="ADF&G", dist="Tweedie, year", area = "full") 
sppred.adfg.iid.tw.s <- plot_map(dat=pred.adfg.iid.tw.s$data, column=exp(est), estimate = pred.adfg.iid.tw.s$data$est, survey="ADF&G", dist="Tweedie, year", area = "reduced") 
sppred.adfg.iid.tw.d <- plot_map(dat=pred.adfg.tw.d$data, column=exp(est), estimate = pred.adfg.tw.d$data$est, survey="ADF&G", dist="Tweedie, year + depth", area = "full") 
sppred.adfg.iid.tw.s.d <- plot_map(dat=pred.adfg.tw.s.d$data, column=exp(est), estimate = pred.adfg.tw.s.d$data$est, survey="ADF&G", dist="Tweedie, year + depth", area = "reduced") 
sppred.adfg.iid.dg.d <- plot_map(dat=pred.adfg.dg.d$data, column=exp(est1), estimate = pred.adfg.dg.d$data$est1, survey="ADF&G", dist="DG, year + depth", area = "full") 
sppred.adfg.iid.dg.s.d <- plot_map(dat=pred.adfg.dg.s.d$data, column=exp(est1), estimate = pred.adfg.dg.s.d$data$est1, survey="ADF&G", dist="DG, year + depth", area = "reduced") 
# NOAA plots
sppred.noaa.iid.tw <- plot_map(dat=pred.noaa.iid.tw$data, column=exp(est), estimate = pred.noaa.iid.tw$data$est, survey="NOAA Norton Sound", dist="Tweedie, year", area = "full") 
sppred.noaa.iid.tw.s <- plot_map(dat=pred.noaa.iid.tw.s$data, column=exp(est), estimate = pred.noaa.iid.tw.s$data$est, survey="NOAA Norton Sound", dist="Tweedie, year", area = "reduced") 
sppred.noaa.iid.dg <- plot_map(dat=pred.noaa.iid.dg$data, column=exp(est1), estimate = pred.noaa.iid.dg$data$est1, survey="NOAA Norton Sound", dist="DG, year", area = "full") 
sppred.noaa.iid.dg.s <- plot_map(dat=pred.noaa.iid.dg.s$data, column=exp(est1), estimate = pred.noaa.iid.dg.s$data$est1, survey="NOAA Norton Sound", dist="DG, year", area = "reduced")
sppred.noaa.iid.dl <- plot_map(dat=pred.noaa.iid.dl$data, column=exp(est1), estimate = pred.noaa.iid.dl$data$est1, survey="NOAA Norton Sound", dist="DL, year", area = "full") 
sppred.noaa.iid.dl.s <- plot_map(dat=pred.noaa.iid.dl.s$data, column=exp(est1), estimate = pred.noaa.iid.dl.s$data$est1, survey="NOAA Norton Sound", dist="DL, year", area = "reduced")
sppred.noaa.iid.tw.d <- plot_map(dat=pred.noaa.iid.tw.d$data, column=exp(est), estimate = pred.noaa.iid.tw.d$data$est, survey="NOAA Norton Sound", dist="Tweedie, year + depth", area = "full") 
sppred.noaa.iid.tw.s.d <- plot_map(dat=pred.noaa.iid.tw.s.d$data, column=exp(est), estimate = pred.noaa.iid.tw.s.d$data$est, survey="NOAA Norton Sound", dist="Tweedie, year + depth", area = "reduced") 
sppred.noaa.iid.dg.d <- plot_map(dat=pred.noaa.iid.dg.d$data, column=exp(est1), estimate = pred.noaa.iid.dg.d$data$est1, survey="NOAA Norton Sound", dist="DG, year + depth", area = "full") 
sppred.noaa.iid.dg.s.d <- plot_map(dat=pred.noaa.iid.dg.s.d$data, column=exp(est1), estimate = pred.noaa.iid.dg.s.d$data$est1, survey="NOAA Norton Sound", dist="DG, year + depth", area = "reduced") 
sppred.noaa.iid.dl.d <- plot_map(dat=pred.noaa.iid.dl.d$data, column=exp(est1), estimate = pred.noaa.iid.dl.d$data$est1, survey="NOAA Norton Sound", dist="DL, year + depth", area = "full") 
sppred.noaa.iid.dl.s.d <- plot_map(dat=pred.noaa.iid.dl.s.d$data, column=exp(est1), estimate = pred.noaa.iid.dl.s.d$data$est1, survey="NOAA Norton Sound", dist="DL, year + depth", area = "reduced") 
# NBS plots
sppred.nbs.iid.tw <- plot_map(dat=pred.nbs.iid.tw$data, column=exp(est), estimate = pred.nbs.iid.tw$data$est, survey="NOAA NBS", dist="Tweedie, year", area = "full")
sppred.nbs.iid.tw.s <- plot_map(dat=pred.nbs.iid.tw.s$data, column=exp(est), estimate = pred.nbs.iid.tw.s$data$est, survey="NOAA NBS", dist="Tweedie, year", area = "reduced")
sppred.nbs.iid.tw.d <- plot_map(dat=pred.nbs.iid.tw.d$data, column=exp(est), estimate = pred.nbs.iid.tw.d$data$est, survey="NOAA NBS", dist="Tweedie, year + depth", area = "full") 
sppred.nbs.iid.tw.s.d <- plot_map(dat=pred.nbs.iid.tw.s.d$data, column=exp(est), estimate = pred.nbs.iid.tw.s.d$data$est, survey="NOAA NBS", dist="Tweedie, year + depth", area = "reduced") 

# save plots
# ADFG plots
ggsave(file.path(plotdir.ns, "sppred_adfg_iid_tw.png"), plot = sppred.adfg.iid.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_adfg_iid_tw_s.png"), plot = sppred.adfg.iid.tw.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_adfg_iid_tw_d.png"), plot = sppred.adfg.iid.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_adfg_iid_tw_s_d.png"), plot = sppred.adfg.iid.tw.s.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_adfg_iid_dg_d.png"), plot = sppred.adfg.iid.dg.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_adfg_iid_dg_s_d.png"), plot = sppred.adfg.iid.dg.s.d, height = 5, width = 7, units = "in")
# NOAA plots
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_tw.png"), plot = sppred.noaa.iid.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_tw_s.png"), plot = sppred.noaa.iid.tw.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dg.png"), plot = sppred.noaa.iid.dg, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dg_s.png"), plot = sppred.noaa.iid.dg.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dl.png"), plot = sppred.noaa.iid.dl, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dl_s.png"), plot = sppred.noaa.iid.dl.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_tw_d.png"), plot = sppred.noaa.iid.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_tw_s_d.png"), plot = sppred.noaa.iid.tw.s.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dg_d.png"), plot = sppred.noaa.iid.dg.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dg_s_d.png"), plot = sppred.noaa.iid.dg.s.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dl_d.png"), plot = sppred.noaa.iid.dl.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_noaa_iid_dl_s_d.png"), plot = sppred.noaa.iid.dl.s.d, height = 5, width = 7, units = "in")
# NBS plots
ggsave(file.path(plotdir.ns, "sppred_nbs_iid_tw.png"), plot = sppred.nbs.iid.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_nbs_iid_tw_s.png"), plot = sppred.nbs.iid.tw.s, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_nbs_iid_tw_d.png"), plot = sppred.nbs.iid.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "sppred_nbs_iid_tw_s_d.png"), plot = sppred.nbs.iid.tw.s.d, height = 5, width = 7, units = "in")

# *************************************************************************************************
# generate index ----
# *************************************************************************************************

# ADFG indices
index.adfg.iid.tw <- get_index(pred.adfg.iid.tw, area = pred.adfg.iid.tw$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.adfg.iid.tw.s <- get_index(pred.adfg.iid.tw.s, area = pred.adfg.iid.tw.s$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.adfg.iid.tw.d <- get_index(pred.adfg.tw.d, area = pred.adfg.tw.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.adfg.iid.tw.s.d <- get_index(pred.adfg.tw.s.d, area = pred.adfg.tw.s.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.adfg.iid.dg.d <- get_index(pred.adfg.dg.d, area = pred.adfg.dg.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.adfg.iid.dg.s.d <- get_index(pred.adfg.dg.s.d, area = pred.adfg.dg.s.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))

# NOAA indices
index.noaa.iid.tw <- get_index(pred.noaa.iid.tw, area = pred.noaa.iid.tw$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.tw.s <- get_index(pred.noaa.iid.tw.s, area = pred.noaa.iid.tw.s$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dg <- get_index(pred.noaa.iid.dg, area = pred.noaa.iid.dg$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dg.s <- get_index(pred.noaa.iid.dg.s, area = pred.noaa.iid.dg.s$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dl <- get_index(pred.noaa.iid.dl, area = pred.noaa.iid.dl$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dl.s <- get_index(pred.noaa.iid.dl.s, area = pred.noaa.iid.dl.s$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.tw.d <- get_index(pred.noaa.iid.tw.d, area = pred.noaa.iid.tw.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.tw.s.d <- get_index(pred.noaa.iid.tw.s.d, area = pred.noaa.iid.tw.s.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dg.d <- get_index(pred.noaa.iid.dg.d, area = pred.noaa.iid.dg.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dg.s.d <- get_index(pred.noaa.iid.dg.s.d, area = pred.noaa.iid.dg.s.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dl.d <- get_index(pred.noaa.iid.dl.d, area = pred.noaa.iid.dl.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.noaa.iid.dl.s.d <- get_index(pred.noaa.iid.dl.s.d, area = pred.noaa.iid.dl.s.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))

# NBS indices
index.nbs.iid.tw <- get_index(pred.nbs.iid.tw, area = pred.nbs.iid.tw$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.nbs.iid.tw.s <- get_index(pred.nbs.iid.tw.s, area = pred.nbs.iid.tw.s$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.nbs.iid.tw.d <- get_index(pred.nbs.iid.tw.d, area = pred.nbs.iid.tw.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))
index.nbs.iid.tw.s.d <- get_index(pred.nbs.iid.tw.s.d, area = pred.nbs.iid.tw.s.d$data$Area.km2, bias_correct = TRUE) %>% mutate(cv = sqrt(exp(se^2) - 1))

# save index values
# ADFG indices
write.csv(index.adfg.iid.tw, paste0(here::here(), "/NSRKC/output/adfg_index_iid_tw.csv"))
write.csv(index.adfg.iid.tw.s, paste0(here::here(), "/NSRKC/output/adfg_index_iid_tw_s.csv"))
write.csv(index.adfg.iid.tw.d, paste0(here::here(), "/NSRKC/output/adfg_index_iid_tw_d.csv"))
write.csv(index.adfg.iid.tw.s.d, paste0(here::here(), "/NSRKC/output/adfg_index_iid_tw_s_d.csv"))
write.csv(index.adfg.iid.dg.d, paste0(here::here(), "/NSRKC/output/adfg_index_iid_dg_d.csv"))
write.csv(index.adfg.iid.dg.s.d, paste0(here::here(), "/NSRKC/output/adfg_index_iid_dg_s_d.csv"))
# NOAA indices
write.csv(index.noaa.iid.tw, paste0(here::here(), "/NSRKC/output/noaa_index_iid_tw.csv"))
write.csv(index.noaa.iid.tw.s, paste0(here::here(), "/NSRKC/output/noaa_index_iid_tw_s.csv"))
write.csv(index.noaa.iid.dg, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dg.csv"))
write.csv(index.noaa.iid.dg.s, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dg_s.csv"))
write.csv(index.noaa.iid.dl, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dl.csv"))
write.csv(index.noaa.iid.dl.s, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dl_s.csv"))
write.csv(index.noaa.iid.tw.d, paste0(here::here(), "/NSRKC/output/noaa_index_iid_tw_d.csv"))
write.csv(index.noaa.iid.tw.s.d, paste0(here::here(), "/NSRKC/output/noaa_index_iid_tw_s_d.csv"))
write.csv(index.noaa.iid.dg.d, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dg_d.csv"))
write.csv(index.noaa.iid.dg.s.d, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dg_s_d.csv"))
write.csv(index.noaa.iid.dl.d, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dl_d.csv"))
write.csv(index.noaa.iid.dl.s.d, paste0(here::here(), "/NSRKC/output/noaa_index_iid_dl_s_d.csv"))
# NBS indices
write.csv(index.nbs.iid.tw, paste0(here::here(), "/NSRKC/output/nbs_index_iid_tw.csv"))
write.csv(index.nbs.iid.tw.s, paste0(here::here(), "/NSRKC/output/nbs_index_iid_tw_s.csv"))
write.csv(index.nbs.iid.tw.d, paste0(here::here(), "/NSRKC/output/nbs_index_iid_tw_d.csv"))
write.csv(index.nbs.iid.tw.s.d, paste0(here::here(), "/NSRKC/output/nbs_index_iid_tw_s_d.csv"))

# plot index
ggplot(index.adfg.iid.tw, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.adfg.iid.tw.s, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.adfg.iid.tw.sims, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.adfg.tw.d, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

ggplot(index.adfg.dg.d, aes(Year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)

# create data frames containing multiple indices
index.adfg.iid.tw.all <- rbind(
  index.adfg.iid.tw %>% mutate(model = "TW, year, full"),
  index.adfg.iid.tw.s %>% mutate(model = "TW, year, reduced"),
  index.adfg.tw.d %>% mutate(model = "TW, year + depth, full"),
  index.adfg.tw.s.d %>% mutate(model = "TW, year + depth, reduced"),
  index.adfg.dg.d %>% mutate(model = "DG, year + depth, full"),
  index.adfg.dg.s.d %>% mutate(model = "DG, year + depth, reduced")
)

index.noaa.iid.all <- rbind(
  index.noaa.iid.tw %>% mutate(model = "TW, year, full"),
  index.noaa.iid.tw.s %>% mutate(model = "TW, year, reduced"),
  index.noaa.iid.dg %>% mutate(model = "DG, year, full"),
  index.noaa.iid.dg.s %>% mutate(model = "DG, year, reduced"),
  index.noaa.iid.dl %>% mutate(model = "DL, year, full"),
  index.noaa.iid.dl.s %>% mutate(model = "DL, year, reduced"),
  index.noaa.iid.tw.d %>% mutate(model = "TW, year + depth, full"),
  index.noaa.iid.tw.s.d %>% mutate(model = "TW, year + depth, reduced"),
  index.noaa.iid.dg.d %>% mutate(model = "DG, year + depth, full"),
  index.noaa.iid.dg.s.d %>% mutate(model = "DG, year + depth, reduced"),
  index.noaa.iid.dl.d %>% mutate(model = "DL, year + depth, full"),
  index.noaa.iid.dl.s.d %>% mutate(model = "DL, year + depth, reduced")
)

index.nbs.iid.all <- rbind(
  index.nbs.iid.tw %>% mutate(model = "year, full"),
  index.nbs.iid.tw.s %>% mutate(model = "year, reduced"),
  index.nbs.iid.tw.d %>% mutate(model = "year + depth, full"),
  index.nbs.iid.tw.s.d %>% mutate(model = "year + depth, reduced")
)

# ************************************************************************************************
# compare predicted index to observations ----
# ************************************************************************************************

obs.pred.adfg <- abund.adfg.mod %>%
  mutate(obs = 1000 * obs_index) %>% 
  rename(Year = year) %>% 
  left_join(index.adfg.iid.tw.all) 

obs.pred.noaa <- abund.noaa.mod %>%
  mutate(obs = 1000 * obs_index) %>% 
  rename(Year = year) %>% 
  left_join(index.noaa.iid.all)

obs.pred.nbs <- abund.nbs.mod %>%
  mutate(obs = 1000 * obs_index) %>% 
  rename(Year = year) %>% 
  left_join(index.nbs.iid.all)

# index plotting function
index_plot <- function(dat, filt, survey){
    dat %>%
    filter(model %in% filt) %>%
    ggplot()+
    geom_line(aes(x = Year, y = est/1000, color = model))+
    geom_ribbon(aes(x = Year, y = est/1000, ymin = lwr/1000, ymax = upr/1000, fill = model), alpha = 0.2) +
    geom_point(aes(x = Year, y = obs_index, shape = "estimate"), color = "grey20")+
    geom_errorbar(aes(x = Year, ymin = obs_l95, ymax = obs_u95), width = 0, color = "grey20")+
    xlab('Year') + 
    ylab('Red king crab abundance estimate (1000s of crab)') +
    scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, 26000)) +
    coord_cartesian(ylim = c(0, NA))+
    scale_color_manual(values = cbpalette) +
    scale_fill_manual(values = cbpalette) +
    labs(fill = "Model-based", color = "Model-based", shape = "Design-based") +
    theme(legend.position = c(0.84,0.75)) +
    ggtitle(paste0(survey, " trawl survey"))
}

# ADFG index plots ----
adfg.fit.iid.tw <- index_plot(obs.pred.adfg, c("TW, year, full", "TW, year, reduced"), "ADF&G")
adfg.fit.iid.tw.d <- index_plot(obs.pred.adfg, c("TW, year + depth, full", "TW, year + depth, reduced"), "ADF&G")
adfg.fit.iid.dg.d <- index_plot(obs.pred.adfg, c("DG, year + depth, full", "DG, year + depth, reduced"), "ADF&G")

ggsave(file.path(plotdir.ns, "nsrkc_adfg_iid_tw_fit.png"), plot = adfg.fit.iid.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_adfg_iid_tw_d_fit.png"), plot = adfg.fit.iid.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_adfg_iid_dg_d_fit.png"), plot = adfg.fit.iid.dg.d, height = 5, width = 7, units = "in")

# NOAA index plots ----
noaa.fit.iid.tw <- index_plot(obs.pred.noaa, c("TW, year, full", "TW, year, reduced"), "NOAA Norton Sound")
noaa.fit.iid.dg <- index_plot(obs.pred.noaa, c("DG, year, full", "DG, year, reduced"), "NOAA Norton Sound")
noaa.fit.iid.dl <- index_plot(obs.pred.noaa, c("DL, year, full", "DL, year, reduced"), "NOAA Norton Sound")
noaa.fit.iid.tw.d <- index_plot(obs.pred.noaa, c("TW, year + depth, full", "TW, year + depth, reduced"), "NOAA Norton Sound")
noaa.fit.iid.dg.d <- index_plot(obs.pred.noaa, c("DG, year + depth, full", "DG, year + depth, reduced"), "NOAA Norton Sound")
noaa.fit.iid.dl.d <- index_plot(obs.pred.noaa, c("DL, year + depth, full", "DL, year + depth, reduced"), "NOAA Norton Sound")
  
ggsave(file.path(plotdir.ns, "nsrkc_noaa_iid_tw_fit.png"), plot = noaa.fit.iid.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_noaa_iid_dg_fit.png"), plot = noaa.fit.iid.dg, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_noaa_iid_dl_fit.png"), plot = noaa.fit.iid.dl, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_noaa_iid_tw_d_fit.png"), plot = noaa.fit.iid.tw.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_noaa_iid_dg_d_fit.png"), plot = noaa.fit.iid.dg.d, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_noaa_iid_dl_d_fit.png"), plot = noaa.fit.iid.dl.d, height = 5, width = 7, units = "in")

# NBS index plots ----
nbs.fit.iid.tw <- index_plot(obs.pred.nbs, c("year, full", "year, reduced"), "NOAA Northern Bering Sea")
nbs.fit.iid.tw.d <- index_plot(obs.pred.nbs, c("year + depth, full", "year + depth, reduced"), "NOAA Northern Bering Sea")

ggsave(file.path(plotdir.ns, "nsrkc_nbs_iid_tw_fit.png"), plot = nbs.fit.iid.tw, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "nsrkc_nbs_iid_tw_d_fit.png"), plot = nbs.fit.iid.tw.d, height = 5, width = 7, units = "in")

# ************************************************************************************************
# plot CVs for design-based versus model-based indices ----
# ************************************************************************************************

# CV plotting function
cv_index_plot <- function(dat, survey) {
  ggplot(dat, aes(x = Year, y = cv, group = type, color = type)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = cbpalette) +
    guides(color = guide_legend(title = "Estimate")) +
    ylab("Coefficient of variation") +
    theme(legend.position = c(0.8,0.8)) +
    ggtitle(paste0(survey, " trawl survey"))
}

# ADFG
comp.cv.adfg <- abund.adfg.mod %>%
  dplyr::select(c(year, obs_cv)) %>%
  mutate(type = "design-based") %>%
  rename("Year" = year, "cv" = obs_cv) %>%
  rbind(index.adfg.iid.tw %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year, full"))  %>%
  rbind(index.adfg.iid.tw.s %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year, reduced")) %>%
  rbind(index.adfg.iid.tw.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year + depth, full"))  %>%
  rbind(index.adfg.iid.tw.s.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year + depth, reduced")) %>%
  rbind(index.adfg.iid.dg.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DG, year + depth, full"))  %>%
  rbind(index.adfg.iid.dg.s.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DG, year + depth, reduced"))

cv.ind.adfg <- cv_index_plot(comp.cv.adfg, "ADF&G")

ggsave(file.path(plotdir.ns, "cv_ind_adfg.png"), plot = cv.ind.adfg, height = 5, width = 7, units = "in")

# NOAA
comp.cv.noaa <- abund.noaa.mod %>%
  dplyr::select(c(year, obs_cv)) %>%
  mutate(type = "design-based") %>%
  rename("Year" = year, "cv" = obs_cv) %>%
  rbind(index.noaa.iid.tw %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year, full"))  %>%
  rbind(index.noaa.iid.tw.s %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year, reduced")) %>%
  rbind(index.noaa.iid.dg %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DG, year, full"))  %>%
  rbind(index.noaa.iid.dg.s %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DG, year, reduced")) %>%
  rbind(index.noaa.iid.dl %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DL, year, full"))  %>%
  rbind(index.noaa.iid.dl.s %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DL, year, reduced")) %>%
  rbind(index.noaa.iid.tw.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year + depth, full"))  %>%
  rbind(index.noaa.iid.tw.s.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "Tweedie, year + depth, reduced")) %>%
  rbind(index.noaa.iid.dg.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DG, year + depth, full"))  %>%
  rbind(index.noaa.iid.dg.s.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DG, year + depth, reduced")) %>%
  rbind(index.noaa.iid.dl.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DL, year + depth, full"))  %>%
  rbind(index.noaa.iid.dl.s.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "DL, year + depth, reduced"))

cv.ind.noaa <- cv_index_plot(comp.cv.noaa %>% filter(str_detect(type, "depth", negate = TRUE)), "NOAA Norton Sound")
cv.ind.noaa.d <- cv_index_plot(comp.cv.noaa %>% filter(str_detect(type, "depth") | type == "design-based"), "NOAA Norton Sound")

ggsave(file.path(plotdir.ns, "cv_ind_noaa.png"), plot = cv.ind.noaa, height = 5, width = 7, units = "in")
ggsave(file.path(plotdir.ns, "cv_ind_noaa_d.png"), plot = cv.ind.noaa.d, height = 5, width = 7, units = "in")

# NBS
comp.cv.nbs <- abund.nbs.mod %>%
  dplyr::select(c(year, obs_cv)) %>%
  mutate(type = "design-based") %>%
  rename("Year" = year, "cv" = obs_cv) %>%
  rbind(index.nbs.iid.tw %>% dplyr::select(c(Year, cv)) %>% mutate(type = "year, full"))  %>%
  rbind(index.nbs.iid.tw.s %>% dplyr::select(c(Year, cv)) %>% mutate(type = "year, reduced")) %>%
  rbind(index.nbs.iid.tw.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "year + depth, full"))  %>%
  rbind(index.nbs.iid.tw.s.d %>% dplyr::select(c(Year, cv)) %>% mutate(type = "year + depth, reduced"))

cv.ind.nbs <- cv_index_plot(comp.cv.nbs, "NOAA Northern Bering Sea")

ggsave(file.path(plotdir.ns, "cv_ind_nbs_tw.png"), plot = cv.ind.nbs, height = 5, width = 7, units = "in")

# ************************************************************************************************
# calculate Ratio, TimeOut, and Magnitude ----
# these are metrics for comparing model- and design-based indices
# https://doi.org/10.1016/j.fishres.2024.107009
# ************************************************************************************************

# ratio: take mean of (model-based abundance index in a given year/design-based index in a given year)

obs.pred.all <- rbind(obs.pred.adfg, obs.pred.noaa, obs.pred.nbs)

compare.indices <- obs.pred.all %>%
  group_by(fleet, model) %>%
  arrange_at(.vars = vars(fleet, model)) %>%
  arrange(match(fleet, c("ADFG_Trawl", "NMFS_Trawl", "NBS_Trawl"))) %>%
  mutate(est_index = est/1000) %>%
  # calculate Ratio
  #mutate(ann.ratio = est_index/obs_index) %>%
  #mutate(Ratio = mean(ann.ratio)) %>%
  mutate(mean.est.index = mean(est_index), mean.obs.index = mean(obs_index)) %>%
  mutate(Ratio = mean.est.index/mean.obs.index) %>%
  # calculate CV ratio
  #mutate(ann.cv.ratio = cv/obs_cv) %>%
  #mutate(CV.Ratio = mean(ann.cv.ratio)) %>%
  mutate(mean.est.cv = mean(cv), mean.obs.cv = mean(obs_cv)) %>%
  mutate(CV.Ratio = mean.est.cv/mean.obs.cv) %>%
  # calculate TimeOut
  mutate(out = case_when(
    est_index > obs_u95 | est_index < obs_l95 ~ 1,
    .default = 0
  )) %>%
  mutate(TimeOut = sum(out)/n()) %>%
  # calculate Magnitude
  mutate(mag.ann = case_when(
    out == 1 ~ (est_index - obs_index)/obs_index,
    .default = 0
  )) %>%
  mutate(Magnitude = sum(mag.ann)) %>%
  # calculate sum of squared deviations
  mutate(ssd.est = sum((est_index - mean(est_index))^2)) %>%
  mutate(ssd.obs = sum((obs_index - mean(obs_index))^2)) %>%
  mutate(ssd.ratio = ssd.est/ssd.obs) %>%
  # calculate sample variance
  mutate(var.est = var(est_index)) %>%
  mutate(var.obs = var(obs_index)) %>%
  mutate(var.ratio = var.est/var.obs) %>%
  # select fields needed
  dplyr::select(fleet, model, Ratio, TimeOut, Magnitude, CV.Ratio, var.ratio) %>%
  unique() %>%
  mutate(grid = case_when(
    grepl("full", model) == TRUE ~ "full",
    grepl("reduced", model) == TRUE ~ "reduced"
  )) %>%
  mutate(effects = case_when(
    str_detect(model, "depth") == TRUE ~ "year + depth",
    str_detect(model, "depth") == FALSE ~ "year"
  )) %>%
  mutate(model = case_when(
    str_detect(model, "TW") == TRUE ~ "Tweedie",
    str_detect(model, "DG") == TRUE ~ "delta-gamma",
    str_detect(model, "DL") == TRUE ~ "delta-lognormal",
    fleet == "NBS_Trawl" ~ "Tweedie"
  )) %>%
  arrange_at(.vars = vars(fleet, model, effects)) %>%
  arrange(match(fleet, c("ADFG_Trawl", "NMFS_Trawl", "NBS_Trawl")))

write.csv(compare.indices, paste0(here::here(), "/NSRKC/output/ns_compare_indices.csv"))

# ************************************************************************************************
# plot VAST index ----
# ************************************************************************************************

vast.index.adfg.iid.dg.30kn <- read.csv(paste0(here::here(), "/NSRKC/VAST output/NSRKC/NSRKC_MALE_GE64_30kts_V2025_Epsilon2_Omega1_Omega2_disabled/Index.csv")) %>%
  mutate(lwr = Estimate - 1.96*Std..Error.for.Estimate, upr = Estimate + 1.96*Std..Error.for.Estimate) %>%
  rename(Year = Time) %>%
  mutate(vast.est = Estimate * 1000, vast.lwr = lwr * 1000, vast.upr = upr * 1000)
vast.index.adfg.iid.dg.50kn <- read.csv(paste0(here::here(), "/NSRKC/VAST output/NSRKC/NSRKC_MALE_GE64_50kts_V2025_Epsilon2_Omega1_Omega2_disabled/Index.csv")) %>%
  mutate(lwr = Estimate - 1.96*Std..Error.for.Estimate, upr = Estimate + 1.96*Std..Error.for.Estimate) %>%
  rename(Year = Time) %>%
  mutate(vast.est = Estimate * 1000, vast.lwr = lwr * 1000, vast.upr = upr * 1000)
vast.index.adfg.iid.dg.100kn <- read.csv(paste0(here::here(), "/NSRKC/VAST output/NSRKC/NSRKC_MALE_GE64_100kts_V2025_Epsilon2_Omega1_Omega2_disabled/Index.csv")) %>%
  mutate(lwr = Estimate - 1.96*Std..Error.for.Estimate, upr = Estimate + 1.96*Std..Error.for.Estimate) %>%
  rename(Year = Time) %>%
  mutate(vast.est = Estimate * 1000, vast.lwr = lwr * 1000, vast.upr = upr * 1000)


ggplot(vast.index.adfg.iid.dg.30kn, aes(Year, vast.est)) + geom_line() +
  geom_ribbon(aes(ymin = vast.lwr, ymax = vast.upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (numbers)') +
  scale_y_continuous(label=scales::comma)


vast.index.adfg.iid.dg <- rbind(
  vast.index.adfg.iid.dg.30kn %>% mutate(model = "IID, delta gamma, 30 kn"),
  vast.index.adfg.iid.dg.50kn %>% mutate(model = "IID, delta gamma, 50 kn"),
  vast.index.adfg.iid.dg.100kn %>% mutate(model = "IID, delta gamma, 100 kn")
)

vast.obs.pred.adfg <- abund.crab.adfg %>%
  left_join(vast.index.adfg.iid.dg) 


vast.adfg.fit.iid.dg <- ggplot(vast.obs.pred.adfg)+
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

ggsave(file.path(plotdir.ns, "vast_nsrkc_iid_dg_fit.png"), plot = vast.adfg.fit.iid.dg, height = 5, width = 7, units = "in")

