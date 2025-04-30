### PURPOSE ----------------------------------------------------------------------
# To generate and compare plots of model-based indices of biomass for EBS snow crab for males >=95mm
# and mature females. Time range is 1980-present.

# Author: Emily Ryznar

# TO DOs:
# 1)

### LOAD LIBRARIES/FUNCTIONS/DATA --------------------------------------------------------
source("./SNOW/Scripts/snow_load_libs_functions.R")


### Load survey observations ----------
m.surv
mf.surv

### EBS ----
years <- c(1980:2019, 2021:2024)

# Filter pred_grid by stock, transform to UTM, replicate by number of years
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", years) %>%
  rename(lon = X, lat = Y)

## Load and process index data ----
## Males sdmTMB 
data <- snow.male95.cpue
category <- "Male95"
region <- "EBS"

# 50 - DG
m.index.bio50 <- read.csv(paste0(dir, "Output/Indices/snow_EBS_Male95_50_DG_biomassindex.csv")) %>%
                  rename(biomass = est, Year = year) %>%
                  mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                  mutate(knots = 50, category = category, region = region, dist = "DG")


# 120 - DG
m.index.bio120 <- read.csv(paste0(dir, "Output/Indices/snow_EBS_Male95_120_DG_biomassindex.csv")) %>%
                    rename(biomass = est, Year = year) %>%
                    mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
                    mutate(knots = 120, category = category, region = region, dist = "DG")

## Mature females sdmTMB  
data <- snow.matfem.cpue
category <- "Mature female"
region <- "EBS"

# 50 - DG
mf.index.bio50 <- read.csv(paste0(dir, "Output/Indices/snow_EBS_Mature female_50_DG_biomassindex.csv")) %>%
  rename(biomass = est, Year = year) %>%
  mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
  mutate(knots = 50, category = category, region = region, dist = "DG")


# 120 - DG
mf.index.bio120 <- read.csv(paste0(dir, "Output/Indices/snow_EBS_Mature female_120_DG_biomassindex.csv")) %>%
  rename(biomass = est, Year = year) %>%
  mutate(biomass = biomass/1000, lwr = lwr/1000, upr= upr/1000)%>% 
  mutate(knots = 120, category = category, region = region, dist = "DG")

## Bind all indices

ind.dat <- rbind(m.index.bio50, m.index.bio120, mf.index.bio50, mf.index.bio120)

## Bind survey dat
surv.dat <- rbind(m.surv, mf.surv) %>%
            filter(REGION == "EBS")


## Plot indices ----
ggplot()+
  geom_ribbon(ind.dat, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(knots)), alpha = 0.4) +
  geom_line(ind.dat, mapping = aes(Year, biomass, color = as.factor(knots)))+
  geom_point(surv.dat,
             mapping = aes(YEAR, BIOMASS_MT), color = "grey20", size = 0.75)+
  geom_errorbar(surv.dat,
                mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI), color = "grey20", width = 0)+
  facet_wrap(~factor(category, levels = c("Male95", "Mature female")), scales = "free_y", nrow = 3)+
  theme_bw()+
  ylab("Biomass (tons)")+
  ggtitle("EBS snow estimated biomass") +
  scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
  scale_fill_manual(values = c("salmon",  "turquoise"), labels = c("50", "120"), name = "Knots")+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        title = element_text(size = 16))

ggsave(plot = bio.ind.plot.EBS, "./BAIRDI/Figures/TannerEBS.biomass.index.png", height= 9, width = 7.5, units = "in")

