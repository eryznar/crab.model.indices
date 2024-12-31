dir2 <- paste0(dir, "Data/VAST/")

# Tweedie
imfem.TW <- rbind(read.csv(paste0(dir2, "VAST_imfem_TW_abund50.csv")) %>%
                mutate(knots = 50, type = "abundance"),
               read.csv(paste0(dir2, "VAST_imfem_TW_abund120.csv")) %>%
                 mutate(knots = 120, type = "abundance"),
               read.csv(paste0(dir2, "VAST_imfem_TW_bio50.csv")) %>%
                 mutate(knots = 50, type = "biomass"),
               read.csv(paste0(dir2, "VAST_imfem_TW_bio120.csv")) %>%
                 mutate(knots = 120, type = "biomass")) %>%
          mutate(matsex = "Immature Female", family = "Tweedie")

matfem.TW <- rbind(read.csv(paste0(dir2, "VAST_matfem_TW_abund50.csv")) %>%
                    mutate(knots = 50, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_matfem_TW_abund120.csv")) %>%
                    mutate(knots = 120, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_matfem_TW_bio50.csv")) %>%
                    mutate(knots = 50, type = "biomass"),
                  read.csv(paste0(dir2, "VAST_matfem_TW_bio120.csv")) %>%
                    mutate(knots = 120, type = "biomass")) %>%
              mutate(matsex = "Mature Female", family = "Tweedie")


male.TW <- rbind(read.csv(paste0(dir2, "VAST_male_TW_abund50.csv")) %>%
                     mutate(knots = 50, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_male_TW_abund120.csv")) %>%
                     mutate(knots = 120, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_male_TW_bio50.csv")) %>%
                     mutate(knots = 50, type = "biomass"),
                   read.csv(paste0(dir2, "VAST_male_TW_bio120.csv")) %>%
                     mutate(knots = 120, type = "biomass")) %>%
            mutate(matsex = "Male", family = "Tweedie")

# Delta gamma
imfem.DG <- rbind(read.csv(paste0(dir2, "VAST_imfem_DG_abund750.csv")) %>%
                    mutate(knots = 750, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_imfem_DG_bio750.csv")) %>%
                    mutate(knots = 750, type = "biomass")) %>%
            mutate(matsex = "Immature Female", family = "Delta gamma")

matfem.DG <- rbind(read.csv(paste0(dir2, "VAST_matfem_DG_abund750.csv")) %>%
                    mutate(knots = 750, type = "abundance"),
                  read.csv(paste0(dir2, "VAST_matfem_DG_bio750.csv")) %>%
                    mutate(knots = 750, type = "biomass")) %>%
              mutate(matsex = "Mature Female", family = "Delta gamma")

male.DG <- rbind(read.csv(paste0(dir2, "VAST_male_DG_abund750.csv")) %>%
                     mutate(knots = 750, type = "abundance"),
                   read.csv(paste0(dir2, "VAST_male_DG_bio750.csv")) %>%
                     mutate(knots = 750, type = "biomass")) %>%
            mutate(matsex = "Male", family = "Delta gamma")

# Bind all
VAST.24 <- rbind(imfem.TW, matfem.TW, male.TW, imfem.DG, matfem.DG, male.DG) 

VAST.abund <- VAST.24 %>%
  filter(type == "abundance") %>%
  mutate(Estimate = Estimate/1e6, Std..Error.for.Estimate = Std..Error.for.Estimate/1e6) %>%
  rename(SE = Std..Error.for.Estimate, Year = Time)

VAST.bio <- VAST.24 %>%
  filter(type == "biomass") %>%
  mutate(Estimate = Estimate/1000, Std..Error.for.Estimate = Std..Error.for.Estimate/1000) %>%
  rename(SE = Std..Error.for.Estimate, Year = Time)

# Read in model indices
TW.abund <- read.csv(paste0(dir, "Output/Tweedie_abund_index.csv")) %>%
  filter(knots != 90) %>% 
  mutate(model = "sdmTMB")
TW.bio <- read.csv(paste0(dir, "Output/Tweedie_bio_index.csv")) %>%
  filter(knots != 90) %>%
  mutate(model = "sdmTMB")

# Observations
tan.obs %>%
  group_by(Year, type, matsex) %>%
  reframe(value = sum(value),
          CI = sum(CI)) -> tan.obs3

# Plot tweedie
VAST.abund %>% 
  filter(family == "Tweedie", Stratum == "Stratum_3") %>% 
  mutate(model = "VAST",
         lwr = exp(log(Estimate) + qnorm(0.025) * SE), 
         upr = exp(log(Estimate) + qnorm(0.975) * SE)) -> VAST.aa
VAST.bio %>% 
  filter(family == "Tweedie", Stratum == "Stratum_3") %>%
  mutate(model = "VAST",
         lwr = exp(log(Estimate) + qnorm(0.025) * SE), 
         upr = exp(log(Estimate) + qnorm(0.975) * SE)) -> VAST.bb

ggplot()+
  geom_ribbon(TW.abund, mapping = aes(x = Year, ymin = lwr, ymax = upr, fill = model), alpha = 0.4)+
  geom_line(TW.abund, mapping = aes(Year, abundance, color = model))+
  geom_point(tan.obs3 %>% filter(type == "abundance"),
             mapping = aes(Year, value), color = "grey20", size = 0.75)+
  geom_errorbar(tan.obs3 %>% filter(type == "abundance"),
                mapping = aes(x = Year, ymin = value - CI, ymax = value+CI), color = "grey20", width = 0)+
  geom_errorbar(VAST.aa, mapping = aes(x = Year, ymin = Estimate - SE, ymax = Estimate + SE, color = model))+
  geom_line(VAST.aa, mapping = aes(Year, Estimate, color = model))+
  facet_grid(factor(knots, levels = c(50, 120)) ~factor(matsex, levels = c("Male", "Mature Female", "Immature Female")), scales = "free_y")+
  theme_bw()+
  ylab("Abundance (millions)")+
  ggtitle("EBS Tanner estimated abundance") +
  # scale_color_manual(values = c("salmon", "turquoise"), labels = c("50", "120"), name = "Knots")+
  # scale_fill_manual(values = c("salmon", "turquoise"), labels = c("50","120"), name = "Knots")+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        title = element_text(size = 16))
