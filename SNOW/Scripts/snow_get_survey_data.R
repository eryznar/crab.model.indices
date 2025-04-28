source("./SNOW/Scripts/snow_load_libs_functions.R")


# Connect to Oracle via AFSC database (NOAA AFSC users only)
channel <- "API"

# Pull specimen data ----
species <- "SNOW"
specimen_data_EBS <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = c(1975:2024),
                                             channel = channel) 

specimen_data_NBS <- crabpack::get_specimen_data(species = species,
                                                 region = "NBS",
                                                 years = c(1975:2024),
                                                 channel = channel) 
# saveRDS(specimen_data_NBS, "Y:/KOD_Research/Ryznar/Crab functional maturity/opilio.functional.maturity/Data/snow_specimen_NBS.rda")
# 
# spec.dat.EBSNBS <- rbind(specimen_data_EBS$specimen, specimen_data_NBS$specimen)
# haul.dat.EBSNBS <- rbind(specimen_data_EBS$haul %>% dplyr::select(!BOTTOM_TYPE), specimen_data_NBS$haul)
# 
# saveRDS(spec.dat.EBSNBS, paste0(dir, "Data/snow_survey_specimen_EBSNBS.rda"))
# saveRDS(haul.dat.EBSNBS, paste0(dir, "Data/snow_survey_haul_EBSNBS.rda"))

# Get all male cpue ----
male_cpue_EBS  <- calc_cpue(crab_data = specimen_data_EBS, 
                               sex = "male", 
                               region = "EBS", 
                               district = "ALL",
                               bin_1mm = TRUE,
                               species = "SNOW") 

male_cpue_NBS  <- calc_cpue(crab_data = specimen_data_NBS, 
                                    sex = "male", 
                                    region = "NBS", 
                                    district = "ALL",
                                    bin_1mm = TRUE,
                                    species = "SNOW") 

saveRDS(setDT(rbind(male_cpue_EBS, male_cpue_NBS)), paste0(dir, "Data/snow_survey_cpue_male_EBSNBS.rda"))

# Get mature female cpue ----
matfem_cpue_EBS  <- calc_cpue(crab_data = specimen_data_EBS, 
                            sex = "female", 
                            crab_category = "mature_female",
                            region = "EBS", 
                            district = "ALL",
                            species = "SNOW") 

matfem_cpue_NBS  <- calc_cpue(crab_data = specimen_data_NBS, 
                              sex = "female", 
                              crab_category = "mature_female",
                              region = "NBS", 
                              district = "ALL",
                              species = "SNOW")

saveRDS(setDT(rbind(matfem_cpue_EBS, matfem_cpue_NBS)), paste0(dir, "Data/snow_survey_cpue_matfem_EBSNBS.rda"))


# Get all male biomass ----
male_bio_EBS  <- calc_bioabund(crab_data = specimen_data_EBS, 
                            sex = "male", 
                            region = "EBS", 
                            district = "ALL",
                            bin_1mm = TRUE,
                            species = "SNOW") 

male_cpue_NBS  <- calc_cpue(crab_data = specimen_data_NBS, 
                            sex = "male", 
                            region = "NBS", 
                            district = "ALL",
                            bin_1mm = TRUE,
                            species = "SNOW") 

saveRDS(setDT(rbind(male_cpue_EBS, male_cpue_NBS)), paste0(dir, "Data/snow_survey_cpue_male_EBSNBS.rda"))
