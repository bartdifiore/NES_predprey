#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sdmTMB)
library(sdmTMBextra)
library(ggeffects)
library(patchwork)
source(here::here("Code/sdmTMB_validation_Functions.R"))
library(gmRi)


# Paths to Box folder
proj_box_path<- cs_path(box_group = "Mills Lab", subfolder = "Projects/NSF_NECC")

# Base map land info
region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

#----------------------------------
## Get data
#----------------------------------
macs <- readRDS(here::here("Data/Derived/all_model_data_mac.rds")) %>% 
  drop_na(year, total_biomass) %>% # Why does the year column have NA's???? 
  rename(mac_biomass = total_biomass)

dogs <- readRDS(here::here("Data/Derived/all_model_data_predators.rds")) %>% 
  drop_na(year) %>% #Why does the year column have NA's????
  rename(dog_biomass = total_biomass)

df <- left_join(macs, dogs)

year_filt <- c(seq(1970, 2016, by = 1),2018,2019,2021,2022,2023, 2024)

red_mod_data <- df |>
  dplyr::filter(year %in% year_filt) 

summary(red_mod_data)

# Cutting depth at 400 m
red_mod_data <- red_mod_data |>
  dplyr::filter(Depth <= 400)
summary(red_mod_data)

# Scale/center covariates
# Get means and sds
column_means <- colMeans(red_mod_data[, c("Depth", "BT_seasonal")], na.rm = TRUE)
column_sds <- apply(red_mod_data[, c("Depth", "BT_seasonal")], 2, sd, na.rm = TRUE)

# Scale the data
mod_data <- red_mod_data |>
  mutate(across(
    c(Depth, BT_seasonal),
    ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
    .names = "{.col}_scaled"
  ))
summary(mod_data)

mod_data <- mod_data %>%
  mutate(time = as.numeric(case_when(season == "Spring" ~ paste(year, "25", sep = "."), 
                                     season == "Summer" ~ paste(year, "50", sep = "."), 
                                     season == "Fall" ~ paste(year, "75", sep = "."))), 
         scaled_time = (time - mean(time))/10) %>%
  as.data.frame()

mod_data$dog_biomass_log_scaled <- scale(log(mod_data$dog_biomass + 1))

#----------------------------------
## Make mesh
#----------------------------------
# This is definitely something we will want to come back to after we have gotten things up and running through to model inferences. 
mod_data <- mod_data %>%
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") 

sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("X", "Y"),cutoff = 30) # added for speed...
plot(sdmTMB_mesh)

#------------------------------------------------------------
## Fit basic species distribution model w/out covariates
#------------------------------------------------------------

fit_wdogs <- sdmTMB(
  mac_biomass ~ season + survey + scaled_time,
  data = mod_data,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = F,
  # share_range = TRUE,
  spatiotemporal = "IID",
  spatial_varying = ~dog_biomass_log_scaled*scaled_time,
  time = "time",
  # time_varying = ~ 0 + year_season_int,
  # time_varying_type = "ar1",
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE,
  do_fit = TRUE
)

sanity(fit_wdogs)
fit_wdogs
tidy(fit_wdogs, effects = "ran_pars")

ggeffects::ggpredict(fit_wdogs, terms = "scaled_time")

saveRDS(fit_wdogs, paste0(proj_box_path, "Predator_prey/Github_LFS/Fitted_Mods/mac_fit_wdogs.rds"))
