#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

glorys_df <- readRDS(here::here("Data/Derived/glorys_grid.rds"))

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = "United States of America", returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-75.8, -56.2)

#-------------------------------------------------------------
## Merge biological and covariate data for prey
#-------------------------------------------------------------
catch_data<- readRDS(here::here("Data/Derived/mackerel.rds"))
str(catch_data)

# Need total biomass by tow/life classm
mac_df_bio <- catch_data |>
  mutate(total_weight_at_length = number_at_length * weight_at_length) |>
  group_by(scientific_name, trawl_id, longitude, latitude, season, year, survey, date) |>
  summarize("total_biomass" = sum(total_weight_at_length)) |>
  ungroup()
summary(mac_df_bio)

env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)
env_tows <- unique(env_data$ID)
summary(env_data)

all(mac_df_bio$trawl_id %in% env_tows)

all_mod_data <- mac_df_bio |>
  # Adjust the "true" NA's before adding implicit NA values
  replace_na(list(total_biomass = 99999)) |>
  full_join(env_data, by = c("longitude" = "DECDEG_BEGLON", "latitude" = "DECDEG_BEGLAT", "trawl_id" = "ID", "year" = "EST_YEAR", "season" = "season", "date" = "DATE", "survey" = "survey")) |>
  # Make NAs 0, and then change 99999 fill to NAs
  replace_na(list(total_biomass = 0)) |>
  mutate(total_biomass = na_if(total_biomass, 99999))
summary(all_mod_data)

write_rds(all_mod_data, "Data/Derived/all_model_data_mac.rds", compress = "gz")

#---------------------------------------------------
## Merge biological and covariate data for predators
#---------------------------------------------------
catch_data<- readRDS(here::here("Data/Derived/combined_and_filtered_predators.rds"))
str(catch_data)

# Need total biomass by tow/life class
pred_df_bio <- catch_data |>
  mutate(total_weight_at_length = number_at_length * weight_at_length) |>
  group_by(trawl_id, longitude, latitude, season, year, survey, date) |>
  summarize("total_biomass" = sum(total_weight_at_length, na.rm = T))
summary(pred_df_bio)

env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)

all_mod_data <- pred_df_bio |>
  # Adjust the "true" NA's before adding implicit NA values
  replace_na(list(total_biomass = 99999)) |>
  full_join(env_data, by = c("longitude" = "DECDEG_BEGLON", "latitude" = "DECDEG_BEGLAT", "trawl_id" = "ID", "year" = "EST_YEAR", "season" = "season", "date" = "DATE", "survey" = "survey")) |>
  # Make NAs 0, and then change 99999 fill to NAs
  replace_na(list(total_biomass = 0)) |>
  mutate(total_biomass = na_if(total_biomass, 99999))
summary(all_mod_data)

write_rds(all_mod_data, "Data/Derived/all_model_data_predators.rds", compress = "gz")
