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


#------------------------------------------------------------
## Plot effect of dogfish on mackerel
#------------------------------------------------------------

fit_wdogs <- readRDS(paste0(proj_box_path, "Predator_prey/Github_LFS/Fitted_Mods/mac_fit_wdogs.rds"))
pred_mod <- readRDS(paste0(proj_box_path, "Predator_prey/Github_LFS/Fitted_Mods/dog_fit_base.rds"))


# In order to have a continuous prediction surface, I need to use the predictions from the predator model to build the prediction grid for this model.

prediction_grid <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds")) %>%
  select(-c(BT_seasonal, SST_seasonal)) %>%
  filter(year %in% year_filt) %>%
  rename(longitude = x, latitude = y) %>%
  mutate(time = as.numeric(case_when(season == "Spring" ~ paste(year, "25", sep = "."), 
                                     season == "Summer" ~ paste(year, "50", sep = "."), 
                                     season == "Fall" ~ paste(year, "75", sep = "."))), 
         scaled_time = (time - mean(time))/10) %>% 
  filter(Depth < 400) %>% 
  select(-Depth, -Year_Season) %>%
  as_tibble() %>%
  add_utm_columns(ll_names = c("longitude", "latitude")) %>%
  mutate(survey = "NEFSC") %>% 
  filter(!time %in% c(2024.50, 2024.75))


pred_yhat <- predict(pred_mod, newdata = prediction_grid)

prediction_grid$pred_yhat <- exp(pred_yhat$est) 
prediction_grid$dog_biomass_log_scaled <- scale(log(prediction_grid$pred_yhat + 1))


fitlp_preds <- predict(fit_wdogs, newdata = prediction_grid)

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)



fitlp_preds2 <- fitlp_preds %>%
  as_tibble() %>%
  mutate(temp_id = paste(longitude, latitude, sep = "_"))

formerge <- data.frame(temp_id = unique(fitlp_preds2$temp_id))
formerge$id <- as.factor(1:length(formerge$temp_id))

fitlp_preds2 %>%
  rename(interaction = `zeta_s_dog_biomass_log_scaled:scaled_time`) %>%
  filter(time == 2023.75, season == "Fall") %>% 
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = interaction)) +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_fill_viridis_c() +
  theme_minimal()+
  labs(
    fill = "Interaction", 
    title = "Temporal effect of spiny dogfish on mackerel", x = "", y = ""
  )
ggsave("Figures/temporal_effects.png")


fitlp_preds2 %>%
  left_join(formerge) %>% 
  filter(id %in% sample(levels(id),9)) %>%
  ggplot(aes(x = time, y = est))+
  geom_line(aes(color = as.numeric(id), linetype = season), show.legend = F)+
  facet_wrap(~id)

fitlp_preds2 %>%
  left_join(formerge) %>% 
  filter(id %in% sample(levels(id),9)) %>%
  ggplot(aes(x = dog_biomass_log_scaled, y = est))+
  geom_line(aes(color = as.numeric(id), linetype = season), show.legend = F)+
  facet_wrap(~id)


obj <- predict(fit_wdogs, return_tmb_object = TRUE)

# This contains the spatially varying slopes
svc_matrix <- obj$report$zeta_s_A[, , 1]  # Drop 3rd dim
colnames(svc_matrix) <- c("dog_biomass_log_scaled", "scaled_time", "interaction")
svc_df <- cbind(
  mod_data[, c("X", "Y")],  # Assuming these are your coordinates
  svc_matrix
)
library(ggplot2)

ggplot(svc_df, aes(x = X, y = Y, fill = interaction)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Spatially Varying Coefficient: dog_biomass_log_scaled × scaled_time",
    fill = "Slope"
  ) +
  theme_minimal()


ggplot(svc_df, aes(x = X, y = Y, color = interaction)) +
  geom_point(size = 1) +
  coord_fixed() +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "SVC: dog_biomass_log_scaled × scaled_time",
    color = "Slope"
  ) +
  theme_minimal()





fit_delta <- sdmTMB(
  mac_biomass ~ season + survey + scaled_time*dog_biomass_log_scaled*BT_seasonal_scaled,
  data = mod_data,
  family = delta_gamma(link1 = "logit", link2 = "log"),
  mesh = sdmTMB_mesh,
  spatial = "on",
  spatiotemporal = "IID",
  time = "time",
  anisotropy = FALSE,
  #spatial_varying = ~ 0,  # if you're not using SVCs here
  silent = FALSE
)

sanity(fit_delta)































mod_data


fit_wdogs_fixed <- sdmTMB(
  mac_biomass ~ season + survey + scaled_time * dog_biomass_log_scaled * BT_seasonal_scaled,
  data = mod_data,
  spatial = "on",
  spatiotemporal = "IID",
  time = "time",
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE,
  do_fit = TRUE
)


fit_dog_svc <- sdmTMB(
  mac_biomass ~ season + survey + dog_biomass_log_scaled * scaled_time * BT_seasonal_scaled,
  data = mod_data,
  spatial = "on",
  spatiotemporal = "IID",
  spatial_varying = ~ dog_biomass_log_scaled:scaled_time:BT_seasonal_scaled,
  time = "time",
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE
)
sanity(fit_dog_svc)
tidy(fit_dog_svc, effects = "fixed")
tidy(fit_dog_svc, effects = "ran_pars")

write_rds(fit_dog_svc, "Data/Derived/fit_dog_svc.rds", compress = "gz")


# Confirm which index corresponds to the 3-way interaction term in zeta_s_A:
# e.g., names and dim(obj$report$zeta_s_A)

obj <- predict(fit_dog_svc, return_tmb_object = TRUE)

str(obj$report$zeta_s_A)  # check names if available

# Extract the 3-way interaction SVC (assuming it's the first or a known slice)
svc_3way <- obj$report$zeta_s_A[, , 1]  # adjust index if needed

# Create dataframe with coordinates and svc values
svc_df <- mod_data %>%
  select(X, Y) %>%
  mutate(interaction_3way = svc_3way)  # confirm slicing for correct dimension

# Plot
ggplot(svc_df, aes(x = X, y = Y, color = interaction_3way)) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(option = "plasma") +
  coord_fixed() +
  labs(
    title = "Spatially Varying 3-Way Interaction: dogfish × time × temperature",
    color = "Slope"
  ) +
  theme_minimal()


# # | Color                                                         | Meaning                                                                                                                                                                                                                                                                              |
# | ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
#   | **Negative values (e.g., cool colors like blues, purples)**   | **At these locations, warming amplifies the negative effect of dogfish on mackerel over time.** <br> In other words, higher temperatures intensify dogfish predation impact on mackerel, causing a stronger decline or slower growth of mackerel biomass where dogfish are abundant. |
#   | **Near zero (neutral colors, e.g., yellows or light colors)** | Temperature does not significantly modify the dogfish × time effect here. <br> Dogfish impact on mackerel over time is stable regardless of temperature changes.                                                                                                                     |
#   | **Positive values (warm colors like reds, oranges)**          | At these locations, warming *reduces* or even *reverses* the negative dogfish effect on mackerel trends. <br> This could mean warming weakens predation pressure or enhances mackerel resilience despite dogfish presence.                                                           |
# Find extreme negative and positive values



svc_df <- mod_data %>%
  # select(X, Y) %>%
  mutate(interaction_3way = svc_3way, 
         est = obj$report$eta_i) %>%  # confirm slicing for correct dimension
  mutate(type = case_when(
    interaction_3way <= quantile(interaction_3way, 0.05) ~ "strong negative",
    interaction_3way >= quantile(interaction_3way, 0.95) ~ "strong positive",
    TRUE ~ NA_character_
  )) %>% 
  drop_na(type)

ggplot(svc_df, aes(y = est, x = dog_biomass_log_scaled))+
  geom_point(aes(color= type))



svc_extremes <- svc_df %>%
  mutate(type = case_when(
    interaction_3way <= quantile(interaction_3way, 0.05) ~ "strong negative",
    interaction_3way >= quantile(interaction_3way, 0.95) ~ "strong positive",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(type))

# Suppose you pick the top hotspot and coldspot by value
hotspot_coords <- svc_extremes %>% filter(type == "hotspot") %>% slice_max(interaction_3way, n = 1)
coldspot_coords <- svc_extremes %>% filter(type == "coldspot") %>% slice_min(interaction_3way, n = 1)

# Set up a grid of years (e.g., standardized times from your model)
time_grid <- unique(mod_data$time)

# Example: fixed values for dogfish and BT
fixed_dog <- seq(min(mod_data$dog_biomass_log_scaled), max(mod_data$dog_biomass_log_scaled), length.out = 30)  # mean
fixed_temp <- 0 # mean

# Combine into prediction data frames
make_trend_df <- function(coords, label) {
  expand.grid(
    time = time_grid,
    dog_biomass_log_scaled = fixed_dog,
    BT_seasonal_scaled = fixed_temp,
    season = "Spring",  # or whatever default
    survey = "NEFSC",   # match model factors
    X = coords$X,
    Y = coords$Y
  ) %>% mutate(location = label)
}

hot_df <- make_trend_df(hotspot_coords, "hotspot")
cold_df <- make_trend_df(coldspot_coords, "coldspot")

# Combine both
pred_df <- bind_rows(hot_df, cold_df)
pred_df$scaled_time <- (pred_df$time - mean(pred_df$time))/10
preds <- predict(fit_dog_svc, newdata = pred_df, return_tmb_object = FALSE)
pred_df$predicted <- preds$est

pred_df %>% 
  filter(time == 2023.25) %>%
  ggplot(aes(x = dog_biomass_log_scaled, y = predicted, color = location)) +
  geom_line(size = 1.2) +
  labs(
    title = "Predicted mackerel biomass trends",
    subtitle = "Hotspot vs. Coldspot based on 3-way interaction strength",
    x = "Scaled Time",
    y = "Predicted Biomass (on link scale or backtransformed)"
  ) +
  theme_minimal()

pred_df %>% 
  ggplot(aes(x = dog_biomass_log_scaled, y = predicted, color = time)) +
  geom_point()+
  labs(
    title = "Predicted mackerel biomass trends",
    subtitle = "Hotspot vs. Coldspot based on 3-way interaction strength",
    x = "Dogfish biomass",
    y = "Predicted Biomass (on link scale or backtransformed)"
  ) +
  facet_wrap(~location)+
  theme_minimal()
