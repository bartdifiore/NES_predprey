
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
saveRDS(fit_delta, paste0(proj_box_path, "Predator_prey/Github_LFS/Fitted_Mods/mac_fit_delta.rds"))



library(dplyr)
library(ggplot2)
library(viridis)  # for color scales

# 1. Marginal effect plots -------------------------------------

time_vals <- unique(mod_data$time)

time_vals <- c(seq(1970.25,2024.25, by = 5), 2019.25, 2024.25)[c(seq(1970.25,2024.25, by = 5), 2019.25, 2024.25) != 2020.25]

# Create newdata grid

temp_vals <- mod_data %>% filter(time %in% time_vals) %>% 
  group_by(time, scaled_time) %>%
  summarize(cold = quantile(BT_seasonal_scaled, 0.1), 
         median = quantile(BT_seasonal_scaled, 0.5), 
         hot = quantile(BT_seasonal_scaled, 0.9)) %>% 
  pivot_longer(cold:hot) %>%
  rename(temp_cat = name, BT_seasonal_scaled = value) %>% 
  mutate(season = "Spring", 
         survey = "NEFSC")
  

newdata <- crossing(
  temp_vals,
  dog_biomass_log_scaled = seq(min(mod_data$dog_biomass_log_scaled), max(mod_data$dog_biomass_log_scaled), length.out = 50)
)



pred_occ <- predict(fit_delta, newdata = newdata, model = 1, re_form = NA)
pred_occ_df <- pred_occ %>% 
  rename(occ = est1, biom = est2)


# Plot presence probability
p1 <- pred_occ_df %>%
  ggplot(aes(x = dog_biomass_log_scaled, y = occ,
                              color = temp_cat)) +
  geom_line(size = 1.2) +
  scale_color_viridis_d(name = "Bottom Temp\n(scaled)") +
  labs(
    x = "Dogfish Biomass (log scaled)",
    y = "Predicted Mackerel Presence Probability",
    title = "Marginal Effects on Mackerel Presence"
  ) +
  facet_wrap(~time)+
  theme_minimal()

# Plot biomass (conditional on presence)
p2 <- pred_occ_df %>% 
  ggplot(aes(x = dog_biomass_log_scaled, y = biom,
                               color = temp_cat)) +
  geom_line(size = 1.2) +
  scale_color_viridis_d(name = "Bottom Temp\n(scaled)") +
  labs(
    x = "Dogfish Biomass (log scaled)",
    y = "Predicted Mackerel Biomass (given presence)",
    title = "Marginal Effects on Mackerel Biomass"
  ) +
  facet_wrap(~time)+
  theme_minimal()

print(p1)
print(p2)

ggsave("Figures/dogfish_effects_mackerel_presence.png", p1)
ggsave("Figures/dogfish_effects_mackerel_biomass.png", p2)


# -------------------------------------------------------------

# 2. Spatial prediction maps -----------------------------------

# Use unique spatial locations from your data for prediction
spatial_grid <- mod_data %>%
  distinct(x, y, dog_biomass_log_scaled, BT_seasonal_scaled, scaled_time, season, survey)

# Predict occurrence probability across space
spatial_occ <- predict(fit_delta, newdata = spatial_grid, model = 1, re_form = NA)
spatial_grid$presence_prob <- spatial_occ$est

# Predict positive biomass across space
spatial_biom <- predict(fit_delta, newdata = spatial_grid, model = 2, re_form = NA)
spatial_grid$biomass <- spatial_biom$est

# Calculate expected biomass = presence_prob * biomass
spatial_grid$expected_biomass <- spatial_grid$presence_prob * spatial_grid$biomass

# Plot presence probability map
p_occ_map <- ggplot(spatial_grid, aes(x = x, y = y, fill = presence_prob)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Presence Probability") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial Prediction: Mackerel Presence Probability")

# Plot biomass map (given presence)
p_biom_map <- ggplot(spatial_grid, aes(x = x, y = y, fill = biomass)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Positive Biomass") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial Prediction: Mackerel Biomass (Given Presence)")

# Plot expected biomass map
p_exp_map <- ggplot(spatial_grid, aes(x = x, y = y, fill = expected_biomass)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Expected Biomass") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial Prediction: Expected Mackerel Biomass")

print(p_occ_map)
print(p_biom_map)
print(p_exp_map)

# -------------------------------------------------------------

# 3. Time trends at cold/hot spots with low/high dogfish ------

cold_temp <- quantile(mod_data$BT_seasonal_scaled, 0.1)
hot_temp <- quantile(mod_data$BT_seasonal_scaled, 0.9)
low_dog <- quantile(mod_data$dog_biomass_log_scaled, 0.1)
high_dog <- quantile(mod_data$dog_biomass_log_scaled, 0.9)

time_seq <- seq(min(mod_data$scaled_time), max(mod_data$scaled_time), length.out = 50)

trend_df <- expand.grid(
  scaled_time = time_seq,
  BT_seasonal_scaled = c(cold_temp, hot_temp),
  dog_biomass_log_scaled = c(low_dog, high_dog),
  season = "Spring",
  survey = "NEFSC"
)

# Predict occurrence and biomass
trend_occ <- predict(fit_delta, newdata = trend_df, model = 1, re_form = NA)
trend_biom <- predict(fit_delta, newdata = trend_df, model = 2, re_form = NA)

trend_df$presence_prob <- trend_occ$est
trend_df$biomass <- trend_biom$est
trend_df$expected_biomass <- trend_df$presence_prob * trend_df$biomass

# Create descriptive labels
trend_df <- trend_df %>%
  mutate(
    Temp = ifelse(BT_seasonal_scaled == cold_temp, "Cold Temp", "Hot Temp"),
    Dogfish = ifelse(dog_biomass_log_scaled == low_dog, "Low Dogfish", "High Dogfish"),
    Group = paste(Dogfish, Temp)
  )

ggplot(trend_df, aes(x = scaled_time, y = expected_biomass, color = Group)) +
  geom_line(size = 1.2) +
  labs(
    x = "Scaled Time",
    y = "Expected Mackerel Biomass",
    color = "Condition",
    title = "Predicted Mackerel Biomass Trends\nat Hot/Cold Temp and Low/High Dogfish Biomass"
  ) +
  theme_minimal()
