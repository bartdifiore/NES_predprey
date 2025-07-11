library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)
library(sf)
library(gmRi)
library(rnaturalearth)
library(rnaturalearthdata)

#-----------------------------------------
## Get data
#-----------------------------------------

# Paths to Box folder
proj_box_path<- cs_path(box_group = "Mills Lab", subfolder = "Projects/NSF_NECC")

pred_mod <- readRDS(paste0(proj_box_path, "Predator_prey/Github_LFS/Fitted_Mods/mac_fit_base.rds"))
mac_mod <- readRDS(paste0(proj_box_path, "Predator_prey/Github_LFS/Fitted_Mods/dog_fit_base.rds"))


#-----------------------------------------
## Get prediction grid
#-----------------------------------------

year_filt <- c(seq(1970, 2016, by = 1),2018,2019,2021,2022,2023, 2024)

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

#-------------------------------------------------------
## Generate spatiotemporal predictions for each model
#-------------------------------------------------------

mac_yhat <- predict(mac_mod, newdata = prediction_grid)
pred_yhat <- predict(pred_mod, newdata = prediction_grid)

source("Code/Carroll_2019_functions.R")

# Combine into one df

df_overlap <- pred_yhat %>% 
  arrange(season, year, time, Date, longitude, latitude) %>%
  ungroup() %>%
  mutate(id = 1:n()) %>%
  select(id, season, year, time, est, longitude, latitude) %>%
  rename(predator_est = est) %>%
  st_drop_geometry() %>%
  left_join(mac_yhat %>% 
              arrange(season, year, time, Date, longitude, latitude) %>%
              ungroup() %>%
              mutate(id = 1:n()) %>%
              select(id,season, year, time, est, longitude, latitude) %>%
              rename(prey_est = est) %>%
              st_drop_geometry) #%>% 
  # mutate(temp = paste(latitude, longitude, sep = "_")) %>% 
  # group_by(temp) %>% 
  # mutate(id = 1:n())

#coords <- unique(df_overlap$temp)

# grid <- data.frame(temp = coords) %>% 
#   separate(temp, into = c("latitude", "longitude"), sep = "_") %>% 
#   mutate(id = 1:n())%>%
#   sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
#   st_set_crs(4326) %>%
#   st_buffer(dist = 1/12)

#-------------------------------------------------------
## Visualize overlaps using RBG mapping technique
#-------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(scales)


# Step 1: Rescale within each year

df_scaled <- df_overlap %>%
  mutate(
    # Option 1: rescale 0–1 for visualization
    predator_scaled = rescale(predator_est),
    prey_scaled     = rescale(prey_est),
    rgb_scaled = rgb(red = predator_scaled, green = 0, blue = prey_scaled), 
    # rgb_scaled = rgb(red = predator_scaled, green = prey_scaled, blue = 0), 
    decade = floor(time/10)*10)
  # ) %>%
  # left_join(grid) %>% 
  # st_as_sf()


region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

# Step 2: Plot as raster with layering over time

out <- ggplot(df_scaled) +
  geom_tile(aes(fill = rgb_scaled, x = longitude, y = latitude), alpha = 0.95) +  # NOTE: geom_tile() also works
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_fill_identity() +
  # coord_fixed() +
  labs(
    title = "Cumulative Overlap of Dogfish (Red) and Mackerel (Blue)",
    subtitle = "Purple = persistent overlap; alpha layering across years", x = "", y = ""
  ) +
  facet_wrap(~decade)+
  theme_minimal()

ggsave("Figures/dogfish_mackerel_rbgoverlap.png", out)




#------------------------------------------------------------------------------
## Estimate overlap metrics for each cell (e.g. do not integrate across space)
#------------------------------------------------------------------------------

grid <- df_overlap %>%
  distinct(latitude, longitude) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>%
  st_buffer(dist = 1/12) %>%
  mutate(cell_id = 1:n(), 
         area = st_area(geometry))

plot(grid)

df_overlap2 <- df_overlap %>% 
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>% 
  st_join(grid) %>%
  mutate(area_km2 = as.numeric(area)/1000000)

overlap <- df_overlap2 %>%
  group_by(id) %>%
  mutate(prey_est_response = exp(prey_est), 
         predator_est_response = exp(predator_est))

overlap$prey_bin <- ifelse(overlap$prey_est_response <= quantile(overlap$prey_est_response, 0.05), 0, 1)
overlap$predator_bin <- ifelse(overlap$predator_est_response <= quantile(overlap$predator_est_response, 0.05), 0, 1)

overlap$co_occurance <- ifelse(overlap$prey_bin == 1 & overlap$predator_bin == 1, 1, 0)

overlaps <- overlap %>% 
  group_by(year, season) %>%
  mutate(p_prey = prey_est_response/sum(prey_est_response), 
         p_predator = predator_est_response/sum(predator_est_response), 
         asymmalpha = p_prey*p_predator/sum(p_prey^2), 
         loc_colloc = p_prey*p_predator/(sqrt(sum(p_predator^2)*sum(p_prey^2))), 
         biomass_overlap = (prey_est_response/max(prey_est_response))*(predator_est_response/max(predator_est_response))/sum(prey_est_response/max(prey_est_response)), 
         schoeners = 1 - 0.5 * (abs(p_prey - p_predator)), 
         bhatta = sqrt(p_prey*p_predator), 
         AB = (predator_est_response - mean(predator_est_response))*(prey_est_response - mean(prey_est_response))/(mean(predator_est_response)*mean(prey_est_response)), 
         pred_prey_ratio = predator_est_response/prey_est_response) %>% 
  pivot_longer(cols = c(co_occurance, asymmalpha:pred_prey_ratio), names_to = "overlap_metric", values_to = "value" )



write_rds(overlaps, "Data/Derived/overlap_metrics.rds", compress = "gz")


#------------------------------------------------------------------------------
## Estimate overlap metrics for larger spatial areas
#------------------------------------------------------------------------------

# Build a 1 degree grid, overlay the observations, estimate metrics 

grid.1 <- st_bbox(overlap) %>%
  st_as_sfc() %>%
  st_make_grid(cellsize = 1) %>%
  st_sf() %>%
  mutate(cell_id = 1:n(), 
         area = st_area(geometry), 
         area_km2 = as.numeric(area/1000000)) %>%
  select(-area) %>%
  st_make_valid()

ggplot() +
  geom_sf(data = grid.1) +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)

# Merge the data to the 1 degree grid

df_overlap3 <- df_overlap %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_join(grid.1)

# Estimate the metrics

overlaps_1deg <- df_overlap3 %>% 
  st_drop_geometry() %>%
  group_by(season, year, cell_id, area_km2) %>%
  mutate(p_res = exp(predator_est), # Here I just use the predicted estimates on the response scale. However, should we use get_index() to estimate relative biomas in each spatial unit???
         l_res = exp(prey_est)) # %>%
# mutate(p_res = exp(predator_est), # Here I just use the predicted estimates on the response scale. However, should we use get_index() to estimate relative biomas in each spatial unit???
#        l_res = exp(lobster_est)
# )

overlaps_1deg$prey_bin <- ifelse(overlaps_1deg$l_res <= quantile(overlaps_1deg$l_res, 0.05), 0, 1)
overlaps_1deg$predator_bin <- ifelse(overlaps_1deg$p_res <= quantile(overlaps_1deg$p_res, 0.05), 0, 1)

out_1_deg <- overlaps_1deg %>% 
  group_by(season, year, cell_id) %>%
  summarize(area_overlap = area_overlapfn(l_res, p_res, area_km2),
            range_overlap = range_overlapfn(l_res, p_res, area_km2),
            asymmalpha = asymmalpha_overlapfn(l_res, p_res), 
            loc_colloc = loc_collocfn(l_res, p_res), 
            biomass_overlap = biomass_overlapfn(l_res, p_res),
            hurlbert = hurlbert_overlapfn(l_res, p_res, area_km2),
            schoeners = schoeners_overlapfn(l_res, p_res), 
            bhatta = bhatta_coeffn(l_res, p_res), 
            AB = AB_overlapfn(l_res, p_res)) %>%
  pivot_longer(cols = c(area_overlap:AB), names_to = "overlap_metric", values_to = "value" ) %>%
  left_join(grid.1)

summary(out_1_deg)

write_rds(out_1_deg, "Data/Derived/overlap_metrics_1deg.rds")  
