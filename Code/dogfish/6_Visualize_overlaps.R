library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#------------------------------------------------------------
## Get data
#------------------------------------------------------------

df <- readRDS("Data/Derived/overlap_metrics_1deg.rds") %>%
  st_as_sf()

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

df %>% 
  st_drop_geometry() %>%
  group_by(overlap_metric) %>% 
  summarize(low = quantile(value, 0.025, na.rm = T), 
            mean = mean(value, na.rm = T), 
            high = quantile(value, 0.975, na.rm = T) 
  )


#------------------------------------------------------------
## Visualize overlap metrics
#------------------------------------------------------------


df %>%
  st_make_valid() %>%
  filter(year %in% c(1975,1985,1995,2005, 2015, 2023), 
         season == "Fall", 
         overlap_metric == "biomass_overlap") %>%
  ggplot() +
  geom_sf(aes(fill = value), color = "transparent") +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  scale_fill_viridis_c(trans = "log") +
  facet_wrap(~year)+
  theme_minimal()


df %>%
  st_make_valid() %>%
  filter(year %in% c(1975,1985,1995,2005, 2015, 2023), 
         season == "Fall", 
         overlap_metric == "asymmalpha") %>%
  ggplot() +
  geom_sf(aes(fill = value), color = "transparent") +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  scale_fill_gradient(low = "white", high = scales::muted("red")) +
  facet_wrap(~year)+
  labs(title = "Asymetrical alpha")+
  theme_minimal()
ggsave("Figures/asymetricala_paneled_by_year.png")

df %>%
  filter(cell_id == 149, 
         season == "Fall", 
         overlap_metric == "asymmalpha") %>%
  ggplot() +
  geom_line(aes(x = year, y = value))

metrics <- unique(df$overlap_metric)[3:9]
labels <- c("Asymetrical alpha", "Local colocation", "Biomass overlap", "Hurlbert's overlap", "Schoener's D", "Bhattacharyya's coefficient", "AB ratio")
out <- list()
for(i in 1:length(metrics)){
  temp <- df %>%
    filter(overlap_metric == metrics[i],
           year %in% c(1975,1985,1995,2005, 2015, 2023)) %>%
    ggplot() +
    geom_sf(aes(fill = value), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    scale_fill_viridis_c() +
    labs(title = labels[i])+
    facet_grid(season~year)+
    theme_minimal()
    
  ggsave(paste("Figures/", labels[i], "_fig.png", sep = ""), plot = temp)
}


# What I'm interested in visualizing is the relative change through time. So for each grid cell has the value of asymmetrical alpha increased or decreased since 1993? To do this we will need to fit a lm() to each time series (or simply estimate a correlation coefficient) and then visualize the slope of this relationship. 

cells <- unique(df$cell_id)[is.na(unique(df$cell_id)) == F]

grid.1 <- df %>% 
  select(geometry, cell_id) %>%
  group_by(cell_id) %>% 
  summarize(first(cell_id)) %>% 
  select(-`first(cell_id)`)


one <- df %>% 
  ungroup() %>%
  filter(overlap_metric == "schoeners") %>%
  filter(cell_id %in% sample(cell_id, 3))

ggplot(one, aes(x = time, y = value))+
  geom_point()+
  facet_wrap(~cell_id + season, scales = "free")


forplot <- df %>% 
  drop_na(cell_id) %>%
  filter(!overlap_metric %in% c("area_overlap", "range_overlap")) %>%
  st_drop_geometry() %>%
  group_by(cell_id, season, overlap_metric) %>% 
  nest() %>% 
  mutate(lm_obj = map(data, ~lm(value ~ year, data = .))) %>%
  mutate(lm_tidy = map(lm_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(cell_id, season, overlap_metric, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  filter(term == "year") %>% 
  select(cell_id, season, overlap_metric, estimate) %>%
  left_join(grid.1)


metrics <- unique(forplot$overlap_metric)
labels <- c("Asymetrical alpha", "Local colocation", "Biomass overlap", "Hurlbert's overlap", "Schoener's D", "Bhattacharyya's coefficient", "AB ratio")
out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    filter(season == "Spring") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/Temporaltrends_in_overlap_spring.png")  

ggsave("Figures/Temporaltrends_schoenersD_spring.png", out[[5]])

out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    filter(season == "Fall") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/Temporaltrends_in_overlap_fall.png")

out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    filter(season == "Summer") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/Temporaltrends_in_overlap_summer.png")


forplot %>% 
  filter(overlap_metric == "biomass_overlap", 
         season == "Fall") %>% 
  ggplot()+
  geom_sf(aes(fill = estimate, geometry = geometry), color = "transparent") +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
  # scale_fill_viridis_c() +
  scale_fill_gradient2(high = scales::muted("red"),
                       mid = "white",
                       low = scales::muted("blue")) +
  labs(title = "Biomass overlap")+
  theme_minimal()
ggsave("Figures/biomass_trends.png")  

#-----------------------------------------------------------
# Are changes in overlap correlated with temperature? 
#-----------------------------------------------------------

# Sample one cell to determine relationship. Then apply to all cells. 

one <- df %>% 
  ungroup() %>%
  filter(time > 1992.9, 
         overlap_metric == "schoeners") %>%
  filter(cell_id %in% sample(cell_id, 3))

ggplot(one, aes(x = mean_sstemp, y = value))+
  geom_point()+
  facet_wrap(~cell_id + season, scales = "free")

one %>%
  filter(season == "Fall") %>%
ggplot(aes(x = time))+
  geom_line(aes(y = value), color = "darkgreen")+
  geom_line(aes(y = scales::rescale(mean_sstemp)), color = "darkblue")+
  facet_wrap(~cell_id)

# So really important to deal with season here. Because when looking across the year fall is always >>> spring for both species. So shows positive relationships across all spatial locations. Here I run separate regressions for each season, and plot.


forplot <- df %>% 
  drop_na(cell_id, mean_sstemp) %>%
  filter(!overlap_metric %in% c("area_overlap", "range_overlap")) %>%
  st_drop_geometry() %>%
  group_by(cell_id, season, overlap_metric) %>% 
  nest() %>% 
  mutate(lm_obj = map(data, ~lm(value ~ mean_sstemp, data = .))) %>%
  mutate(lm_tidy = map(lm_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(cell_id, season, overlap_metric, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  filter(term == "mean_sstemp") %>% 
  select(cell_id, season, overlap_metric, estimate) %>%
  left_join(grid.1)


metrics <- unique(forplot$overlap_metric)
labels <- c("Asymetrical alpha", "Local colocation", "Biomass overlap", "Hurlbert's overlap", "Schoener's D", "Bhattacharyya's coefficient", "AB ratio")
out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    filter(season == "Fall") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/temperature_effects_on_overlap_fall.png")  


out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    filter(season == "Spring") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/temperature_effects_on_overlap_spring.png") 

out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    filter(season == "Summer") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/temperature_effects_on_overlap_summer.png") 









df %>% 
  ggplot(aes(x = mean_btemp, y = value))+
  geom_point()+
  facet_wrap(~overlap_metric, scales = "free")+
  theme_minimal()

df %>% 
  filter(overlap_metric == "schoeners") %>%
  filter(cell_id %in% sample(cell_id, 5)) %>%
  ggplot(aes(x = time, group = cell_id))+
  geom_line(aes(y = value), color = "green")+
  geom_line(aes(y = mean_btemp), color = "darkred")+
  theme_minimal()




library(dplyr)
library(mgcv)

# Filter to one overlap metric
df_filtered <- df %>%
  filter(overlap_metric == "schoeners", 
         time > 1992) %>% 
  st_drop_geometry() %>%
  drop_na(mean_btemp, mean_sstemp, cell_id) %>%
  mutate(cell_id = as.factor(cell_id)) %>% 
  filter(value < 1)


# Fit GAM with a smooth effect of temperature and a random intercept by cell_id
model <- gam(value ~ s(mean_btemp) + s(cell_id, bs = "re"), data = df_filtered)
model <- gamm4(value ~ s(mean_btemp) + s(mean_sstemp), data = df_filtered, family =  random = ~(1|cell_id))
model <- sdmTMB(value ~ s(mean_btemp) + s(mean_sstemp) + (1|cell_id),
                data = df_filtered,
                family = Beta(link = "logit"), 
                spatial = "off",
                spatiotemporal = "off",
                #time = "time"
                )

# Check summary
summary(model)



# Step 1: Create a new data frame for predictions over a range of btemp
btemp_seq <- seq(min(df_filtered$mean_btemp, na.rm = TRUE),
                 max(df_filtered$mean_btemp, na.rm = TRUE),
                 length.out = 100)

# Hold sstemp constant at median
sstemp_fixed <- median(df_filtered$mean_sstemp, na.rm = TRUE)

# Use the most common cell_id to satisfy model structure
ref_cell <- sample(model$data$cell_id, 5)

newdata_btemp <- expand.grid(
  mean_btemp = btemp_seq,
  mean_sstemp = sstemp_fixed,
  cell_id = as.factor(ref_cell)
)

# Step 2: Predict from model
pred_btemp <- predict(model, newdata = newdata_btemp, se_fit = TRUE)


sstemp_seq <- seq(min(df_filtered$mean_sstemp, na.rm = TRUE),
                  max(df_filtered$mean_sstemp, na.rm = TRUE),
                  length.out = 100)

btemp_fixed <- median(df_filtered$mean_btemp, na.rm = TRUE)

newdata_sstemp <- expand.grid(
  mean_btemp = btemp_fixed,
  mean_sstemp = sstemp_seq,
  cell_id = as.factor(ref_cell)
)

pred_sstemp <- predict(model, newdata = newdata_sstemp, se_fit = TRUE)

library(ggplot2)

# Plot for mean_btemp
ggplot(pred_btemp, aes(x = mean_btemp, y = plogis(est), group = cell_id)) +  # plogis for inverse-logit
  geom_line() +
  geom_ribbon(aes(ymin = plogis(est - 1.96 * est_se),
                  ymax = plogis(est + 1.96 * est_se)), alpha = 0.2) +
  labs(x = "Mean Bottom Temperature",
       y = "Predicted Overlap (on response scale)",
       title = "Marginal Effect of Bottom Temp") +
  theme_minimal()

# Plot for mean_sstemp
ggplot(pred_sstemp, aes(x = mean_sstemp, y = plogis(est), group = cell_id)) +
  geom_line() +
  geom_ribbon(aes(ymin = plogis(est - 1.96 * est_se),
                  ymax = plogis(est + 1.96 * est_se)), alpha = 0.2) +
  labs(x = "Mean Surface Temperature",
       y = "Predicted Overlap (on response scale)",
       title = "Marginal Effect of Surface Temp") +
  theme_minimal()































