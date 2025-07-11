year_filt <- c(seq(1970, 2016, by = 1),2018,2019,2021,2022,2023, 2024)

bt_temps <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds")) %>%
  select(-c(SST_seasonal)) %>%
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
  filter(!time %in% c(2024.50, 2024.75)) %>%
  select(longitude, latitude, season, year, BT_seasonal, time, Date) %>% 
  arrange(season, year, time, Date, longitude, latitude) %>%
  ungroup() %>%
  mutate(id = 1:n())


df_overlap %>% 
  left_join(bt_temps) %>%
  mutate(overlap = predator_est + prey_est) %>%
  ggplot(aes(x = BT_seasonal, y = overlap))+
  geom_point()

df_overlap %>% 
  left_join(bt_temps) %>%
  mutate(overlap = predator_est + prey_est) %>%
  ggplot(aes(x = time, y = overlap))+
  geom_point()

df_overlap %>% 
  mutate(space_id = paste(latitude, longitude), 
         overlap = predator_est + prey_est, 
         dominance = predator_est - prey_est) %>%
  filter(space_id %in% sample(space_id, 10)) %>%
  ggplot(aes(x = time, y = dominance))+
  geom_line(aes(group = space_id))+
  facet_wrap(~space_id)

library(mgcv)

# Step 1: Create overlap metric (log predator + log prey)
df_model <- df_overlap %>%
  left_join(bt_temps) %>%
  mutate(overlap = predator_est + prey_est, 
         decade = floor(time/10)*10) %>%
  filter(!is.na(BT_seasonal))

# Step 2: Fit GAM with smooth for temperature
gam_fit <- gam(overlap ~ s(BT_seasonal,by = as.factor(decade), k = 4), data = df_model)

# Step 3: Summarize model
summary(gam_fit)

plot(ggeffects::ggpredict(gam_fit, terms = ~BT_seasonal+decade), add.data = T)

# Step 4: Plot effect of temperature on overlap
plot_df <- df_model %>%
  mutate(predicted = predict(gam_fit, newdata = .))

ggplot(plot_df, aes(x = BT_seasonal, y = overlap)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  labs(
    x = "Bottom Temperature (seasonal)",
    y = "Predator–Prey Overlap (log scale)",
    title = "Effect of Temperature on Predator–Prey Overlap"
  ) +
  theme_minimal()

