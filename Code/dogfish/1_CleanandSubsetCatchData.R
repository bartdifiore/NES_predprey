#----------------------------------
## Libraries
#----------------------------------

library(tidyverse)
library(sf)
library(gmRi)
source("Code/enhance_r_funcs.R")

pred_subset <- "Squalus acanthias"

#----------------------------------
## Load & filter predators
#----------------------------------
ma <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  filter(scientific_name %in% pred_subset) %>%
  mutate(date = lubridate::date(date))

dfo <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds")  %>%
  filter(scientific_name %in% pred_subset)

nefsc <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>% 
  mutate(scientific_name = str_to_sentence(scientific_name)) %>%
  filter(scientific_name %in% pred_subset)

me <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  filter(scientific_name %in% pred_subset) %>% 
  mutate(date = lubridate::date(date))

pred_df <- rbind(ma, dfo, nefsc, me) %>%
  mutate(season = str_to_sentence(season))

write_rds(pred_df, "Data/Derived/combined_and_filtered_predators.rds", compress = "gz")


#----------------------------------
## Load and filter prey
#----------------------------------
ma_mac <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  filter(scientific_name == "Scomber scombrus") |>
  mutate(
    date = lubridate::date(date)
  )

dfo_mac <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds") %>%
  filter(scientific_name == "Scomber scombrus") %>%
  mutate(
    date = lubridate::date(date)
  )

nefsc_mac <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>%
  mutate(scientific_name = str_to_sentence(scientific_name)) %>%
  filter(scientific_name == "Scomber scombrus") |>
  mutate(
    date = lubridate::date(date)
  )

me_mac <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  filter(scientific_name == "Scomber scombrus") |>
  mutate(
    date = lubridate::date(date)
  )

mac_df <- rbind(ma_mac, dfo_mac, nefsc_mac, me_mac) %>%
  mutate(season = str_to_sentence(season))

# Some quick exploration
# Tows
mac_df %>% 
  ungroup() %>%
  distinct(trawl_id, survey) %>%
  group_by(survey) %>%
  summarize(n = n())

# Total observed
mac_df %>% 
  group_by(survey) %>%
  summarize(total_observed = sum(number_at_length, na.rm = T))

mac_df %>% 
  group_by(survey, year) %>% 
  summarize(total_observed = sum(number_at_length)) %>% 
  ggplot(aes(x = year, y = total_observed))+
  geom_line(aes(color = survey))+
  theme_classic()

mac_df %>% 
  group_by(survey, year) %>% 
  summarize(total_observed = sum(number_at_length)) %>% 
  ggplot(aes(x = year, y = total_observed))+
  geom_line(aes(color = survey))+
  facet_wrap(~survey + life_class, scales = "free", ncol = 2)+
  theme_classic()


ggplot(mac_df, aes(x = weight_at_length))+
  geom_histogram()+
  facet_wrap(~survey, scales = "free", ncol = 2) 


write_rds(mac_df, "Data/Derived/mackerel.rds", compress = "gz")



#----------------------------------
## Initial visualizations
#----------------------------------

mac_df %>% 
  bind_rows(pred_df) %>%
  group_by(scientific_name, year, season) %>% 
  summarize(total = sum(number_at_length, na.rm = T)) %>% 
  group_by(scientific_name) %>%
  mutate(total_scaled = total/max(total, na.rm = T)) %>%
  ggplot(aes(x = year, y = total_scaled))+
  geom_line(aes(color = scientific_name))+
  facet_wrap(~season)













