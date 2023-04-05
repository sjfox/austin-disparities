#### Make figures for manuscript
library(tidyverse)
library(lubridate)
library(cowplot)
library(sf)
theme_set(theme_cowplot())

source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')
source('R/load-results.R')


# Seroprevalence conversion to Austin -------------------------------------
tx_pop <- read_csv('processed-data/population-data/tx_pop_age.csv')
read_csv('raw-data/Nationwide_Commercial_Laboratory_Seroprevalence_Survey.csv') %>% 
  filter(Site == 'TX') %>% 
  rename(pop = `Catchment population`,
         date_range = `Date Range of Specimen Collection`) %>% 
  separate(date_range, into = c('date_begin', 'date_end'), sep = ' - ') %>% 
  mutate(date_begin = trimws(date_begin),
         date_end = trimws(date_end)
  ) %>% 
  separate(date_begin, into = c('begin_day', 'begin_year'), sep = ', ') %>% 
  separate(date_end, into = c('end_day', 'end_year'), sep = ', ') %>% 
  mutate(begin_year = ifelse(is.na(begin_year), end_year, begin_year),
         begin_date = paste0(begin_day, ', ', begin_year)) %>% 
  mutate(begin_date = mdy(begin_date)) %>% 
  mutate(infection_date = begin_date - days(7)) %>% 
  select(state = Site, 
         begin_date,
         round = Round,
         `0-17_mean` = `Rate (%) [0-17 Years Prevalence]`,
         `0-17_lower` = `Lower CI [0-17 Years Prevalence]`,
         `0-17_upper` = `Upper CI [0-17 Years Prevalence]`,
         `18-49_mean` = `Rate (%) [18-49 Years Prevalence]`,
         `18-49_lower` = `Lower CI [18-49 Years Prevalence]`,
         `18-49_upper` = `Upper CI [18-49 Years Prevalence]`,
         `50-64_mean` = `Rate (%) [50-64 Years Prevalence]`,
         `50-64_lower` = `Lower CI [50-64 Years Prevalence]`,
         `50-64_upper` = `Upper CI [50-64 Years Prevalence]`,
         `65+_mean` = `Rate (%) [65+ Years Prevalence]`,
         `65+_lower` = `Lower CI [65+ Years Prevalence]`,
         `65+_upper` = `Upper CI [65+ Years Prevalence]`) %>% 
  gather(key, value, -state, -begin_date, -round) %>% 
  separate(key, into = c('age_group', 'metric'), sep = '_') %>% 
  spread(metric, value) %>% 
  mutate(mean = ifelse(mean == 777, NA, mean)) %>% 
  left_join(tx_pop, by = 'age_group') %>% 
  mutate_at(vars(lower:upper), ~./100*pop) -> tx_seroprevalence


age_zip_inf_df %>% 
  unnest(inf_timing) %>% 
  # mutate(zip_age = paste0(zip,'_',age_group_paper)) %>% 
  select(id, zip, age_group_paper, date, new_infections) %>% 
  spread(zip, new_infections) %>% 
  mutate(new_inf = rowSums(.[c(-1,-2,-3)], na.rm = T)) %>% 
  select(id, age_group_paper, date, new_inf) %>% 
  group_by(id, age_group_paper) %>% 
  mutate(cum_inf = cumsum(new_inf)) -> age_inf_cum_ts

age_inf_cum_ts %>% 
  ungroup() %>% 
  select(age_group_paper, id, date, cum_inf) %>% 
  spread(id, cum_inf) %>% 
  group_by(age_group_paper) %>% 
  fill(names(.)) %>%
  mutate_at(.vars = vars(-date,-age_group_paper), ~ifelse(is.na(.), 0, .)) %>% 
  gather(id, cum_inf, -date,-age_group_paper) %>% 
  group_by(age_group_paper, date) %>% 
  summarize(meaninf = mean(cum_inf),
            lowerinf = quantile(cum_inf, probs = 0.025, type = 3),
            upperinf = quantile(cum_inf, probs = 0.975, type = 3)) -> age_inf_cum_ts



age_inf_cum_ts %>% 
  left_join(age_cov_df, by = 'age_group_paper') %>% 
  mutate(meaninf = meaninf/population,
         lowerinf = lowerinf/population,
         upperinf = upperinf/population) %>% 
  left_join(tx_seroprevalence %>% 
              mutate(age_group_paper = age_group), by = c('age_group_paper', 'date' = 'begin_date')) %>% 
  filter(date>=first_infection_date, date < final_infection_date) %>% 
  ggplot(aes(date, meaninf)) + 
  theme(strip.background = element_rect(fill=NA)) +
  geom_ribbon(aes(ymin = lowerinf, ymax=upperinf),alpha = .1) +
  geom_line() +
  geom_point(aes(y = mean/pop),color='darkred') +
  geom_errorbar(aes(ymin = lower/pop, ymax = upper/pop), color='darkred') +
  facet_wrap(~age_group_paper) +
  background_grid(major = 'xy', minor = 'y')+
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y') +
  labs(x =NULL, y='Cumulative infections (%)') ->age_sero_comparison
age_sero_comparison
save_plot('figs/age_sero_comparison.png',
          plot=age_sero_comparison, 
          bg = 'white', 
          base_height = 6, 
          base_asp = 1.8)

age_inf_cum_ts %>% 
  left_join(age_cov_df, by = 'age_group_paper') %>% 
  mutate(meaninf = meaninf/population,
         lowerinf = lowerinf/population,
         upperinf = upperinf/population) %>% 
  left_join(tx_seroprevalence %>% 
              mutate(age_group_paper = age_group), by = c('age_group_paper', 'date' = 'begin_date')) %>% 
  filter(date>=first_infection_date, date < final_infection_date) %>% 
  mutate(ratio = meaninf/(mean/pop)) %>% 
  filter(!is.na(ratio)) %>% 
  ggplot(aes(date, ratio, color = age_group_paper)) + 
    geom_point() +
    geom_smooth(se=F) +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_brewer(type = 'qual', palette = 2) +
    labs(x = NULL, y = 'Mean relative infection risk', color = NULL) +
    background_grid(major = 'xy')+
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y') ->age_tx_ratios
age_tx_ratios
save_plot('figs/age_tx_ratios.png',
          plot=age_tx_ratios, 
          bg = 'white', 
          base_height = 5, 
          base_asp = 1.6)


## Estimate final ratios
get_norm_sd <- function(mean, lwr, upr) {
  ## Conservatively takes the maximum standard deviation
  ## 1.96 is critical value to get standard deviation
  ## exponentiate to get variance
  (pmax((mean-lwr), (upr-mean))/1.96)
}

get_tx_inf_samps <- function(inf_mean, inf_sd, nsamps = 1000){
  round(rnorm(nsamps, mean = inf_mean, sd = inf_sd))
}

age_zip_inf_df %>% 
  select(zip, age_group_paper, est_infections) %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(id = seq_along(est_infections)) %>% 
  group_by(age_group_paper, id) %>% 
  summarize(travis_infections = sum(est_infections), .groups = 'drop') %>% 
  left_join(age_cov_df, by = 'age_group_paper') %>% 
  mutate(travis_inf_rate = travis_infections/population) %>% 
  select(age_group_paper, id, travis_inf_rate) -> trav_age_inf


tx_seroprevalence %>% 
  filter(begin_date == '2021-05-27') %>% 
  mutate(final_sd = get_norm_sd(mean, lower, upper)) %>% 
  mutate(tx_age_inf_samps = map2(mean, final_sd, get_tx_inf_samps, nsamps = max(trav_age_inf$id))) %>% 
  unnest(tx_age_inf_samps) %>% 
  mutate(tx_age_inf_rate = tx_age_inf_samps/pop) %>% 
  select(age_group, tx_age_inf_rate) %>% 
  group_by(age_group) %>% 
  mutate(id = seq_along(tx_age_inf_rate)) %>% 
  ungroup() -> tx_age_inf

trav_age_inf %>% 
  left_join(tx_age_inf, by = c('age_group_paper' = 'age_group', 'id')) %>% 
  mutate(ratio = 1 - travis_inf_rate/tx_age_inf_rate) %>% 
  group_by(age_group_paper) %>% 
  summarize(ratiomean = mean(ratio),
            ratiolo = quantile(ratio, probs = 0.025),
            ratiohi = quantile(ratio, probs = 0.975))
