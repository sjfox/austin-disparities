#### Make figures for manuscript
library(tidyverse)
library(lubridate)
library(cowplot)
library(sf)
theme_set(theme_cowplot())

source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')
source('R/load-results.R')


# Seroprevalence for Texas -------------------------------------

read_csv('raw-data/Nationwide_Commercial_Laboratory_Seroprevalence_Survey.csv', guess_max = 2000) %>% 
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
         overall_sero = `Rate (%) [Cumulative Prevalence]`,
         lb_sero = `Lower CI [Cumulative Prevalence]`,
         ub_sero = `Upper CI [Cumulative Prevalence]`)  -> tx_seroprevalence



age_zip_df %>% group_by(age_group_paper) %>% summarize(sum(n_admits))

all_inf_ts %>% 
  filter(date == max(date)) %>% 
  select(id, cum_inf = cum_inf) %>% 
  mutate(rr = samp_conditional_reporting_rate(detected_cases = travis_total_cases,
                                              estimated_infections = cum_inf,
                                              num_samps = length(cum_inf))) %>% 
  ungroup() %>% 
  summarize(cuminf_mean = mean(cum_inf)/travis_pop,
            cuminf_lo = quantile(cum_inf, probs = 0.025)/travis_pop,
            cuminf_hi = quantile(cum_inf, probs = 0.975)/travis_pop,
            rr_mean = mean(rr),
            rr_lo = quantile(rr, probs = 0.025),
            rr_hi = quantile(rr, probs = 0.975)
  )

# Fig 1 time-series data and comparison sero and inf estimate ---------------------------
age_zip_inf_df %>% 
  filter(n_admits != 0) %>%
  unnest(admission_dates) %>% 
  count(admission_dates) %>% 
  ggplot(aes(admission_dates, n/travis_pop*1000000)) + 
  geom_col(width = 1) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = 'Hospital admissions (per 1M)') +
  background_grid(major = 'y', minor = 'y') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y', 
               limits = c(ymd(first_infection_date), ymd(final_infection_date))) -> f1_hosp_admit_plot

all_inf_ts %>% 
  ungroup() %>% 
  select(id, date, cum_inf) %>% 
  spread(id, cum_inf) %>% 
  fill(names(.)) %>%
  mutate_at(.vars = vars(-date), ~ifelse(is.na(.), 0, .)) %>% 
  gather(id, cum_inf, -date) %>% 
  group_by(date) %>% 
  summarize(meaninf = mean(cum_inf),
            lowerinf = quantile(cum_inf, probs = 0.025, type = 3),
            upperinf = quantile(cum_inf, probs = 0.975, type = 3)) %>% 
  filter(date>=first_infection_date, date < final_infection_date) %>% 
  mutate_at(.vars = vars(meaninf:upperinf), ~./travis_pop) %>% 
  ggplot(aes(date, meaninf)) +
  geom_ribbon(aes(ymin = lowerinf, ymax = upperinf), color = NA, alpha = .3) +
  geom_line() +
  geom_point(data = tx_seroprevalence, aes(x=begin_date, y = overall_sero/100), inherit.aes = FALSE, color = 'darkred') +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  geom_errorbar(data = tx_seroprevalence, aes(x=begin_date, ymin = lb_sero/100, ymax = ub_sero/100), inherit.aes = FALSE, color = 'darkred') +
  labs(x = NULL, y = 'Cumulative infections (%)') +
  background_grid(major = 'y', minor = 'y') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y',
               limits = c(ymd(first_infection_date), ymd(final_infection_date))) -> f1_cum_inf_comparison

plot_grid(f1_hosp_admit_plot, f1_cum_inf_comparison, nrow = 2, labels = c('A', 'B'), align = 'v') -> f1_combined
save_plot('figs/fig1_combined.png', f1_combined, base_height = 6.5, base_asp = 1.3,bg = "white")  
save_plot('figs/fig1.tiff', 
          f1_combined, 
          base_height = 6.5, 
          base_asp = 1.3, 
          bg = "white")  


all_inf_ts %>% 
  ungroup() %>% 
  select(id, date, cum_inf) %>% 
  spread(id, cum_inf) %>% 
  fill(names(.)) %>%
  mutate_at(.vars = vars(-date), ~ifelse(is.na(.), 0, .)) %>% 
  gather(id, cum_inf, -date) %>% 
  group_by(date) %>% 
  summarize(meaninf = mean(cum_inf),
            lowerinf = quantile(cum_inf, probs = 0.025, type = 3),
            upperinf = quantile(cum_inf, probs = 0.975, type = 3)) %>% 
  filter(date == '2020-09-23') %>% 
  mutate_at(.vars = vars(meaninf:upperinf), ~./travis_pop)

tx_seroprevalence %>% 
  filter(begin_date == '2020-09-23')

save.image(file='processed-data/fig1-image.rda')
