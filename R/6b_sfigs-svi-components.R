#### Make figures for manuscript
library(tidyverse)
library(lubridate)
library(cowplot)
library(sf)
library(lme4)
library(effects)
library(sjPlot)
theme_set(theme_cowplot())

source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')
source('R/load-results.R')

plot_grid <- cowplot::plot_grid
save_plot <- cowplot::save_plot

# fig 4: Correlations SVI and infections/reporting rate ----------------------

age_zip_inf_df %>% 
  select(zip, age_group_paper, est_infections) %>% 
  unnest(est_infections) %>% 
  group_by(zip,age_group_paper) %>% 
  mutate(id = seq_along(est_infections)) %>% 
  group_by(zip,id) %>% 
  summarize(infections = sum(est_infections)) %>% 
  ungroup() %>% 
  left_join(zip_cov_df, by = 'zip') %>% 
  mutate(cases = ifelse(infections < cases, infections, cases)) %>% ## Makes it so maximum reporting rate is 100%
  mutate(rr = samp_conditional_reporting_rate(cases, infections, n()),
         inf_rate = infections/population) %>% 
  select(id, zip, inf_rate, cases, rr, svi, infections, population)-> zip_corr_long_df
  
zip_corr_long_df %>% 
  group_by(zip) %>% 
  summarize(svi = unique(svi),
            inf_mean = mean(inf_rate),
            inf_lo = quantile(inf_rate, probs = 0.025),
            inf_hi = quantile(inf_rate, probs = 0.975),
            rr_mean = mean(rr),
            rr_lo = quantile(rr, probs = 0.025),
            rr_hi = quantile(rr, probs = 0.975)
  ) %>% 
  left_join(get_SVI_per_ZIP() %>% 
            mutate(ZIP = as.character(ZIP)) %>% 
            select(-SVI), by = c('zip' = 'ZIP')) -> zip_corr_short_df
zip_corr_short_df %>% 
  select(-RPL_THEME1, -RPL_THEME2, -RPL_THEME3, -RPL_THEME4, -ZIP_COVERED, -zip,-svi) ->zip_corr_short_df

clean_keys <- tibble(key = colnames(zip_corr_short_df)[c(-1,-2, -3,-4,-5,-6)],
                     cleaned_key = c('Below poverty (%)',
                                     'Unemployment (%)',
                                     'Per capita income ($)',
                                     'No high school diploma (%)',
                                     'Age 65+ (%)',
                                     'Age <17 (%)',
                                     'Noninstitutional disability (%)',
                                     'Single parent houshold (%)',
                                     'Minority population (%)',
                                     'Limited english language (%)',
                                     'Housing >10 units (%)',
                                     'Mobile home (%)',
                                     'Crowded housing (%)',
                                     'No vehicle (%)',
                                     'Group quarters (%)'
                                     ))

zip_corr_short_df %>% 
  select(-starts_with('rr')) %>%
  gather(key, value, -inf_mean, -inf_lo, -inf_hi) %>% 
  left_join(clean_keys, by = 'key') %>% 
  ggplot(aes(value, inf_mean)) +
    geom_errorbar(aes(ymin = inf_lo, ymax = inf_hi), alpha=.5) +  
    geom_point() +
    facet_wrap(~cleaned_key, scales = 'free_x', nrow = 5,
               strip.position = 'bottom') +
    background_grid(major = 'xy') +
    scale_y_continuous(labels = scales::percent) +
    theme(panel.spacing.x = unit(6, "mm"),
          strip.background = element_rect(fill=NA),
          strip.placement = "outside")+
    labs(x = NULL, y = 'Cumulative infections (%)') ->inf_rate_svi_components
inf_rate_svi_components
save_plot('figs/sfig-inf_rate_svi_components.png', 
          inf_rate_svi_components, 
          base_height = 12, 
          base_asp = 0.8, bg= 'white')

zip_corr_short_df %>% 
  select(-starts_with('inf')) %>%
  gather(key, value, -rr_mean, -rr_lo, -rr_hi) %>% 
  left_join(clean_keys, by = 'key') %>% 
  ggplot(aes(value, rr_mean)) +
  geom_errorbar(aes(ymin = rr_lo, ymax = rr_hi), alpha=.5) +  
  geom_point() +
  facet_wrap(~cleaned_key, scales = 'free_x', nrow = 5,
             strip.position = 'bottom') +
  background_grid(major = 'xy') +
  scale_y_continuous(labels = scales::percent) +
  theme(panel.spacing.x = unit(6, "mm"),
        strip.background = element_rect(fill=NA),
        strip.placement = "outside")+
  labs(x = NULL, y = 'Reporting rate (%)') ->rep_rate_svi_components
rep_rate_svi_components
save_plot('figs/sfig-rep_rate_svi_components.png', 
          rep_rate_svi_components, 
          base_height = 12, 
          base_asp = 0.8, bg= 'white')

## Get correlations
zip_corr_long_df %>%
  left_join(get_SVI_per_ZIP() %>%
              mutate(ZIP = as.character(ZIP)) %>%
              select(-SVI), by = c('zip' = 'ZIP')) %>%
  select(-cases, -infections, -population,
         -RPL_THEME1, -RPL_THEME2, -RPL_THEME3,
         -RPL_THEME4, -ZIP_COVERED, -zip) %>%
  select(-rr) %>%
  gather(key, value, -inf_rate, -id) %>%
  group_by(key, id) %>%
  summarize(corr = cor(value, inf_rate)) %>%
  group_by(key) %>%
  summarize(cor_mean = mean(corr),
            cor_lo = quantile(corr, probs = 0.025),
            cor_hi = quantile(corr,probs = 0.975)) %>%
  left_join(clean_keys, by = 'key') %>% 
  arrange(cor_mean)

zip_corr_long_df %>%
  left_join(get_SVI_per_ZIP() %>%
              mutate(ZIP = as.character(ZIP)) %>%
              select(-SVI), by = c('zip' = 'ZIP')) %>%
  select(-cases, -infections, -population,
         -RPL_THEME1, -RPL_THEME2, -RPL_THEME3,
         -RPL_THEME4, -ZIP_COVERED, -zip) %>%
  select(-inf_rate) %>%
  gather(key, value, -rr, -id) %>%
  group_by(key, id) %>%
  summarize(corr = cor(value, rr)) %>%
  group_by(key) %>%
  summarize(cor_mean = mean(corr),
            cor_lo = quantile(corr, probs = 0.025),
            cor_hi = quantile(corr,probs = 0.975)) %>%
  left_join(clean_keys, by = 'key') %>% 
  arrange(cor_mean)


