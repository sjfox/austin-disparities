#### Make figures for manuscript
library(tidyverse)
library(lubridate)
library(cowplot)
library(sf)
theme_set(theme_cowplot())

source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')
source('R/load-results.R')


## Hospitalization rate
age_zip_inf_df %>% 
  group_by(age_group_paper) %>% 
  summarize(n_admits = sum(n_admits)) %>% 
  left_join(age_cov_df) %>% 
  mutate(hosp_rate = n_admits/population*100000)-> age_hosp_rate

## Compare 65+ to everyone else
age_65plus_pop <- age_cov_df %>% 
  mutate(age_65plus = ifelse(age_group_paper == '65+', '65+', '<65')) %>%
  group_by(age_65plus) %>% 
  summarize(population = sum(population),
            cases = sum(cases))

age_zip_inf_df %>% 
  mutate(age_65plus = ifelse(age_group_paper == '65+', '65+', '<65')) %>% 
  group_by(age_65plus) %>% 
  summarize(n_admits = sum(n_admits)) %>% 
  left_join(age_65plus_pop) %>% 
  mutate(hosp_rate = n_admits/population*100000)

age_hosp_rate %>% 
  ggplot(aes(age_group_paper, hosp_rate, fill = age_group_paper)) + 
  geom_col() + 
  geom_hline(yintercept = age_hosp_rate %>% 
               summarize(avg_rate = sum(n_admits)/sum(population)*100000) %>% 
               pull(avg_rate),
             lty = 2) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = 'Hospitalizations per 100k') +
  scale_color_brewer(type = 'qual', palette = 2) + 
  scale_fill_brewer(type = 'qual', palette = 2) +
  theme(legend.position = 'none') + 
  background_grid(major = 'y', minor = 'y') -> age_hosp_rate_plot
age_hosp_rate_plot

## Reported case rate
age_cov_df %>%
  mutate(case_rate = cases/population*100000)

age_cov_df %>%
  mutate(case_rate = cases/population*100000) %>% 
  ggplot(aes(age_group_paper, case_rate, fill = age_group_paper)) + 
  geom_col() + 
  geom_hline(yintercept = age_cov_df %>% 
               summarize(avg_rate = sum(cases)/sum(population)*100000) %>% 
               pull(avg_rate),
             lty = 2) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = 'Reported cases per 100k') +
  scale_color_brewer(type = 'qual', palette = 2) + 
  scale_fill_brewer(type = 'qual', palette = 2) + 
  theme(legend.position = 'none') + 
  background_grid(major = 'y', minor = 'y') -> age_case_rate_plot
age_case_rate_plot

## Overall infection plot
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(age_group_paper, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  summarize(meaninf = mean(inf),
            lowerinf = quantile(inf, probs = 0.025),
            upperinf = quantile(inf, probs = 0.975)) %>% 
  left_join(age_cov_df) %>% 
  mutate_at(vars(meaninf:upperinf), ~./population) ->age_infection_rate
age_infection_rate


## 65+ vs everyone else infection rate
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  mutate(age_65plus = ifelse(age_group_paper == '65+', '65+', '<65')) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(age_65plus, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  summarize(meaninf = mean(inf),
            lowerinf = quantile(inf, probs = 0.025),
            upperinf = quantile(inf, probs = 0.975)) %>% 
  left_join(age_65plus_pop) %>% 
  mutate_at(vars(meaninf:upperinf), ~./population) -> age65plus_inf_rate

age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  summarize(meaninf = mean(inf),
            lowerinf = quantile(inf, probs = 0.025),
            upperinf = quantile(inf, probs = 0.975)) %>% 
  mutate_at(vars(meaninf:upperinf), ~./travis_pop) ->overall_infection_rate


age_infection_rate %>% 
  ggplot(aes(age_group_paper, meaninf, color = age_group_paper)) + 
  geom_hline(yintercept = overall_infection_rate$meaninf, lty = 2) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = lowerinf, ymax = upperinf), size = 0.75) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),limits = c(0, 0.35)) +
  labs(x = NULL, y = 'Cumulative infections (%)') +
  scale_color_brewer(type = 'qual', palette = 2) + 
  scale_fill_brewer(type = 'qual', palette = 2) + 
  theme(legend.position = 'none') + 
  background_grid(major = 'y', minor = 'y') -> age_cum_inf_plot
age_cum_inf_plot

## Reporting rate plot
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(age_group_paper, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  left_join(age_cov_df) %>% 
  mutate(rr = samp_conditional_reporting_rate(detected_cases = cases, 
                                              estimated_infections = inf, 
                                              num_samps = n())) %>% 
  summarize(meanrr = mean(rr),
            lowerrr = quantile(rr, probs = 0.025),
            upperrr = quantile(rr, probs = 0.975)) -> age_rr
age_rr

##65+ vs everyone else
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  mutate(age_65plus = ifelse(age_group_paper == '65+', '65+', '<65')) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(age_65plus, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  left_join(age_65plus_pop) %>% 
  mutate(rr = samp_conditional_reporting_rate(detected_cases = cases, 
                                              estimated_infections = inf, 
                                              num_samps = n())) %>% 
  summarize(meanrr = mean(rr),
            lowerrr = quantile(rr, probs = 0.025),
            upperrr = quantile(rr, probs = 0.975)) -> age65plus_rr
age65plus_rr

age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  mutate(rr = samp_conditional_reporting_rate(detected_cases = sum(age_cov_df$cases), 
                                              estimated_infections = inf, 
                                              num_samps = n())) %>% 
  summarize(meanrr = mean(rr),
            lowerrr = quantile(rr, probs = 0.025),
            upperrr = quantile(rr, probs = 0.975)) -> overall_rr

age_rr %>%
  ggplot(aes(age_group_paper, meanrr, color = age_group_paper)) + 
  geom_hline(yintercept = overall_rr$meanrr, lty = 2) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = lowerrr, ymax = upperrr), size = 0.75) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),limits = c(0, 0.6)) +
  labs(x = NULL, y = 'Reporting rate (%)') +
  scale_color_brewer(type = 'qual', palette = 2) + 
  scale_fill_brewer(type = 'qual', palette = 2) + 
  theme(legend.position = 'none') + 
  background_grid(major = 'y', minor = 'y') -> age_rr_plot
age_rr_plot

## First panel of age comparison
plot_grid(age_hosp_rate_plot, age_case_rate_plot, age_cum_inf_plot, age_rr_plot, 
          nrow = 2, align = 'hv', labels = 'AUTO') -> cum_age_metrics_plot

## Time-series of age infections
age_zip_inf_df %>% 
  unnest(inf_timing) %>% 
  filter(date>=first_infection_date, date < final_infection_date) %>% 
  group_by(age_group_paper, id, date) %>% 
  summarize(new_inf = sum(new_infections)) %>% 
  padr::pad(group = c('age_group_paper', 'id'), break_above = 3000000) %>% 
  mutate(new_inf = ifelse(is.na(new_inf), 0, new_inf)) %>% 
  group_by(age_group_paper, date) %>% 
  summarize(meaninf = mean(new_inf),
            lowerinf = quantile(new_inf, probs = 0.025),
            upperinf = quantile(new_inf, probs = 0.975)) -> age_inf_ts

age_inf_ts %>% 
  left_join(age_cov_df) %>% 
  mutate_at(vars(meaninf:upperinf), ~./population*100000) %>% 
  filter(date>'2020-02-01', date < '2021-05-01') %>% 
  ggplot(aes(date, meaninf, color = age_group_paper, fill = age_group_paper)) + 
  geom_ribbon(aes(ymin = lowerinf, ymax = upperinf), alpha = .2, color = NA) +
  geom_line(size = 1) +
  scale_y_continuous() +
  labs(x = NULL, y = 'Infections per 100k', color = NULL, fill = NULL) +
  scale_color_brewer(type = 'qual', palette = 2) + 
  scale_fill_brewer(type = 'qual', palette = 2) + 
  # theme_cowplot(font_size = 21)   +
  background_grid(major = 'xy', minor = 'y') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y') -> age_inf_ts_plot
age_inf_ts_plot


age_zip_inf_df %>% 
  unnest(inf_timing) %>% 
  filter(date >= '2020-03-01') %>% 
  mutate(wave = case_when(date>=ymd('2020-03-01') & date < ymd('2020-05-01') ~ 'Spring',
                          date>=ymd('2020-06-01') & date < ymd('2020-08-01') ~ 'Summer',
                          date>=ymd('2020-12-01') & date < ymd('2021-02-01') ~ 'Winter',
                          T ~ 'neither')) %>% 
  group_by(age_group_paper, id, wave) %>% 
  summarize(new_inf = sum(new_infections)) %>% 
  group_by(wave, id) %>% 
  mutate(frac_inf = new_inf/sum(new_inf)) %>% 
  ungroup() %>% 
  group_by(age_group_paper, wave) %>% 
  summarize(meaninf = mean(new_inf),
            lowerinf = quantile(new_inf, probs = 0.025),
            upperinf = quantile(new_inf, probs = 0.975),
            meanfracinf = mean(frac_inf),
            lowerfracinf = quantile(frac_inf, probs = 0.025),
            upperfracinf = quantile(frac_inf, probs = 0.975)) -> age_inf_wave
age_inf_wave %>% 
  filter(wave != 'neither') %>% 
  left_join(age_cov_df %>% mutate(popfrac = population/sum(population)) %>% select(age_group_paper,popfrac)) %>% 
  ggplot(aes(wave, meanfracinf, fill = age_group_paper)) +
  geom_col(alpha = .5) +
  geom_errorbar(aes(ymin = lowerfracinf, ymax = upperfracinf), width = .5) +
  facet_wrap(~age_group_paper,  nrow = 1) +
  theme(strip.background = element_rect(fill=NA), 
        axis.text.x = element_text(angle = 45, hjust =1,vjust=1)) +
  labs(x = NULL, y = 'Fraction of all infections', fill = NULL, color = NULL) +
  background_grid(major = 'y', minor = 'y') +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), expand = c(0,0), limits = c(0,0.75)) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  geom_hline(aes(yintercept = popfrac, color = age_group_paper), lty = 1) +
  scale_color_brewer(type = 'qual', palette = 2) -> age_wave_summary_plot
age_wave_summary_plot

# ## Weekly infection breakdown for each age group
age_zip_inf_df %>%
  unnest(inf_timing) %>%
  filter(date >= '2020-03-01') %>%
  mutate(week = floor_date(date, unit = 'week')) %>%
  group_by(age_group_paper, id, week) %>%
  summarize(new_inf = sum(new_infections)) %>%
  group_by(week, id) %>%
  mutate(frac_inf = new_inf/sum(new_inf)) %>%
  ungroup() %>%
  group_by(age_group_paper, week) %>%
  summarize(meaninf = mean(new_inf),
            lowerinf = quantile(new_inf, probs = 0.025),
            upperinf = quantile(new_inf, probs = 0.975),
            meanfracinf = mean(frac_inf),
            lowerfracinf = quantile(frac_inf, probs = 0.025),
            upperfracinf = quantile(frac_inf, probs = 0.975)) -> age_inf_week
# 
# age_inf_week %>% 
#   filter(week>=first_infection_date, week < final_infection_date) %>% 
#   left_join(age_cov_df %>% mutate(fracpop = population/sum(population))) %>% 
#   ggplot(aes(week, meanfracinf, fill = age_group_paper)) +
#   geom_col(width = 7) +
#   labs(x = NULL, y = 'Proportion of infections', color = NULL, fill = NULL) +
#   scale_color_brewer(type = 'qual', palette = 2) + 
#   scale_fill_brewer(type = 'qual', palette = 2) + 
#   scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0,1.01)) +
#   background_grid(major = 'xy') +
#   scale_x_date(date_breaks = '2 month', date_labels = '%b-%y',
#                limits = c(ymd(first_infection_date), ymd(final_infection_date))) -> week_age_plot
# week_age_plot
# plot_grid(age_inf_ts_plot + theme(legend.position = 'none'), 
#           week_age_plot + theme(legend.position = 'none'), 
#           align = 'v', nrow = 2, labels = c('E', 'F')) %>% 
#   plot_grid(get_legend(week_age_plot), nrow = 1, rel_widths = c(1,.12)) -> age_ts_plots


# plot_grid(age_inf_ts_plot + theme(legend.position = 'none'),
#           age_wave_summary_plot + theme(legend.position = 'none'),
#           align = 'v', nrow = 2, labels = c('E', 'F'), axis = 'trbl') %>%
#   plot_grid(get_legend(age_wave_summary_plot), nrow = 1, rel_widths = c(1,.12)) -> age_ts_plots

# plot_grid(cum_age_metrics_plot, age_ts_plots, labels = c(NA, NA), nrow = 1, rel_widths = c(1,1)) -> fig2
# save_plot('figs/fig2.png', fig2, base_height = 6.5, base_asp = 2.2,bg = "white")


# Alternate version of figure 2 -------------------------------------------
plot_grid(age_hosp_rate_plot, age_case_rate_plot, age_cum_inf_plot, age_rr_plot, 
          nrow = 1, align = 'hv', labels = 'AUTO') -> cum_age_metrics_plot_alt
plot_grid(age_inf_ts_plot + theme(legend.position = 'none'),
          age_wave_summary_plot + theme(legend.position = 'none'),
          align = 'v', nrow = 1, labels = c('E', 'F'), axis = 'trbl')  -> age_ts_plots_alt
plot_grid(cum_age_metrics_plot_alt, age_ts_plots_alt, labels = c(NA, NA), nrow = 2, rel_heights = c(1,1)) %>%
  plot_grid(get_legend(age_wave_summary_plot), nrow = 1, rel_widths = c(1,.08)) -> fig2_alt
save_plot('figs/fig2.png', fig2_alt, base_height = 6.5, base_asp = 2.2,bg = "white")
save_plot('figs/fig2.tiff', fig2_alt, base_height = 6.5, base_asp = 2.2,bg = "white")


## Children wave stats
age_inf_wave %>% 
  filter(wave != 'neither') %>% 
  left_join(age_cov_df %>% 
              mutate(popfrac = population/sum(population)) %>% 
              select(age_group_paper,popfrac))

save.image(file='processed-data/fig2-image.rda')
