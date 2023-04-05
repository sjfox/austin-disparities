#### Make figures for manuscript
library(sjPlot)
library(rmapzen)
library(merTools)
library(lubridate)
library(cowplot)
library(sf)
library(tidyverse)
library(effects)

theme_set(theme_cowplot())

plot_grid <- cowplot::plot_grid
save_plot <- cowplot::save_plot


source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')
source('R/load-results.R')

# Supplemental figures ----------------------------------------------------
(aph_zip_age_tests %>% summarize(cases = sum(cases)) %>% pull(cases))/travis_total_cases

## Time-series for other datasets
travis_nyt %>% 
  filter(date>=first_infection_date, date < final_infection_date) %>% 
  mutate(new_cases = c(0,diff(cases))/travis_pop*100000) %>% 
  mutate(new_cases = ifelse(new_cases<0,0,new_cases)) %>% 
  ggplot(aes(date, new_cases)) + 
  geom_col(width = 1) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = 'Reported cases (per 100k)') +
  background_grid(major = 'y', minor = 'y') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y', 
               limits = c(ymd(first_infection_date), ymd(final_infection_date))) -> nyt_case_plot
nyt_case_plot

travis_nyt %>% 
  filter(date>=first_infection_date, date < final_infection_date) %>% 
  mutate(new_deaths = c(0,diff(deaths))/travis_pop*1000000) %>% 
  mutate(new_deaths = ifelse(new_deaths<0,0,new_deaths)) %>% 
  ggplot(aes(date, new_deaths)) + 
  geom_col(width = 1) +
  background_grid(major = 'y', minor = 'y') +
  labs(x = NULL, y = 'Reported deaths (per 1M)') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y', 
               limits = c(ymd(first_infection_date), ymd(final_infection_date))) -> nyt_death_plot
nyt_death_plot


plot_grid(nyt_case_plot, nyt_death_plot, ncol = 1, labels = 'AUTO', align = 'v') -> nyt_travis_plot
save_plot('figs/s1_nyt-travis-plot.png', nyt_travis_plot, base_height = 6, base_asp = 1.3,bg = "white")

## Fraction infected by age for each month...
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
age_inf_week %>% 
  filter(week>=first_infection_date, week < (ymd(final_infection_date)-days(21))) %>% 
  left_join(age_cov_df %>% mutate(fracpop = population/sum(population))) %>% 
  ggplot(aes(week, meanfracinf/fracpop)) +
  geom_point() +
  geom_errorbar(aes(ymin = lowerfracinf/fracpop, ymax = upperfracinf/fracpop)) +
  facet_wrap(~age_group_paper) +
  geom_hline(yintercept = 1, lty = 2) +
  theme(strip.background = element_rect(fill=NA), strip.text = element_text(size = 16), axis.line = element_blank())  +
  background_grid(major = 'xy') +
  annotate("segment", x=structure(-Inf, class = "Date"), xend=structure(-Inf, class = "Date"), y=-Inf, yend=Inf, size = 1)+
  annotate("segment", x=structure(-Inf, class = "Date"), xend=structure(Inf, class = "Date"), y=-Inf, yend=-Inf, size = 1) +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y', expand = c(0,0)) +
  labs(x = NULL, y = 'Relative infection rate') -> weekly_age_relative_inf_plot

save_plot('figs/s2_weekly-age-relative-inf.png', weekly_age_relative_inf_plot, base_height = 6, base_asp = 2,bg = "white")

## Compare age trends for reported cases and hospitalizations
age_cov_df %>% 
  unnest(data) %>% 
  mutate(cum_cases = ifelse(cum_cases == 0, NA_integer_, cum_cases)) %>%
  mutate(cum_cases = ifelse(date == '2020-11-01' & age_group_paper == '65+', NA_integer_, cum_cases)) %>% 
  mutate(cum_cases = ifelse(date == '2020-10-21' & age_group_paper == '65+', NA_integer_, cum_cases)) %>% 
  mutate(cum_cases = ifelse(date == '2020-07-04' & age_group_paper == '18-49', NA_integer_, cum_cases)) %>% 
  mutate(cum_cases = ifelse(date == '2020-10-21' & age_group_paper == '50-64', NA_integer_, cum_cases)) %>% 
  mutate(cum_cases = ifelse(date == '2020-11-01' & age_group_paper == '50-64', NA_integer_, cum_cases)) %>% 
  mutate(cum_cases = ifelse(date == '2021-03-24', NA_integer_, cum_cases)) %>%
  group_by(age_group_paper) %>% 
  mutate(new_cases = diff(c(0, cum_cases))) %>% 
  mutate(new_cases = ifelse(new_cases <0, 0, new_cases)) %>% 
  ungroup() %>% 
  filter(date!= '2020-04-10') %>%  ## Remove first date, because data don't start at 0
  group_by(age_group_paper) %>% 
  padr::pad() %>% 
  mutate(cases_7day = slider::slide_dbl(new_cases/population*100000, .f = mean, .before =3, .after = 3, .complete = T, na.rm=T)) %>% 
  ggplot(aes(date, cases_7day, color = age_group_paper)) + 
  geom_line(size = 0.75) +
  background_grid(major = 'xy') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y', expand = c(0,0)) +
  labs(x = NULL, y = 'Reported cases (per 100k)', color = NULL) +
  scale_color_brewer(type = 'qual', palette = 2) -> reported_age_cases_plot
save_plot('figs/s3_reported-cases-age.png', reported_age_cases_plot, base_height = 4, base_asp = 1.8,bg = "white")

age_zip_inf_df %>% 
  filter(n_admits !=0) %>% 
  select(age_group_paper, admission_dates) %>% 
  unnest(admission_dates) %>% 
  group_by(admission_dates, age_group_paper) %>% 
  summarize(admits = n()) %>% 
  left_join(age_cov_df) %>% 
  mutate(admits = admits/population*1000000) %>% 
  group_by(age_group_paper) %>% 
  padr::pad() %>% 
  ggplot(aes(admission_dates, admits)) + 
  facet_wrap(~age_group_paper)+
  geom_col(width = 1) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = 'Hospital admissions (per 1M)') +
  background_grid(major = 'xy', minor = 'y') +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%y', 
               limits = c(ymd(first_infection_date), ymd(final_infection_date))) +
  annotate("segment", x=structure(-Inf, class = "Date"), xend=structure(-Inf, class = "Date"), y=-Inf, yend=Inf, size = 1)+
  annotate("segment", x=structure(-Inf, class = "Date"), xend=structure(Inf, class = "Date"), y=-Inf, yend=-Inf, size = 1) +
  theme(strip.background= element_rect(fill=NA), 
        axis.line = element_blank(), 
        strip.text = element_text(size = 16)) -> age_hosp_admit_ts_plot
save_plot('figs/s4-age-hosp-counts-ts.png', age_hosp_admit_ts_plot, base_height = 6, base_asp = 1.8,bg = "white")

## Age IHR by zip code
load('processed-data/fig3-image.rda')
age_zip_df %>% 
  filter(age_group_paper == '0-17') %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = ihrmean, 
                plot_title = '0-17 IHR', 
                scales::label_percent(accuracy = .01)) -> zip_ihr_017
age_zip_df %>% 
  filter(age_group_paper == '18-49') %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = ihrmean, 
                plot_title = '18-49 IHR', 
                scales::label_percent(accuracy = .1)) -> zip_ihr_1849
age_zip_df %>% 
  filter(age_group_paper == '50-64') %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = ihrmean, 
                plot_title = '50-64 IHR', 
                scales::label_percent(accuracy = .1)) -> zip_ihr_5064
age_zip_df %>% 
  filter(age_group_paper == '65+') %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = ihrmean, 
                plot_title = '65+ IHR', 
                scales::label_percent(accuracy = .1)) -> zip_ihr_65


plot_grid(zip_ihr_017, zip_ihr_1849,zip_ihr_5064, zip_ihr_65, nrow = 2, align = 'hv') -> age_zip_ihr_map
save_plot('figs/s5_age-zip-ihr-map.png', age_zip_ihr_map, base_height = 5, base_asp = 1.7,bg = "white")

## Age-specific cumulative infection plots
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(population, zip, age_group_paper) %>% 
  summarize(meaninf = mean(est_infections),
            lowerinf = quantile(est_infections, probs = 0.025),
            upperinf = quantile(est_infections, probs = 0.975)) %>% 
  mutate_at(vars(meaninf:upperinf), ~./population) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) -> age_zip_inf

age_zip_inf %>% 
  filter(age_group_paper == '0-17') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = meaninf, 
                plot_title = '0-17 Infection rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_cuminf_017
age_zip_inf %>% 
  filter(age_group_paper == '18-49') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = meaninf, 
                plot_title = '18-49 Infection rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_cuminf_1849
age_zip_inf %>% 
  filter(age_group_paper == '50-64') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = meaninf, 
                plot_title = '50-64 Infection rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_cuminf_5064
age_zip_inf %>% 
  filter(age_group_paper == '65+') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = meaninf, 
                plot_title = '65+ Infection rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_cuminf_65


plot_grid(zip_cuminf_017, zip_cuminf_1849,zip_cuminf_5064, zip_cuminf_65, nrow = 2, align = 'hv') -> age_zip_cuminf_map
save_plot('figs/s6_age-zip-cuminf-map.png', age_zip_cuminf_map, base_height = 5, base_asp = 1.7,bg = "white")

## Get the age ZIP reporting rate maps

age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  left_join(aph_zip_age_tests, by = c('age_group_paper', 'zip')) %>% 
  mutate(est_infections = ifelse(est_infections < cases, cases, est_infections)) %>% 
  mutate(reporting_rate = samp_conditional_reporting_rate(cases, est_infections, n())) %>% 
  group_by(population, zip, age_group_paper) %>% 
  summarize(rr_mean = mean(reporting_rate),
            rr_lower = quantile(reporting_rate, probs = 0.025),
            rr_upper = quantile(reporting_rate, probs = 0.975)) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) -> age_zip_rr_df

age_zip_rr_df %>% 
  filter(age_group_paper == '0-17') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = rr_mean, 
                plot_title = '0-17 Reporting rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_rr_017
age_zip_rr_df %>% 
  filter(age_group_paper == '18-49') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = rr_mean, 
                plot_title = '18-49 Reporting rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_rr_1849
age_zip_rr_df %>% 
  filter(age_group_paper == '50-64') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = rr_mean, 
                plot_title = '50-64 Reporting rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_rr_5064
age_zip_rr_df %>% 
  filter(age_group_paper == '65+') %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = rr_mean, 
                plot_title = '65+ Reporting rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_rr_65


plot_grid(zip_rr_017, zip_rr_1849,zip_rr_5064, zip_rr_65, nrow = 2, align = 'hv') -> age_zip_rr_map
save_plot('figs/s6_age-zip-rr-map.png', age_zip_rr_map, base_height = 5, base_asp = 1.7,bg = "white")


# Look at relationship between infection and reporting rates --------------
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(zip, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  left_join(zip_cov_df) %>% 
  mutate(cases = ifelse(inf < cases, inf, cases)) %>% 
  mutate(rr = samp_conditional_reporting_rate(detected_cases = cases, 
                                              estimated_infections = inf, 
                                              num_samps = n())) -> zip_inf_rr_df

get_nls_fit <- function(data){
  # browser()
  data %>% 
    nls(data = ., formula = rr ~ a*inf^b, start = c(a=10, b=-0.5)) %>% 
    broom::tidy()  
}
## Infections calculated per 100,000 so it can be compared with this paper:
## https://www.nature.com/articles/s41586-020-03095-6
zip_inf_rr_df %>% 
  mutate(inf = inf/population*100000) %>%
  ungroup() %>% 
  select(samp, inf, rr) %>% 
  nest(data = -samp) %>%
  mutate(fitnls = map(data, get_nls_fit)) %>% 
  unnest(fitnls) %>% 
  group_by(term) %>% 
  summarize(estimate_mean = mean(estimate),
            estimate_lo = quantile(estimate, probs=0.025),
            estimate_hi = quantile(estimate, probs = 0.975)) -> nls_rr_inf_est
nls_rr_inf_est

nls_mean_curve <- tibble(inf = seq(0,100000, by = 10),
                         rr = nls_rr_inf_est$estimate_mean[1]*inf^nls_rr_inf_est$estimate_mean[2])  

zip_inf_rr_df %>% 
  mutate(inf = inf/population*100000) %>%
  group_by(zip) %>% 
  summarize(inf_mean = mean(inf),
            inf_lo = quantile(inf, probs = 0.025),
            inf_hi = quantile(inf, probs = 0.975),
            rr_mean = mean(rr),
            rr_lo = quantile(rr, probs = 0.025),
            rr_hi = quantile(rr, probs = 0.975)) %>% 
  ggplot(aes(inf_mean, rr_mean)) + 
    geom_point() +
    geom_errorbar(aes(ymin = rr_lo, ymax = rr_hi), alpha = .5) +
    geom_errorbarh(aes(xmin = inf_lo, xmax = inf_hi), alpha = .5) +
    coord_cartesian(ylim = c(0,1), xlim = c(0, 75000)) +
    geom_line(data = nls_mean_curve, aes(x = inf, y = rr), lty = 2, inherit.aes = FALSE) +
    labs(x = 'Infections per 100,000', y = 'Reporting Rate (%)') +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::comma) +
    background_grid(major = 'xy') -> sfig_rr_inf_relationship
save_plot('figs/sfig_rr_inf_relationship.png', sfig_rr_inf_relationship, base_height = 5, base_asp = 1.3, bg = "white")

## Compare SVI relationship with reported case and hospitalization numbers
rm(list=ls())
library(lme4)
library(broom)
library(broom.mixed)
load('processed-data/fig4-image.rda')

glmer(cases ~ svi + (1|zip), 
      data = zip_cov_df %>% 
        mutate(case_rate = cases/population*100000),
      offset = log(population),
      family = 'poisson') -> case_svi_glmer_mod
confint(case_svi_glmer_mod, parm = 'svi') -> case_svi_glmer_ci

## Relative risk estimates for highest SVI to lowest reporting rate
exp(c(coef(case_svi_glmer_mod)$zip[1,'svi'], case_svi_glmer_ci)*lmer_max_diff)

case_svi_boot <- bootMer(case_svi_glmer_mod, FUN = boot_pred, nsim = 100)
new_dat_glmer %>% 
  bind_cols(case_svi_boot$t %>% t() %>% as_tibble()) %>% 
  gather(samp, pred, -svi, -zip) %>% 
  group_by(svi) %>% 
  summarize(medpred = median(pred),
            lower = quantile(pred, probs = .025),
            upper = quantile (pred, probs = .975)) -> case_svi_pred_summary

zip_cov_df %>% 
  mutate(case_rate = cases/population*100000) %>% 
  ggplot(aes(svi, case_rate)) + 
  geom_ribbon(data = case_svi_pred_summary, aes(x = svi, ymin = lower*100000, ymax = upper*100000), 
              inherit.aes = F, color = NA, fill = 'blue', alpha = .4) +
  geom_line(data = case_svi_pred_summary, aes(x = svi, y = medpred*100000), 
            inherit.aes = F, color = 'blue', size =1) +
  theme(strip.background= element_rect(fill=NA), 
        strip.text = element_text(size = 16)) +
  geom_point() +
  labs(x = 'SVI', y = 'Case rate (per 100k)') +
  theme(legend.position = 'none') +
  background_grid(major = ('xy')) -> case_rate_svi_plot
case_rate_svi_plot

age_zip_inf_df %>% 
  group_by(zip) %>% 
  summarize(admits = sum(n_admits)) %>% 
  left_join(zip_cov_df, by = 'zip') %>% 
  mutate(admit_rate = admits/population*100000) -> hosp_svi_df

glmer(admits ~ svi + (1|zip), 
      data = hosp_svi_df %>% filter(zip != '78701'),
      offset = log(population),
      family = 'poisson') -> hosp_svi_glmer_mod
confint(hosp_svi_glmer_mod, parm = 'svi') -> hosp_svi_glmer_ci

## Relative risk estimates for highest SVI to lowest reporting rate
exp(c(coef(hosp_svi_glmer_mod)$zip[1,'svi'], hosp_svi_glmer_ci)*lmer_max_diff)

hosp_svi_boot <- bootMer(hosp_svi_glmer_mod, FUN = boot_pred, nsim = 100)
new_dat_glmer %>% 
  bind_cols(hosp_svi_boot$t %>% t() %>% as_tibble()) %>% 
  gather(samp, pred, -svi, -zip) %>% 
  group_by(svi) %>% 
  summarize(medpred = median(pred),
            lower = quantile(pred, probs = .025),
            upper = quantile (pred, probs = .975)) -> hosp_svi_pred_summary

hosp_svi_df %>% 
  ggplot(aes(svi, admit_rate)) + 
  geom_ribbon(data = hosp_svi_pred_summary, aes(x = svi, ymin = lower*100000, ymax = upper*100000), 
              inherit.aes = F, color = NA, fill = 'blue', alpha = .4) +
  geom_line(data = hosp_svi_pred_summary, aes(x = svi, y = medpred*100000), 
            inherit.aes = F, color = 'blue', size =1) +
  geom_point() +
  theme(strip.background= element_rect(fill=NA), 
        strip.text = element_text(size = 16)) +
  labs(x = 'SVI', y = 'Hospitalization rate (per 100k)') +
  theme(legend.position = 'none') +
  background_grid(major = ('xy')) -> hosp_rate_svi_plot
hosp_rate_svi_plot


plot_grid(case_rate_svi_plot, hosp_rate_svi_plot, align = 'hv', nrow = 1) -> reported_metric_svi_zip_plot
save_plot('figs/covidburden_svi_zip.png', reported_metric_svi_zip_plot, base_height = 4, base_asp = 2, bg = 'white')


## Comparing age group zip SVI infection relationships
age_zip_inf_df %>% 
  unnest(est_infections) %>%
  left_join(aph_zip_age_tests, by = c('age_group_paper', 'zip')) %>% 
  mutate(est_infections = ifelse(est_infections < cases, cases, est_infections)) %>% 
  left_join(zip_cov_df %>% select(zip, svi), by = 'zip') %>% 
  mutate(cases = ifelse(cases>est_infections,est_infections, cases)) %>% 
  select(zip, age_group_paper, population, est_infections, svi, cases) -> zip_age_svi_glmer_df

glmer(est_infections ~ svi*age_group_paper + (1|zip), 
      data = zip_age_svi_glmer_df,
      offset = log(population),
      family = 'poisson') -> inf_age_svi_glmer_mod
summary(inf_age_svi_glmer_mod)

new_dat_glmer %>% 
  expand_grid(age_group_paper = unique(zip_age_svi_glmer_df$age_group_paper)) -> new_dat_age_glmer

broom::tidy(inf_age_svi_glmer_mod) %>% 
  filter(effect == 'fixed') %>% 
  filter(grepl('svi', term)) %>% 
  select(term, estimate, sd = std.error) %>% 
  mutate(samps = map2(estimate, sd, ~rnorm(1000, mean = .x, sd = .y))) %>% 
  unnest(samps) %>% 
  group_by(term) %>% 
  mutate(samp = seq_along(term)) %>% 
  ungroup() %>% 
  select(term, samp, samps) %>% 
  spread(term, samps) %>% 
  mutate_at(vars(-samp,-svi), .funs = ~.+svi) %>% 
  gather(term, value, -samp) %>% 
  mutate(value = exp(value*lmer_max_diff)) %>% 
  group_by(term) %>% 
  summarize(avg_coef = mean(value),
            lo_coef = quantile(value, probs = .025),
            hi_coef = quantile(value, probs = .975))



new_dat_age_glmer %>% 
  bind_cols(merTools::predictInterval(inf_age_svi_glmer_mod, 
                                      newdata = new_dat_age_glmer, 
                                      which = 'fixed', 
                                      level = 0.95) %>% exp()) %>% 
  mutate_at(vars(fit, upr, lwr), ~ifelse(.>1, 1, .)) -> age_zip_inf_pred_conf

## Reporting rate estimation
glmer(cases ~ svi*age_group_paper + (1|zip), 
      data = zip_age_svi_glmer_df,
      offset = log(est_infections+.5),
      family = 'poisson') -> rr_age_svi_glmer_mod
summary(rr_age_svi_glmer_mod)
broom::tidy(rr_age_svi_glmer_mod) %>% 
  filter(effect == 'fixed') %>% 
  filter(grepl('svi', term)) %>% 
  select(term, estimate, sd = std.error) %>% 
  mutate(samps = map2(estimate, sd, ~rnorm(1000, mean = .x, sd = .y))) %>% 
  unnest(samps) %>% 
  group_by(term) %>% 
  mutate(samp = seq_along(term)) %>% 
  ungroup() %>% 
  select(term, samp, samps) %>% 
  spread(term, samps) %>% 
  mutate_at(vars(-samp,-svi), .funs = ~.+svi) %>% 
  gather(term, value, -samp) %>% 
  mutate(value = exp(value*lmer_max_diff)) %>% 
  group_by(term) %>% 
  summarize(avg_coef = mean(value),
            lo_coef = quantile(value, probs = .025),
            hi_coef = quantile(value, probs = .975))

new_dat_age_glmer %>% 
  bind_cols(merTools::predictInterval(rr_age_svi_glmer_mod, 
                                      newdata = new_dat_age_glmer, 
                                      which = 'fixed', 
                                      level = 0.95) %>% exp()) %>% 
  mutate_at(vars(fit, upr, lwr), ~ifelse(.>1, 1, .)) -> age_zip_rr_pred_conf




zip_age_svi_glmer_df %>% 
  group_by(zip, age_group_paper,svi) %>% 
  summarize(infrate_mean = mean(est_infections/population),
            infrate_lo = quantile(est_infections/population, probs = 0.025),
            infrate_hi = quantile(est_infections/population, probs=0.975),
            rr_mean = mean(cases/est_infections),
            rr_lo = quantile(cases/est_infections, probs = 0.025),
            rr_hi = quantile(cases/est_infections, probs = 0.975)) %>% 
  ungroup() %>% 
  ggplot(aes(svi, infrate_mean)) +
  facet_wrap(~age_group_paper) +
  geom_line(data = age_zip_inf_pred_conf, aes(x = svi, y = fit), color = 'blue') +
  geom_ribbon(data = age_zip_inf_pred_conf, aes(x = svi, ymin = lwr, ymax = upr),
              fill = 'blue', color = NA, alpha = .3, inherit.aes=F) +
  geom_point() +
  geom_errorbar(aes(ymin = infrate_lo, ymax=infrate_hi)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1) +
  scale_y_continuous(labels = scales::percent) +
  theme(strip.background= element_rect(fill=NA), 
        axis.line = element_blank(), 
        strip.text = element_text(size = 16)) +
  labs(x = 'SVI', y = 'Cumulative infections (%)') +
  background_grid(major = ('xy')) ->  zip_age_infrate_svi_plot
zip_age_infrate_svi_plot
save_plot('figs/zip-age-infrate-svi.png', zip_age_infrate_svi_plot, base_height = 6, base_asp = 1.2, bg = 'white')


zip_age_svi_glmer_df %>% 
  group_by(zip, age_group_paper,svi) %>% 
  summarize(infrate_mean = mean(est_infections/population),
            infrate_lo = quantile(est_infections/population, probs = 0.025),
            infrate_hi = quantile(est_infections/population, probs=0.975),
            rr_mean = mean(cases/est_infections),
            rr_lo = quantile(cases/est_infections, probs = 0.025),
            rr_hi = quantile(cases/est_infections, probs = 0.975)) %>% 
  ggplot(aes(svi, rr_mean)) +
  facet_wrap(~age_group_paper) +
  geom_line(data = age_zip_rr_pred_conf, aes(x = svi, y = fit), color = 'blue') +
  geom_ribbon(data = age_zip_rr_pred_conf, aes(x = svi, ymin = lwr, ymax = upr), 
              fill = 'blue', color = NA, alpha = .3, inherit.aes=F) +
  geom_point() +
  geom_errorbar(aes(ymin = rr_lo, ymax=rr_hi)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1) +
  scale_y_continuous(labels = scales::percent) +
  theme(strip.background= element_rect(fill=NA), 
        axis.line = element_blank(), 
        strip.text = element_text(size = 16)) +
  labs(x = 'SVI', y = 'Reporting Rate (%)') +
  background_grid(major = ('xy')) -> zip_age_rr_svi_plot
zip_age_rr_svi_plot
save_plot('figs/zip-age-rr-svi.png', zip_age_rr_svi_plot, base_height = 6, base_asp = 1.2, bg = 'white')


# Comparison between case,hosp,inf models ---------------------------------


exp(c(coef(case_svi_glmer_mod)$zip[1,'svi'], case_svi_glmer_ci)*lmer_max_diff)
exp(c(coef(inf_glmer_mod)$zip[1,'svi'], inf_svi_glmer_ci)*lmer_max_diff)
exp(c(coef(hosp_svi_glmer_mod)$zip[1,'svi'], hosp_svi_glmer_ci)*lmer_max_diff)

age_zip_inf_df %>% 
  group_by(age_group_paper, zip,population) %>% 
  summarize(admits = sum(n_admits)) %>% 
  left_join(zip_cov_df %>% select(zip,svi), by = 'zip') %>% 
  mutate(admit_rate = admits/population*100000) -> hosp_age_svi_df
glmer(admits ~ svi*age_group_paper + (1|zip), 
      data = hosp_age_svi_df,
      offset = log(population),
      family = 'poisson') -> hosp_age_svi_glmer_mod
summary(hosp_age_svi_glmer_mod)

age_zip_inf_df %>%
  select(zip,age_group_paper,population) %>%
  left_join(aph_zip_age_tests, by = c('zip', 'age_group_paper')) %>%
  left_join(zip_cov_df %>% select(zip, svi), by = 'zip') -> case_age_zvi_df
glmer(cases ~ svi*age_group_paper + (1|zip),
      data = case_age_zvi_df,
      offset = log(population),
      family = 'poisson') -> case_age_svi_glmer_mod
summary(case_age_svi_glmer_mod)


bind_rows(broom::tidy(rr_age_svi_glmer_mod) %>% 
              mutate(model = 'Reporting rate'),
            broom::tidy(case_age_svi_glmer_mod) %>%
              mutate(model = 'Cases'),
            broom::tidy(hosp_age_svi_glmer_mod) %>% 
              mutate(model = 'Hospital admissions'),
            broom::tidy(inf_age_svi_glmer_mod) %>% 
              mutate(model = 'Infections')) %>% 
  filter(effect == 'fixed') %>% 
  filter(grepl('svi', term)) %>% 
  select(model, term, estimate, sd = std.error) %>% 
  mutate(samps = map2(estimate, sd, ~rnorm(1000, mean = .x, sd = .y))) %>% 
  unnest(samps) %>% 
  group_by(model, term) %>% 
  mutate(samp = seq_along(term)) %>% 
  ungroup() %>% 
  select(model, term, samp, samps) %>% 
  spread(term, samps) %>% 
  mutate_at(vars(-model,-samp,-svi), .funs = ~.+svi) %>% 
  gather(term, value, -samp, -model) %>% 
  mutate(value = exp(value*lmer_max_diff)) %>% 
  group_by(model, term) %>% 
  summarize(avg_coef = mean(value),
            lo_coef = quantile(value, probs = .025),
            hi_coef = quantile(value, probs = .975)) %>% 
  mutate(cleaned = paste0(round(avg_coef,  digits = 2), ' (',
                          round(lo_coef,  digits = 2), '-',
                          round(hi_coef,  digits = 2), ')'
                          )) %>% 
  select(model, term, cleaned) %>% 
  spread(model, cleaned) %>% 
  select(term, `Reporting rate`,Cases, Infections, `Hospital admissions`) ->age_zip_disparity_results
age_zip_disparity_results

write_csv(age_zip_disparity_results, 'processed-data/age_zip_disparity_results.csv')


# Make case rate map between APH and overall ------------------------------
load('raw-data/spatial-files/zip_geometries.rda')
 # needed for mz_ functions to get roads
mz_set_tile_host_nextzen(key="hxNDKuWbRgetjkLAf_7MUQ")
## Get roads for maps
get_vector_tiles <- function(bbox){
  mz_box=mz_rect(bbox$xmin,bbox$ymin,bbox$xmax,bbox$ymax)
  mz_vector_tiles(mz_box)
}
zcta_geom <- st_union(zip_cov_df %>% 
                        left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
                        pull(geometry))
zcta_bbox <- st_bbox(zcta_geom)
zcta_vector_tiles <- get_vector_tiles(zcta_bbox)
zcta_roads <- as_sf(zcta_vector_tiles$roads) 

I35<- zcta_roads %>%
  mutate(st_transform(geometry, st_crs(zip_geom))) %>% 
  filter(ref== "I 35" | ref=="I 35;US 190" | ref=="I 35;US 290" | ref=="I 35;US 77") %>% # | ref=="I 35-H Business") %>%
  pull(geometry) %>%
  st_union() %>%
  st_transform(st_crs(zcta_geom)) %>%
  st_intersection(zcta_geom)

US183 = zcta_roads %>%
  mutate(st_transform(geometry, st_crs(zip_geom))) %>% 
  # filter(grepl('183', ref)) %>% 
  filter(ref=="US 183" | ref=="US 183;FM 20" | ref=="US 183;TX 80" | ref=="US 90;US 183" | ref== "US 183;US 190;US 281"
         | ref=="US 183;US 190;US 281;FM 580" | ref=="US 183;US 190;US 281;Truck" | ref=="US 183;US 190" | id=="2c2066c3b6fafe2574463154f5cc877d" ) %>%
  pull(geometry) %>%
  st_union() %>%
  st_transform(st_crs(zcta_geom)) %>%
  st_intersection(zcta_geom)

plot_zip_maps <- function(sf_obj, fill_col, plot_title, scales_func){
  fill_col <- enquo(fill_col)
  # browser()
  sf_obj %>% 
    ggplot() +
    geom_sf(aes(geometry=geometry, fill=!!fill_col), 
            size = 0.05, color="black")+ # , show.legend = FALSE
    scale_fill_gradient(low="lightyellow", high="slateblue4", labels = scales_func)+
    guides(fill = guide_colourbar(frame.colour = "black", ticks.colour="black"))+
    geom_sf(data = I35, col = "black", size=0.3)+
    geom_sf(data = US183, col = "black", size=0.3)+
    annotate(geom="text", x=-97.68, y=30.50, label="I-35", size=3)+
    annotate(geom="text", x=-97.80, y=30.50, label="US 183", size=3)+
    labs(title = plot_title, fill = NULL)+
    theme_map() +
    theme(plot.title = element_text(hjust = 0.5,size = 14), plot.title.position = 'plot') +
    NULL
}
  

## Case rate
age_zip_inf_df %>%
  select(zip,age_group_paper,population) %>%
  left_join(aph_zip_age_tests, by = c('zip', 'age_group_paper')) %>%
  left_join(zip_cov_df %>% select(zip, svi), by = 'zip') -> case_age_zvi_df
case_age_zvi_df %>% 
  group_by(zip) %>% 
  summarize(cases = sum(cases)) %>% 
  left_join(zip_cov_df %>% 
              rename(all_cases = cases), by = 'zip') %>% 
  mutate(case_ratio = cases/all_cases) %>% 
  select(zip, cases, all_cases, case_ratio) %>% 
  # gather(key, value, cases:all_cases) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  # filter(key == 'all_cases') %>% 
  ggplot() +
    geom_sf(aes(geometry=geometry, fill=case_ratio), 
            size = 0.05, color="black")+ # , show.legend = FALSE
    scale_fill_gradient(low="lightyellow", high="slateblue4", labels = scales::label_comma())+
    guides(fill = guide_colourbar(frame.colour = "black", ticks.colour="black"))+
    geom_sf(data = I35, col = "black", size=0.3)+
    geom_sf(data = US183, col = "black", size=0.3)+
    annotate(geom="text", x=-97.68, y=30.50, label="I-35", size=3)+
    annotate(geom="text", x=-97.80, y=30.50, label="US 183", size=3)+
    labs(title = NULL, fill = 'Proportion \nof cases')+
    theme_map() +
    # facet_wrap(~key) +
    # theme(plot.title = element_text(hjust = 0.5,size = 14), plot.title.position = 'plot') +
    NULL -> zip_case_ratio_map
zip_case_ratio_map  


case_age_zvi_df %>%
  group_by(zip) %>%
  summarize(cases = sum(cases)) %>%
  left_join(zip_cov_df %>%
              rename(all_cases = cases), by = 'zip') %>%
  mutate(case_ratio = cases/all_cases) %>%
  select(zip, aph_cases = cases, all_cases, case_ratio) %>%
  left_join(zip_cov_df) %>%
  ggplot(aes(svi, case_ratio)) +
  geom_point() +
  stat_smooth(method = 'lm', se = F) +
  labs(x = 'SVI', y = 'Fraction of cases') +
  background_grid(major = 'xy') -> zip_case_ratio_scatter
zip_case_ratio_scatter
plot_grid(zip_case_ratio_map, zip_case_ratio_scatter, labels = 'AUTO', rel_widths = c(1.8,1)) -> zip_case_ratio
save_plot('figs/zip_case_ratio.png', zip_case_ratio, base_height = 4, base_asp = 3, bg = 'white')
# Make the hospitalization rate time-series  plot ---------------------------
# age_zip_inf_df %>% 
#   filter(n_admits >0) %>% 
#   select(zip, admission_dates) %>% 
#   unnest(admission_dates) %>% 
#   mutate(date_group = cut.Date(admission_dates, breaks = "4 week", labels = FALSE)) %>% 
#   group_by(date_group) %>% 
#   mutate(date = min(admission_dates)) %>% 
#   ungroup() %>% 
#   left_join(svi_zip_rank_df, by = 'zip') %>% 
#   group_by(date, date_group, svi_label) %>% 
#   summarize(admits = n()) %>% 
#   left_join(svi_zip_label_pop_df, by = 'svi_label') %>% 
#   mutate(admit_per_cap = admits/pop) %>% 
#   group_by(date) %>% 
#   summarize(hrelrisk = admit_per_cap[svi_label == 'mostvulnerable']/admit_per_cap[svi_label=='leastvulnerable']) %>% 
#   filter(date != '2021-05-31') %>% 
#   ggplot(aes(date, hrelrisk)) + geom_point()
# 
# age_zip_inf_df %>% 
#   filter(n_admits >0) %>% 
#   select(zip, admission_dates) %>% 
#   unnest(admission_dates) %>% 
#   mutate(date_group = cut.Date(admission_dates, breaks = "4 week", labels = FALSE)) %>% 
#   group_by(date_group) %>% 
#   mutate(date = min(admission_dates)) %>% 
#   ungroup() %>% 
#   left_join(svi_zip_rank_df, by = 'zip') %>% 
#   group_by(svi_label) %>% 
#   summarize(admits = n()) %>% 
#   left_join(svi_zip_label_pop_df, by = 'svi_label') %>% 
#   mutate(admit_per_cap = admits/pop) %>%
#   summarize(hrelrisk = admit_per_cap[svi_label == 'mostvulnerable']/admit_per_cap[svi_label=='leastvulnerable'])


read_csv('raw-data/ihr-comparison.csv') %>% 
  gather(key, age, age_lo, age_hi) %>%
  mutate(age = as.integer(age)) %>% 
  group_by(region) %>% 
  arrange(region, age) %>% 
  padr::pad_int(by = 'age', group = 'region') %>% 
  select(-key) %>% 
  tidyr::fill(mean, lo, hi, study) %>% 
  mutate_at(vars(mean,lo,hi), .funs = ~./100) %>% 
  ggplot(aes(age, mean, color = region)) + 
  geom_line() + 
  scale_y_log10(breaks = c(.001, 0.0035, .01, .035, .1, .35),
                labels = scales::label_percent(accuracy = .01)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = region), alpha = .2, color = NA) +
  scale_color_brewer(type = 'qual', palette = 2) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(x = 'Age', y = 'Infection hospitalization rate (%)', fill = '', color = '') +
  background_grid(major = 'xy', minor = 'x') -> ihr_comparison_plot

save_plot('figs/ihr_comparison_plot.png', ihr_comparison_plot, 
          base_height = 4.5, base_asp = 1.6, bg = 'white')



save.image('processed-data/si-image.rda')
