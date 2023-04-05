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
# ## Linear regression for infections and svi
# get_linear_regression <- function(df, response_col){
#   lm(get(response_col)~svi, data = df) -> mod
#   tibble(svi_impact = mod$coefficients['svi'],
#          svi_intercept = mod$coefficients['(Intercept)'],
#          rsquared = summary(mod)$r.squared,
#          pval = summary(mod)$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}) -> results
#   colnames(results) <- paste0(response_col,'_', colnames(results))
#   results
# }
# 
# get_mean_ci <- function(vec,narm_arg=T){
#   return(c(mean = mean(vec, na.rm = narm_arg), 
#            lo = quantile(vec, probs = 0.025, na.rm = narm_arg) %>% unname(),
#            hi = quantile(vec, probs = 0.975, na.rm = narm_arg) %>% unname()))
# }

age_zip_inf_df %>% 
  select(zip,age_group_paper,est_infections) %>% 
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
  select(id, zip, inf_rate, cases, rr, svi, infections, population) -> lmer_df

## Compare model estimates and plots for cumulative infections
glmer(infections ~ svi + (1|zip), 
      data = lmer_df, 
      offset = log(population), 
      family = 'poisson')  -> inf_glmer_mod
confint(inf_glmer_mod, parm = 'svi') -> inf_svi_glmer_ci
effects::effect(term= "svi", mod= inf_glmer_mod) %>% 
  as.data.frame %>% 
  mutate_at(vars(fit:upper), ~./mean(zip_cov_df$population)) -> inf_svi_lmer_df
# summary(inf_glmer_mod)
# sjPlot::plot_model(inf_glmer_mod)
# sjPlot:: tab_model(inf_glmer_mod)
svi_vals_for_lmer_plot <- lmer_df %>% pull(svi) %>% 
  unique() %>% 
  quantile(probs = seq(0,1, by = .1), type = 3) %>% 
  unname()

new_dat_glmer <- lmer_df %>% 
  filter(svi %in% svi_vals_for_lmer_plot) %>% 
  distinct(zip,svi) %>% 
  arrange(svi)
lmer_max_diff <- diff(quantile(zip_cov_df$svi, probs = c(0.25, 0.75))) %>% unname()

## Relative risk estimates for highest SVI to lowest infections
exp(c(coef(inf_glmer_mod)$zip[1,'svi'], inf_svi_glmer_ci)*lmer_max_diff)

boot_pred <- function(mod, 
                      ndata=new_dat_glmer){
  exp(predict(mod, newdata = ndata, re.form = NA))
}

inf_boot <- bootMer(inf_glmer_mod, FUN = boot_pred, nsim = 100)
new_dat_glmer %>% 
  bind_cols(inf_boot$t %>% t() %>% as_tibble()) %>% 
  gather(samp, pred, -svi, -zip) %>% 
  group_by(svi) %>% 
  summarize(medpred = median(pred),
            lower = quantile(pred, probs = .025),
            upper = quantile (pred, probs = .975)) ->inf_pred_summary

## Compare model estimates and plots for the reporting rate
glmer(cases ~ svi + (1|zip), 
      data = lmer_df , 
      offset = log(infections), 
      family = 'poisson')  -> rr_glmer_mod
summary(rr_glmer_mod)
confint(rr_glmer_mod, parm = 'svi') -> rr_svi_glmer_ci
# summary(rr_glmer_mod)
# sjPlot::plot_model(rr_glmer_mod)
# sjPlot:: tab_model(rr_glmer_mod)

## Relative risk estimates for highest SVI to lowest reporting rate
exp(c(coef(rr_glmer_mod)$zip[1,'svi'], rr_svi_glmer_ci)*lmer_max_diff)

rr_boot <- bootMer(rr_glmer_mod, FUN = boot_pred, nsim = 100)
new_dat_glmer %>% 
  bind_cols(rr_boot$t %>% t() %>% as_tibble()) %>% 
  gather(samp, pred, -svi, -zip) %>% 
  group_by(svi) %>% 
  summarize(medpred = median(pred),
            lower = quantile(pred, probs = .025),
            upper = quantile (pred, probs = .975)) -> rr_pred_summary


## Interpretation of poisson regression ideas
# infections ~ P(lambda = exp(b1*svi + log(population) + b0))
# exp(b1*svi)exp(log(population))exp(b0)
# infections/population = ~

## Get summary data for making correlation plots
lmer_df %>% 
  group_by(zip) %>% 
  summarize(svi = unique(svi),
            inf_mean = mean(inf_rate),
            inf_lo = quantile(inf_rate, probs = 0.025),
            inf_hi = quantile(inf_rate, probs = 0.975),
            rr_mean = mean(rr), 
            rr_lo = quantile(rr, probs = 0.025),
            rr_hi = quantile(rr, probs = 0.975)) -> zip_corr_df


zip_corr_df %>% 
  ggplot(aes(svi, inf_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = inf_lo, ymax = inf_hi)) +
  geom_ribbon(data = inf_pred_summary, aes(x = svi, ymin = lower, ymax = upper), 
              inherit.aes = F, color = NA, fill = 'blue', alpha = .4) +
  geom_line(data = inf_pred_summary, aes(x = svi, y = medpred), 
            inherit.aes = F, color = 'blue', size =1) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "SVI", y = 'Cumulative infections (%)') +
  theme(legend.position = 'none') +
  background_grid(major = 'xy') -> cum_inf_corr_plot
cum_inf_corr_plot  

# rr_svi_mod_summary <- rr_svi_lmer_df %>% 
#   group_by(svi) %>% 
#   summarize(medpred = median(pred),
#             lower = quantile(pred, probs = .025),
#             upper = quantile (pred, probs = .975))
zip_corr_df %>% 
  ggplot(aes(svi, rr_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = rr_lo, ymax = rr_hi)) +
  geom_ribbon(data = rr_pred_summary, aes(x = svi, ymin = lower, ymax = upper), 
              inherit.aes = F, color = NA, fill = 'blue', alpha = .4) +
  geom_line(data = rr_pred_summary, aes(x = svi, y = medpred), 
            inherit.aes = F, color = 'blue', size =1) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "SVI", y = 'Reporting rate (%)') +
  theme(legend.position = 'none') +
  background_grid(major = 'xy') -> rr_corr_plot
rr_corr_plot  

plot_grid(cum_inf_corr_plot, rr_corr_plot, nrow = 1, labels = 'AUTO') -> svi_corr_plots
save_plot('figs/svi_corr_plots.png', svi_corr_plots, base_height = 5, base_asp = 2,bg = "white")


## Calculates infection rate and reporting rate for 4 week time periods
## Zip daily reported cases
zip_cov_df %>% 
  unnest(ts_cases) %>% 
  group_by(zip) %>% 
  mutate(new_cases = diff(c(0,cum_cases))) %>% 
  mutate(new_cases = ifelse(new_cases<0,0,new_cases)) -> zip_daily_cases_cov_df


# new_inf_all_zip_days <- age_zip_inf_df %>% 
#   unnest(inf_timing) %>% 
#   select(zip, id, date, new_infections) %>% 
#   filter(date >= first_infection_date, date < final_infection_date) %>% 
#   group_by(id, date, zip) %>% 
#   summarize(new_inf = sum(new_infections)) %>% 
#   group_by(id, zip) %>% 
#   padr::pad(start_val = ymd(first_infection_date), end_val = ymd(final_infection_date), break_above = 22) %>% 
#   ungroup() %>% 
#   mutate(new_inf = ifelse(is.na(new_inf), 0, new_inf))
# new_inf_all_zip_days %>% 
#   left_join(zip_cov_df %>% select(zip, svi,population), by = 'zip') %>% 
#   mutate(time = difftime(date, min(date), units = 'days') %>% as.numeric) -> new_inf_all_zip_days
# 
# 
# glmer(new_inf ~ svi + svi:time + (1|zip), 
#       data = new_inf_all_zip_days %>% filter(zip != '78701'), ## Removes downtown because of outlier 
#       offset = log(population), 
#       family = 'poisson')  -> full_mod
# summary(full_mod)

## Calculate the 4 week infections and new cases
interval_breaks <- seq.Date(from = ymd(first_infection_date), to = ymd(final_infection_date), by = '4 week')
age_zip_inf_df %>% 
  unnest(inf_timing) %>% 
  select(zip, id, date, new_infections) %>% 
  filter(date >= first_infection_date, date < final_infection_date) %>% 
  group_by(id, date, zip) %>% 
  summarize(new_inf = sum(new_infections)) %>% 
  group_by(id, zip) %>% 
  padr::pad(start_val = ymd(first_infection_date), end_val = ymd(final_infection_date), break_above = 22) %>% 
  ungroup() %>% 
  mutate(new_inf = ifelse(is.na(new_inf), 0, new_inf)) %>% 
  mutate(date_group = cut.Date(date, breaks =interval_breaks)) %>% 
  left_join(zip_daily_cases_cov_df %>% select(zip, date, new_cases), by = c('zip', 'date')) %>% 
  group_by(id, date_group, zip) %>% 
  summarize(new_inf = sum(new_inf),
            new_cases = sum(new_cases, na.rm=T)) %>% 
  mutate(new_cases = ifelse(is.na(new_cases), 0, new_cases)) -> zip_4week_infections_cases_df


get_zip_lmer  <- function(df, date_group, response_col, offset_col){
  # browser()
  
  ## Removes zeroes from new infections, because add no information to
  ## model, and they mess up computations (causes errors sometimes...)
  ## See: https://stackoverflow.com/questions/72556746/offsetissue-error-message-occurred-when-fitting-a-poisson-mixed-effect-mode
  ## Also can be thought of as ZIP codes not participating, 
  ## see: https://blog.stata.com/2011/08/22/use-poisson-rather-than-regress-tell-a-friend/
  ## Zeroes only occur in smaller ZIP codes, so doesn't dramatically impact results, just stabilizes estimates
  # browser()
  df <- df %>%
    filter(inf_rate!=0)
  
  glmer(get(response_col) ~ svi + (1|zip), 
        data = df,
        offset = log(get(offset_col)), 
        family = 'poisson')  -> mod
  
  confint(mod, parm = 'svi', method= 'Wald') -> ci
  
  ## Relative risk estimates for highest SVI to lowest reporting rate
  svi_impact_vec <- exp(c(coef(mod)$zip[1,'svi'], ci)*lmer_max_diff)
  tibble(svi_slope_mean = svi_impact_vec[1],
         svi_slope_lo = svi_impact_vec[2],
         svi_slope_hi = svi_impact_vec[3]) -> results
  colnames(results) <- paste0(response_col,'_', colnames(results))
  print(results)
  results
}
zip_4week_infections_cases_df %>%
  filter(!is.na(date_group)) %>% 
  left_join(zip_cov_df %>% select(zip,population,svi),by = 'zip') %>% 
  ungroup() %>% 
  mutate(new_cases = ifelse(new_inf < new_cases, new_inf, new_cases)) %>% 
  mutate(inf_rate = new_inf/population,
         reporting_rate = samp_conditional_reporting_rate(new_cases, new_inf, n())) %>% 
  select(date_group, zip, new_inf, new_cases, inf_rate, reporting_rate, svi, population) %>% 
  nest(data = c(zip, new_inf, new_cases, inf_rate, reporting_rate, svi, population)) %>% 
  mutate(date_group = ymd(as.character(date_group))) %>% 
  mutate(infrate_lms = map2(data, date_group, get_zip_lmer, response_col = 'new_inf', offset_col = 'population'),
         rr_lms = map2(data, date_group, get_zip_lmer, response_col = 'new_cases', offset_col = 'new_inf')) %>% 
  select(-data) %>% 
  unnest(cols = c('rr_lms', 'infrate_lms')) %>% 
  mutate(date_group = ymd(as.character(date_group))) -> lmer_4week_results

## Track the number of zip codes with 0 infection estimates which cause an error in estimation
## Above 10% of ZIP codes gives unreliable estimates.
# zip_4week_infections_cases_df %>%
#   left_join(zip_cov_df %>% select(zip,population,svi),by = 'zip') %>% ungroup() %>%
#   group_by(date_group, zip, svi) %>%
#   summarize(new_inf_rate = mean(new_inf/population),
#             new_inf_rate_lo = quantile(new_inf/population, probs = 0.025),
#             new_inf_rate_hi = quantile(new_inf/population, probs = 0.975),
#             new_inf = mean(new_inf),
#             population = unique(population)) %>% 
#   group_by(date_group) %>% 
#   summarize(sum(new_inf_rate==0)/n())

make_4week_plot <- function(lmer_results, plott = 'inf'){
  # browser()
  if(plott == 'inf'){
    lmer_results %>% 
      filter(!is.na(date_group)) %>% 
      # mutate(coloring = ifelse(date_group == '2021-04-25', 'bad', 'good')) %>% 
      ggplot(aes(date_group, new_inf_svi_slope_mean)) + 
      geom_col(aes(x = admission_dates, y = n*.31), 
               data = age_zip_inf_df %>% 
                 unnest(admission_dates) %>% 
                 count(admission_dates) %>% 
                 mutate(n = slider::slide_dbl(n, mean, .before = 3, .after = 3)), 
               inherit.aes = F, alpha = .2, width = 1) +
      geom_point() +
      geom_errorbar(aes(ymin = new_inf_svi_slope_lo, 
                        ymax = new_inf_svi_slope_hi))+
      # scale_color_manual(values = c('purple', 'black')) +
      labs(y = 'Relative infection risk') -> p
  } else{
    lmer_results %>% 
      filter(!is.na(date_group)) %>% 
      # mutate(coloring = ifelse(date_group == '2021-04-25', 'bad', 'good')) %>% 
      ggplot(aes(date_group, new_cases_svi_slope_mean)) + 
      geom_col(aes(x = admission_dates, y = n*.033), 
               data = age_zip_inf_df %>% 
                 unnest(admission_dates) %>% 
                 count(admission_dates) %>% 
                 mutate(n = slider::slide_dbl(n, mean, .before = 3, .after = 3)), 
               inherit.aes = F, alpha = .2, width = 1) +
      geom_point() +
      geom_errorbar(aes(ymin = new_cases_svi_slope_lo, 
                        ymax = new_cases_svi_slope_hi)) +
      # scale_color_manual(values = c('purple', 'black')) +
      # coord_cartesian(ylim = c(-.1, .02)) +
      labs(y = 'Relative reporting rate') -> p
  }
  # browser()
  p + 
    labs(x = NULL) + 
    geom_hline(yintercept = 1, color = "darkred", size = 1, lty = 2) +
    scale_x_date(date_breaks = '2 month', date_labels = '%b-%y',
                 limits = c(ymd(first_infection_date)-days(14), ymd(final_infection_date))) + 
    theme(legend.position = 'none') +
    background_grid(major ='xy')
}

## 4 Week slope of the line plots
infrate_4week_plot <- make_4week_plot(lmer_4week_results %>% 
                                        filter(!is.na(date_group)), 'inf')
infrate_4week_plot 


rr_4week_plot <- make_4week_plot(lmer_4week_results, 'rr')
rr_4week_plot



plot_grid(cum_inf_corr_plot, infrate_4week_plot,# infrate_inequality_scatter_plot,
          rr_corr_plot, rr_4week_plot,#rr_inequality_scatter_plot,
          nrow = 2, labels = 'AUTO', rel_widths = c(1, 1.5, 1, 1.5), align = 'h')-> fig4_svi_inf_report_plot
fig4_svi_inf_report_plot

save_plot('figs/fig4_svi_inf_report_plot.png', fig4_svi_inf_report_plot, base_height = 7, base_asp = 1.8,bg = "white")
save_plot('figs/fig4.tiff', fig4_svi_inf_report_plot, base_height = 7, base_asp = 1.8,bg = "white")



# SVI relationship stats for ms -------------------------------------------
## Infection rate impact of SVI parameter
inf_glmer_mod %>% sjPlot::tab_model()


## Reporting rate impact of SVI parameter
rr_glmer_mod %>% sjPlot::tab_model()

## Overall SVI grouping difference in infection risk
exp(c(coef(inf_glmer_mod)$zip[1,'svi'], inf_svi_glmer_ci)*lmer_max_diff)

## Overall reporting rate results
exp(c(coef(rr_glmer_mod)$zip[1,'svi'], rr_svi_glmer_ci)*lmer_max_diff)

## example time points for infection risk
# lmer_4week_results


save.image(file='processed-data/fig4-image.rda')


# Potential plot for supplement -------------------------------------------

# zip_4week_infections_cases_df %>%
#   left_join(zip_cov_df %>% select(zip,population,svi),by = 'zip') %>% ungroup() %>%
#   group_by(date_group, zip, svi) %>%
#   summarize(new_inf_rate = mean(new_inf/population),
#             new_inf_rate_lo = quantile(new_inf/population, probs = 0.025),
#             new_inf_rate_hi = quantile(new_inf/population, probs = 0.975),
#             new_inf = mean(new_inf),
#             population = unique(population)) %>% # group_by(date_group) %>% summarize(sum(new_inf))
#   ggplot(aes(svi, new_inf_rate)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = new_inf_rate_lo, ymax = new_inf_rate_hi)) +
#   facet_wrap(~date_group, scales = 'free_y')
