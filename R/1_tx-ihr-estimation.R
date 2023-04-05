## Texas IHRs
library(tidyverse)
library(lubridate)
library(cowplot)
library(reticulate)
use_condaenv("r-reticulate", required=T)
library(writexl)
theme_set(theme_cowplot())
set.seed(235)
get_tx_pop <- function(){
  require(tidycensus)
  ## Get population estimates
  tibble(acs_variable_code = sprintf("B01001_%0.3d", c(3:25, 27:49)), ## 3:25 males, 27:49 females
         age_grouping = rep(c('0-4', '5-9', '10-14', '15-17', 
                              '18-19', '20', '21','22-24',
                              '25-29', '30-34', '35-39', '40-44', 
                              '45-49', '50-54', '55-59', '60-61', 
                              '62-64', '65-66', '67-69', '70-74', 
                              '75-79', '80-84', '85+'), 2)) -> acs_vars
  
  age_dict = tibble(acs_age_group = c('0-4', '5-9', '10-14', '15-17', 
                                      '18-19', '20', '21','22-24',
                                      '25-29', '30-34', '35-39', '40-44', 
                                      '45-49', '50-54', '55-59', '60-61', 
                                      '62-64', '65-66', '67-69', '70-74', 
                                      '75-79', '80-84', '85+'),
                    age_group = c('0-17', '0-17', '0-17', '0-17', 
                                  '18-49', '18-49', '18-49', '18-49',
                                  '18-49', '18-49', '18-49', '18-49',
                                  '18-49', '50-64', '50-64', '50-64',
                                  '50-64', '65+', '65+', '65+',
                                  '65+', '65+', '65+'))
  
  #data_year = 2019 # default is 2018 for the 2014-2018 ACS database
  pop = get_acs(geography="state", variables= acs_vars$acs_variable_code, geometry=FALSE, year = 2019)
  
  # browser()
  pop %>% 
    left_join(acs_vars, by = c('variable' = 'acs_variable_code')) %>% 
    left_join(age_dict, by = c('age_grouping' = 'acs_age_group')) %>%
    filter(NAME == 'Texas') %>% 
    group_by(age_group) %>% 
    summarize(pop = sum(estimate))
}

if(!dir.exists('processed-data/population-data/')){
  dir.create('processed-data/population-data/')
}

if(!file.exists('processed-data/population-data/tx_pop_age.csv')){
  tx_pop <- get_tx_pop()  
  write_csv(tx_pop, 'processed-data/population-data/tx_pop_age.csv')
} else{
  tx_pop <- read_csv('processed-data/population-data/tx_pop_age.csv')
}




## Correspond to infections 7 days before initial specimen collection
## https://jamanetwork.com/journals/jamainternalmedicine/fullarticle/2768834
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

#https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/g62h-syeh
## Read in TX hospital admission data
read_csv('raw-data/COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv') %>% 
  select(state, 
         date, 
         `<18` = previous_day_admission_pediatric_covid_confirmed,
         `18-19` = `previous_day_admission_adult_covid_confirmed_18-19`,
         `20-29` = `previous_day_admission_adult_covid_confirmed_20-29`,
         `30-39` = `previous_day_admission_adult_covid_confirmed_30-39`,
         `40-49` = `previous_day_admission_adult_covid_confirmed_40-49`,
         `50-59` = `previous_day_admission_adult_covid_confirmed_50-59`,
         `60-69` = `previous_day_admission_adult_covid_confirmed_60-69`,
         `70-79` = `previous_day_admission_adult_covid_confirmed_70-79`,
         `80+` = `previous_day_admission_adult_covid_confirmed_80+`) %>% 
  filter(state == 'TX', date > '2020-07-23') %>% ##Removes noisy first days of dataset
  arrange(date) %>% 
  gather(age_group, n, -state, -date) %>% 
  mutate(n = ifelse(age_group == '18-19' & ##These lines clean up issues with hospitalizations in this age group
                      (date == '2020-08-15' |
                         date == '2020-09-14' |
                         date == '2020-09-15' |
                         date == '2020-09-16' 
                      ), 0, n)) %>% 
  mutate(n = ifelse(n<0, 0, n)) %>% 
  mutate(age_group = factor(age_group, 
                            levels = c('<18', '18-19', '20-29', 
                                       '30-39', '40-49', '50-59', 
                                       '60-69', '70-79', '80+'))) -> tx_admits


# Calculate hosp and in from Aug2020 to Jun2021 ---------------------------
num_samps <- 1000
get_age_hosp_admits_delay <- function(hosp_delay, admits) {
  admits %>% 
    filter(date >= ymd('2020-07-29') + days(hosp_delay),
           date < ymd('2021-05-27') + days(hosp_delay)) %>% 
    group_by(age_group) %>% 
    summarize(n_admit = sum(n)) %>% 
    mutate(age_grouping = list('0-17', '18-49', '18-49', '18-49', '18-49', '50-64', c('50-64','65+'), '65+', '65+')) %>% 
    mutate(n_admit = ifelse(age_group == '60-69',round(n_admit/2), n_admit)) %>% 
    unnest(age_grouping) %>% 
    group_by(age_grouping) %>% 
    summarize(n_admit = sum(n_admit)) %>% 
    rename(age_group = age_grouping) 
}
get_norm_sd <- function(mean, lwr, upr) {
  ## Conservatively takes the maximum standard deviation
  ## 1.96 is critical value to get standard deviation
  ## exponentiate to get variance
  (pmax((mean-lwr), (upr-mean))/1.96)
}


hosp_delay_shape = 11.05^2/40.85
hosp_delay_rate = 11.05/40.85
hosp_delays <- round(rgamma(1000, shape = hosp_delay_shape, rate = hosp_delay_rate))

hosp_delays %>% 
  map(get_age_hosp_admits_delay, admits = tx_admits) %>% 
  bind_rows(.id = 'id') -> admit_samps

## Get distributions for generating infection estimates
tx_seroprevalence %>% 
  filter(begin_date == '2021-05-27') %>% 
  left_join(tx_seroprevalence %>% 
              filter(begin_date == '2020-07-29') %>% 
              select(-state, -begin_date, -round, -pop) %>% 
              rename(baseline_lower = lower,
                     baseline_mean = mean,
                     baseline_upper = upper)
  ) %>% 
  mutate(baseline_sd = get_norm_sd(baseline_mean, baseline_lower, baseline_upper),
         final_sd = get_norm_sd(mean, lower, upper)) %>% 
  select(age_group, 
         final_mean=mean, final_sd = final_sd,
         baseline_mean, baseline_sd) -> sero_est_inf_distributions

## Now sample infections for each admission sample
admit_samps %>% 
  left_join(sero_est_inf_distributions, by = 'age_group') %>% 
  mutate(seroestimated_infected = round(rnorm(nrow(admit_samps), mean = final_mean, sd = final_sd) -
                                    rnorm(nrow(admit_samps), mean = baseline_mean, sd = baseline_sd))) %>% 
  mutate(ihr = n_admit/seroestimated_infected) %>% 
  group_by(age_group) %>% 
  summarize(
    admit_lower = quantile(n_admit, probs = 0.025),
    admit_mean = mean(n_admit),
    admit_upper = quantile(n_admit, probs = 0.975),
    ni_lower = quantile(seroestimated_infected, probs = 0.025),
    ni_mean = mean(seroestimated_infected),
    ni_upper = quantile(seroestimated_infected, probs = 0.975),
    ihr_lower = quantile(ihr, probs = 0.025),
    ihr_mean = mean(ihr),
    ihr_upper = quantile(ihr, probs = 0.975)) -> tx_ihr

tx_ihr %>% 
  mutate(n_admits = paste0(scales::comma(admit_mean), ' (',
                         scales::comma(admit_lower),'-',
                         scales::comma(admit_upper),')'
            ),
         new_inf = paste0(scales::comma(ni_mean), ' (',
                          scales::comma(ni_lower),'-',
                          scales::comma(ni_upper),')'
                          ),
         ihr = paste0(scales::percent(ihr_mean, accuracy = .01), ' (',
                      scales::percent(ihr_lower, accuracy = .01),'-',
                      scales::percent(ihr_upper, accuracy = .01),')'
                      )
         ) %>% 
  select(age_group, n_admits, new_inf, ihr) -> tx_ihr_printed

tx_ihr_printed
write_csv(tx_ihr_printed, 'processed-data/tx_ihr_printed.csv')



# Add IHRs to spreadsheet -------------------------------------------------

ihr_comp <- read_csv('raw-data/ihr-comparison.csv')

ihr_comp %>% 
  filter(study != 'fox') %>% 
  bind_rows(tx_ihr %>% 
              mutate(age_lo =c(0,18,50,65),
                     age_hi =c(17,49,64,100),
                     mean = ihr_mean*100,
                     lo = ihr_lower*100,
                     hi = ihr_upper*100,
                     study = 'fox', 
                     region = 'Texas') %>% 
              select(age_lo, age_hi, mean, lo, hi, study, region)) ->new_ihr_comp

write_csv(new_ihr_comp, 'raw-data/ihr-comparison.csv')

  

# Run IHR script with the IHR estimates to get ZIP IHRS ----------------------------
zip_age_groupings <- c("0_0.5", "0.5_4", "5_9", "10_14", 
                       "15_19", "15_19", "15_19", "15_19", "15_19", 
                       "20_24", "25_29", "30_34", "35_39", 
                       "40_44", "45_49", "50_54", "55_59", 
                       "60_64", "65_69", "70_74", "75+")
zip_age_groupings <- factor(zip_age_groupings, levels = unique(zip_age_groupings))

tx_ihr %>% 
  select(AgeGroup = age_group, 
         Mean = ihr_mean,
         LowerBound = ihr_lower,
         UpperBound = ihr_upper) %>% 
  slice(c(1,1,1,1,
          1,1,1,2,2,
          2,2,2,2,2,2,
          3,3,3,
          4,4,4)) %>% 
  mutate(AgeGroup = zip_age_groupings) %>% 
  group_by(AgeGroup) %>% 
  summarize_all(mean) -> ihr_for_zip_conversion

writexl::write_xlsx(list("IHR" = ihr_for_zip_conversion),
                    'ihr-code/estimated-tx-ihr.xlsx')


# Run python script with reticulate ---------------------------------------

setwd('ihr-code/')
source_python('high_risk_covid_to_IHR_age_specific_rho.py')
setwd('../')

file.rename(from = 'ihr-code/zip-age-ihr_estimated.csv', 
            to = 'processed-data/zip-age-ihr_estimated.csv')
file.remove('ihr-code/estimated-tx-ihr.xlsx')
