library(tidyverse)
library(lubridate)
library(readxl)
library(tidycensus)
library(cowplot)
library(sf)
theme_set(theme_cowplot())
source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')

select <- dplyr::select

# Load all data for analysis ---------------------------------------
if(file.exists('processed-data/real_hosp_admits.csv')){
  patient_hosp_data <- read_csv('processed-data/real_hosp_admits.csv')
} else{
  patient_hosp_data <- read_csv('raw-data/fake_hosp_admits.csv')
}

data_cutoff_date <- ymd('2021-06-01')
age_zip_df <- get_age_zip_df(patient_hosp_data = patient_hosp_data, refresh_acs = F)
age_cov_df <- get_age_df(age_zip_df, data_cutoff_date)
zip_cov_df <- get_zip_df(age_zip_df, data_cutoff_date)


age_zip_df %>% 
  left_join(age_zip_df %>% 
              get_infection_estimates(num_samps = 1000), 
            by = c('age_group_paper', 'zip'))   %>% 
  unnest(est_infections) %>% 
  mutate(inf_timing = pmap(tibble(total_inf = est_infections,
                                  zipage_hosp_vec = admission_dates), 
                           get_infection_timing))  %>% 
  group_by_at(vars(-est_infections, -inf_timing)) %>% 
  summarize(est_infections = list(est_infections),
            bind_rows(inf_timing, .id = 'id') %>% nest(data = everything())) %>% 
  rename(inf_timing = data) %>% 
  ungroup() -> age_zip_inf_df

save(age_zip_df, age_cov_df, 
     zip_cov_df, age_zip_inf_df, 
     file = 'processed-data/results_upto_01-06-2021_2022-run.rda')

