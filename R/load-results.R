


load('processed-data/results_upto_01-06-2021_2022-run.rda')

if(file.exists('processed-data/shared-parms.rda')){
  load('processed-data/shared-parms.rda')
} else{
  # Parameters used throughout ----------------------------------------------
  first_infection_date <- '2020-03-01'
  final_infection_date <- '2021-06-01'
  
  # Overall metropolitan region data and stats ------------------------------
  travis_pop <- age_cov_df %>% pull(population) %>% sum()
  ## Cases and deaths
  travis_nyt <- read_csv('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv') %>% 
    filter(county == 'Travis', state == 'Texas')
  travis_nyt %>% filter(date == final_infection_date)
  travis_total_cases <- travis_nyt %>% filter(date == final_infection_date) %>% pull(cases)
  
  ## Reported hospitalizations
  age_zip_inf_df %>% 
    summarize(total_admits = sum(n_admits))
  
  ## Estimate for total infections
  age_zip_inf_df %>% 
    unnest(inf_timing) %>% 
    mutate(zip_age = paste0(zip,'_',age_group_paper)) %>% 
    select(id, zip_age, date, new_infections) %>% 
    spread(zip_age, new_infections) %>% 
    mutate(new_inf = rowSums(.[c(-1,-2)], na.rm = T)) %>% 
    select(id, date, new_inf) %>% 
    group_by(id) %>% 
    mutate(cum_inf = cumsum(new_inf)) -> all_inf_ts
  
  if(file.exists('processed-data/real-aph-tests.csv')){
    aph_zip_age_tests <- read_csv('processed-data/real-aph-tests.csv', col_types = c('ccnn'))
  } else{
    aph_zip_age_tests <- read_csv('raw-data/fake-aph-tests.csv', col_types = c('ccnn'))
  }
  
  save.image(file='processed-data/shared-parms.rda')  
}

