## File processing functions
require(readxl)
require(tidycensus)

# Fitting the beta distribution -------------------------------------------
# Get mean and upper/lower bounds
get_upper_lower <- function(shape1, shape2){
  # browser()
  bounds <- qbeta(p = c(0.025, .975), shape1 = shape1, shape2 = shape2) 
  c(shape1/(shape1+shape2), bounds[1], bounds[2]) 
}

ll_beta <- function(pars, ihr_vec){ # function to minimize to get best param estimates
  ## This fits only the second shape parameter of beta distribution
  ## The first shape parameter is constrained so that 
  ## the mean of distribution matches the mean IHR
  # shape1 <- as.numeric((2 * ihr_vec[1] - ihr_vec[1] * exp(pars[1]) - 1) / (ihr_vec[1] - 1)) # constrain to mode
  shape1 <- as.numeric((ihr_vec[1] * exp(pars[1])) / (1 - ihr_vec[1])) # constrain to mean
  samp_ihr_vec <- get_upper_lower(shape1, exp(pars[1]))
  
  sum((samp_ihr_vec[2] - ihr_vec[2])^2) ## Focus on lower bound
  # sum((samp_ihr_vec[2:3] - ihr_vec[2:3])^2) ## Focus on both bounds
  # sum((samp_ihr_vec[3] - ihr_vec[3])^2) ## Focus only on upper bound
}

get_beta_parms <- function(ihr_mean, ihr_lb, ihr_ub){
  ihr_vec <- c(ihr_mean, ihr_lb, ihr_ub)
  if(any(is.na(ihr_vec))){stop('There are NAs in your IHR')}
  parm_est <- optimize(ll_beta, c(0.000, 25), ihr_vec = ihr_vec)
  shape2 = exp(parm_est$minimum)
  shape1 = as.numeric(ihr_mean*shape2/(1-ihr_mean))
  quantiles = qbeta(c(0.025,0.975), shape1, shape2)
  
  tibble(shape1 = shape1,
         shape2 = shape2,
         mean_fit = shape1 / (shape1 +shape2),
         lower_fit = quantiles[1],
         upper_fit = quantiles[2])
}

get_travis_ihr_age_table <- function(){
  # read_csv('produced_data/ihr_files/COVID_IHR_per_age_group_per_zip_rho6_base_country_France.csv') %>% 
  read_csv('processed-data/zip-age-ihr_estimated.csv') %>% 
    inner_join(read_csv('raw-data/spatial-files/zips_in_austin_rr_msa.csv') %>% 
                 filter(COUNTY == 'Travis'), 
               by = c('ZCTA5' = 'ZIP')) %>% 
    select(zip = ZCTA5,
           metric = Measure, `0_0.5`:`75+`) %>% 
    gather(age_group, value, `0_0.5`:`75+`) %>% 
    spread(metric, value) %>% 
    mutate(ihr_age_group = age_group,
           age_group = ifelse(age_group == '75+', '75_120', age_group))  %>% 
    separate(age_group, into = c('start_age', 'end_age'), sep = "_") %>% 
    mutate(end_age = ifelse(is.na(end_age), 120, end_age),
           end_age = ifelse(end_age == '0.5', 0, end_age),
           start_age = ifelse(start_age == '0.5', 1, start_age)) %>% 
    gather(key, value, start_age:end_age) %>% 
    group_by(ihr_age_group, zip, LowerBound, Mean, UpperBound) %>% 
    arrange(value) %>% 
    mutate(value = as.numeric(value)) %>% 
    padr::pad_int(by = 'value', group = c('ihr_age_group','zip', 'LowerBound', 'Mean', 'UpperBound')) %>% 
    distinct(ihr_age_group, value, zip, LowerBound, Mean, UpperBound) %>% 
    ungroup() %>% 
    mutate(zip = as.character(zip)) %>% 
    select(age = value,
           ihr_age_group,
           zip,
           ihrlower = LowerBound,
           ihrmean = Mean,
           ihrupper = UpperBound) 
}


get_za_population = function(refresh=FALSE){
  if(!dir.exists('processed-data/population-data')){
    dir.create('processed-data/population-data')
  }
  ## Gets the zip and age-grouped population estimates from ACS 2019
  if(!file.exists("processed-data/population-data/travis_10yr_pop_per_zip_master.csv") | 
     refresh){
    tibble(acs_variable_code = sprintf("B01001_%0.3d", c(3:25, 27:49)), ## 3:25 males, 27:49 females
           age_grouping = rep(c('0-4', '5-9', '10-14', '15-17', 
                                '18-19', '20', '21','22-24',
                                '25-29', '30-34', '35-39', '40-44', 
                                '45-49', '50-54', '55-59', '60-61', 
                                '62-64', '65-66', '67-69', '70-74', 
                                '75-79', '80-84', '85+'), 2)) -> acs_vars
    
    age_dict <- tibble(acs_age_group = c('0-4', '5-9', '10-14', '15-17', 
                                         '18-19', '20', '21','22-24',
                                         '25-29', '30-34', '35-39', '40-44', 
                                         '45-49', '50-54', '55-59', '60-61', 
                                         '62-64', '65-66', '67-69', '70-74', 
                                         '75-79', '80-84', '85+'),
                       age_group = c('0-9', '0-9', '10-19', '10-19', 
                                     '10-19', '20-29', '20-29', '20-29',
                                     '20-29', '30-39', '30-39', '40-49',
                                     '40-49', '50-59', '50-59', '60-69',
                                     '60-69', '60-69', '60-69', '70-79',
                                     '70-79', '80+', '80+'))
    
    #data_year = 2019 # default is 2018 for the 2014-2018 ACS database
    pop = get_acs(geography="zcta", variables= acs_vars$acs_variable_code, geometry=FALSE, year = 2019)
    
    pop %>% 
      mutate(GEOID = str_sub(GEOID, start=-5)) %>% 
      inner_join(read_csv('raw-data/spatial-files/zips_in_austin_rr_msa.csv') %>% 
                   mutate(ZIP = as.character(ZIP)) %>% 
                   filter(COUNTY == 'Travis'), by = c('GEOID' = 'ZIP')) %>% 
      left_join(acs_vars, by = c('variable' = 'acs_variable_code')) %>% 
      left_join(age_dict, by = c('age_grouping' = 'acs_age_group')) %>% 
      select(age_group,
             zip = GEOID,
             pop=estimate) %>% 
      group_by(age_group, zip) %>% 
      summarize(population = sum(pop)) %>% 
      ungroup() -> travis_acs_pop_estimates10yr
    
    write_csv(travis_acs_pop_estimates10yr, 
              file = "processed-data/population-data/travis_10yr_pop_per_zip_master.csv")
  } else{ # read in file if it already exists as this takes some time to create
    travis_acs_pop_estimates10yr = read_csv("processed-data/population-data/travis_10yr_pop_per_zip_master.csv", 
                                            col_types = c(zip="c"))
  }
  travis_acs_pop_estimates10yr
} 


get_age_zip_df <- function(patient_hosp_data, 
                           refresh_acs=FALSE){
  # browser()
  ## Need to produce tibble of all the zip and age group combinations, alongside their populations, IHRs, and admission data
  
  age_hosp_dict <- tibble(age = 0:120,
                          age_group10yr = c(rep(c('0-9',  '10-19', '20-29', '30-39', 
                                                  '40-49','50-59', '60-69', '70-79'), each=10),
                                            rep('80+', 41)),
                          age_group_paper = c(rep('0-17', 18), 
                                              rep('18-49', 32), 
                                              rep('50-64', 15),
                                              rep('65+', 56)))
  overlapping_pop_weights <- tibble(age_group10yr = c('10-19', '10-19', '60-69', '60-69'),
                                    age_group_paper = c('0-17', '18-49', '50-64', '65+'),
                                    pop_weight = c(8/10, 2/10, 5/10, 5/10))
  
  ## Process the age data, starts in 10 year increments, and need to summarize by age groups for paper
  get_za_population(refresh = refresh_acs) %>% 
    left_join(age_hosp_dict %>% 
                distinct(age_group10yr, age_group_paper), 
              by = c('age_group' = 'age_group10yr')) %>% 
    left_join(overlapping_pop_weights, 
              by = c('age_group' = 'age_group10yr', 'age_group_paper')) %>% 
    mutate(pop_weight = ifelse(is.na(pop_weight), 1, pop_weight)) %>% 
    group_by(zip, age_group_paper) %>% 
    summarize(population = round(sum(population*pop_weight))) %>% 
    ungroup() -> za_population
  
  ihr_za_table <- get_travis_ihr_age_table() 
  
  ## Create zip/age ihr table that will be used for zipcodes and ages 
  ## Can't average by population, because it doesn't recreate infection trends
  ## Young groups get too much weight, and IHRs are pulled down
  ## This averages ihrs across the ages according to the total hospital admissions by age
  ihr_za_table %>% 
    left_join(patient_hosp_data %>% 
                count(age), by = 'age') %>% 
    filter(!is.na(n)) %>% 
    inner_join(age_hosp_dict, by = 'age') %>% 
    group_by(zip,age_group_paper) %>% 
    summarize_at(.vars = vars(ihrlower:ihrupper), 
                 .funs = ~weighted.mean(x = ., w = n)) %>% ##This averages age groups into final paper age groups based on admission weighting
    ungroup() -> ihr_agezip
  
  correct_date_nulls <- function(dates, curr_zip, patient_hosp_data){
    if(is.null(dates)){
      patient_hosp_data %>% 
        filter(zip == curr_zip) %>% 
        pull(admit_date)
    }else{
      dates
    }
  }
  
  ihr_agezip %>% 
    inner_join(za_population, by = c('zip','age_group_paper')) %>% 
    left_join(patient_hosp_data %>% 
                left_join(age_hosp_dict) %>% 
                mutate(zip = as.character(zip)) %>% 
                group_by(zip, age_group_paper) %>% 
                summarize(n_admits = n(),
                          admission_dates = list(admit_date)), by = c('zip', 'age_group_paper')) %>% 
    mutate(n_admits = ifelse(is.na(n_admits), 0, n_admits),
           admission_dates = map2(admission_dates, zip, 
                                  correct_date_nulls, 
                                  patient_hosp_data = patient_hosp_data)) -> age_zip_df
  
  age_zip_df %>% 
    bind_cols(age_zip_df %>% 
                select(ihr_mean = ihrmean,
                       ihr_lb = ihrlower,
                       ihr_ub = ihrupper) %>% 
                pmap(get_beta_parms) %>% 
                bind_rows())
}


get_age_df <- function(age_zip_df, cutoff_date){
  ## Want to produce tibble of age population, and reported cases for easy use later
  # browser()
  age_zip_df %>% 
    group_by(age_group_paper) %>% 
    summarize(population = sum(population)) %>% 
    left_join(    read_csv('raw-data/Confirmed_Demographics_(Public_View).csv') %>% 
                    select_at(vars(Last_Updat, starts_with('Age_Ttl'))) %>% 
                    rename(date = Last_Updat) %>%
                    mutate(`0-17` = Age_Ttl_Le + Age_Ttl_01 + round(Age_Ttl_10*8/10),
                           `18-49` = round(Age_Ttl_10*2/10) + Age_Ttl_20 + Age_Ttl_30 + Age_Ttl_40,
                           `50-64` = Age_Ttl_50 + round(Age_Ttl_60*1/2),
                           `65+` = round(Age_Ttl_60*1/2) + Age_Ttl_70 + Age_Ttl_80) %>% 
                    select(date, `0-17`:`65+`) %>% 
                    gather(age_group_paper, cum_cases, -date) %>% 
                    mutate(date = as.Date(ymd_hms(date))) %>% 
                    filter(date < cutoff_date, !is.na(cum_cases)) %>% 
                    group_by(age_group_paper) %>% 
                    summarize(cases = max(cum_cases, na.rm=T),
                              tibble(date = date, cum_cases = cum_cases) %>% nest(data = everything())), 
                  by = 'age_group_paper')
}


get_SVI_per_ZIP = function(){
  zips_arr_msa = read_csv("raw-data/spatial-files/zips_in_austin_rr_msa.csv", 
                          col_types = cols(ZIP = col_character())) %>%
    filter(!(ZIP=="76573"))
  
  #### DATA USED FOR SVI CONVERSION ####
  # convert between zip code and census tract as of 09-2020 # https://www.huduser.gov/portal/datasets/usps_crosswalk.html#codebook
  ctract_zip_crosswalk = read_csv("raw-data/spatial-files/ZIP_TRACT_092020.csv") 
  msa_crosswalk = subset(ctract_zip_crosswalk, ZIP %in% zips_arr_msa$ZIP)
  ctract_svi_tx = read_csv("raw-data/spatial-files/CDC_SVI_Texas_2018.csv") %>% # CDC SVI metrics for the state of Texas
    rename(TRACT = FIPS) # rename for join
  msa_crosswalk_svi = msa_crosswalk %>% 
    left_join(ctract_svi_tx, by="TRACT") %>% # join SVI data to the ZIP to TRACT conversion table
    select(ZIP, TRACT, RES_RATIO, RPL_THEMES, everything())
  
  # Calculate the weighted sum of SVI and it's percentage components for each ZIP with residential ratio in each TRACT
  msa_crosswalk_svi %>%
    group_by(ZIP) %>%
    filter(RPL_THEMES>=0) %>%
    summarise(ZIP_COVERED = sum(RES_RATIO),      # proportion of ZIP used in residential areas of Tract
              SVI = sum(RES_RATIO*RPL_THEMES),   # SVI
              POV = sum(RES_RATIO*EP_POV),       # percent below poverty line
              UNEMP = sum(RES_RATIO*EP_UNEMP),   # unemployment rate estimate
              PCI = sum(RES_RATIO*EP_PCI),       # per capita income
              NOHSDP = sum(RES_RATIO*EP_NOHSDP), # no high school diploma
              AGE65 = sum(RES_RATIO*EP_AGE65),   # percent over age 65
              AGE17 = sum(RES_RATIO*EP_AGE17),   # percent 17 and below
              DISABL = sum(RES_RATIO*EP_DISABL), # percent non-institutionalized disabled
              SNGPNT = sum(RES_RATIO*EP_SNGPNT), # percent single parent house hold children under 18
              MINRTY = sum(RES_RATIO*EP_MINRTY), # percent all minorities, excludes white+non-hispanic
              LIMENG = sum(RES_RATIO*EP_LIMENG), # percent greater than 5 with limited english
              MUNIT = sum(RES_RATIO*EP_MUNIT),   # percent housing structures with 10 or more units
              MOBILE = sum(RES_RATIO*EP_MOBILE), # percent mobile homes
              CROWD = sum(RES_RATIO*EP_CROWD),   # percent occupied houses with more people than rooms
              NOVEH = sum(RES_RATIO*EP_NOVEH),   # percent households no vehicle
              GROUPQ = sum(RES_RATIO*EP_GROUPQ), # percent in institutionalized group quarters
              RPL_THEME1 = sum(RES_RATIO*RPL_THEME1),
              RPL_THEME2 = sum(RES_RATIO*RPL_THEME2),
              RPL_THEME3 = sum(RES_RATIO*RPL_THEME3),
              RPL_THEME4 = sum(RES_RATIO*RPL_THEME4)
    ) 
} 


get_zip_df <- function(age_zip_df, cutoff_date){
  ## Want to produce tibble of zip, population, SVI, and reported cases for easy use later
  age_zip_df %>% 
    group_by(zip) %>% 
    summarize(population = sum(population)) %>% 
    left_join(read_csv('raw-data/Austin_Travis_Zip_Code_Counts_(Public_View).csv') %>% 
                filter(is.na(Notes)) %>% 
                select(-Notes, -FID) %>% 
                mutate(Date = as.Date(ymd_hms(Date))) %>% 
                gather(zip, cum_cases, -Date) %>% 
                separate(col = 'zip', into = c('remove', 'zip')) %>% 
                select(date = Date, zip, cum_cases) %>%
                arrange(date) %>% 
                filter(date < cutoff_date, !is.na(cum_cases)) %>% 
                group_by(zip) %>% 
                summarize(cases = max(cum_cases, na.rm=T),
                          tibble(date = date, cum_cases = cum_cases) %>% nest(data = everything())) %>% 
                inner_join(read_csv('raw-data/spatial-files/zips_in_austin_rr_msa.csv') %>% 
                             mutate(ZIP = as.character(ZIP)) %>% 
                             filter(COUNTY == 'Travis') %>% 
                             select(ZIP), by = c('zip' = 'ZIP')) %>% 
                rename(ts_cases = data), by = 'zip') %>% 
    left_join(get_SVI_per_ZIP() %>% 
                select(zip = ZIP,
                       svi = SVI) %>% 
                mutate(zip = as.character(zip)), by = 'zip')
}








