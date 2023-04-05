library(R2jags)
library(rjags)
## Creates sampler for a single zip code
## Have age-based hospitalization and IHR data
## only a single data point for case data and reporting rate

samp_conditional_reporting_rate <- function(detected_cases, 
                                            estimated_infections,
                                            num_samps){
  if(any(estimated_infections - detected_cases < 0)){
    stop('Error: More detected cases than infections.')
  }
  rbeta(num_samps,
        1 + detected_cases, 
        1 + estimated_infections - detected_cases)
}



get_infection_posterior <- function(shape1, 
                                    shape2, 
                                    n_admits, 
                                    population, 
                                    num_samps){
  model_code <- "
    model
    {
      # Likelihood
      H ~ dbinom(mu, round(I))
      
      # Priors
      mu ~ dbeta(mu_shape1, mu_shape2)
      I ~ dunif(minI, maxI)
    }
    "
  
  # Set up the data
  model_data <- list(H = n_admits, minI = n_admits, maxI = population, 
                     mu_shape1 = shape1, mu_shape2 = shape2)
  
  # Choose the parameters to watch
  model_parameters <- c("I")
  # Run the model
  model_run <- jags(
    data = model_data,
    parameters.to.save = model_parameters,
    model.file = textConnection(model_code),
    n.chains = 4, # Number of different starting positions
    n.iter = (num_samps)*2/4+200, # Number of iterations
    n.burnin = 200, # Number of iterations to remove at start
    n.thin = 2
  )
  
  round(model_run$BUGSoutput$sims.list$I[,1])
}

get_infection_estimates <- function(age_zip_df,
                                    num_samps = 1000){
  age_zip_df %>% 
    # head() %>% 
    mutate(est_infections = pmap(.l = tibble(shape1 = shape1,
                                             shape2 = shape2,
                                             n_admits = n_admits,
                                             population = population,
                                             num_samps = num_samps), 
                                 .f = get_infection_posterior) ) %>%
    select(zip, age_group_paper, est_infections) 
    # nest(est_infections = c(est_infections)) 
}

get_infection_timing <- function(total_inf, 
                                 zipage_hosp_vec,
                                 hosp_delay_shape = 11.05^2/40.85,
                                 hosp_delay_rate = 11.05/40.85){
  tibble(infection_date = sample(zipage_hosp_vec, size = total_inf, replace = T) -
                            days((rgamma(total_inf, 
                                        shape = hosp_delay_shape, 
                                        rate = hosp_delay_rate) %>% 
                                   round()))) %>% 
    count(infection_date) %>% 
    rename(date = infection_date,
           new_infections = n)  
}


# estimate_infections_and_timing <- function(num_hospital_admits,
#                                            population_size,
#                                            ihr_shape1,
#                                            ihr_shape2,
#                                            hosp_delay_shape = 11.05^2/40.85,
#                                            hosp_delay_rate = 11.05/40.85,
#                                            hosp_df = NULL,
#                                            num_samps = 1000
#                                            ){
#   
#   # browser()
#   get_infection_estimates(num_hospitalizations = num_hospital_admits, 
#                           ihr_shape1 = ihr_shape1, 
#                           ihr_shape2 = ihr_shape2, 
#                           population_size = population_size, 
#                           num_samps = num_samps) %>% 
#     mutate(samp = seq_along(ihr)) -> infection_df
#   
#   if(sum(num_hospital_admits) > 0){
#     infection_df %>% 
#       mutate(infection_times = map(.x = infections, 
#                                    get_infection_timing,
#                                    hosp_df = hosp_df,
#                                    hosp_delay_shape = hosp_delay_shape,
#                                    hosp_delay_rate = hosp_delay_rate)) -> infection_df
#   } 
#   infection_df
# }

