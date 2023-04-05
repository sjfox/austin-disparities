library(tidyverse)
library(lubridate)
library(readxl)
library(tidycensus)
library(cowplot)
library(sf)
theme_set(theme_cowplot())
source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')

set.seed(23123)

nadmits <- 150
ex_shape1 <- 25
ex_shape2 <- 100

tibble(zip = 'example1',
       age_group_paper = 'example1',
       shape1 = ex_shape1, ## shape confirmation means about 1 in 5 infections are hospitalized
       shape2 = ex_shape2,
       n_admits = nadmits,
       population = 10000) %>% 
  get_infection_estimates(num_samps = 1000) -> ex1_infections
hosp_date_options <- seq.Date(from = ymd('2020-04-01'), 
                              to = ymd('2020-10-01'), by = 'day')
hosp_dates <- sample(hosp_date_options, 
                     prob = dnorm(1:length(hosp_date_options),  ## Makes the distribution of hospital admission normal (epidemic lookking)
                                  mean = length(hosp_date_options)/2, 
                                  sd = length(hosp_date_options)/5),
                     replace = T,
                     size = nadmits)
ex1_infections %>% 
  unnest(est_infections) %>% 
  mutate(inf_dates = map(est_infections, 
                         get_infection_timing, 
                         zipage_hosp_vec = hosp_dates))  -> ex1_timing
  
mindate <- ymd('2020-02-01')
maxdate <- ymd('2020-10-01')

tibble(ihr = rbeta(100000, ex_shape1, ex_shape2)) %>% 
  ggplot(aes(ihr)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    labs(x = 'Infection hospitalization rate',
         y = 'Probability density') +
    background_grid(major = 'xy') +
    scale_y_continuous(expand = c(0,0)) -> ihr_explot
ihr_explot

ex1_infections %>% 
  unnest(est_infections) %>% 
  ggplot(aes(est_infections)) +
  geom_histogram(bins = 50) +
  labs(x = 'Estimated infections',
       y = 'Sample count') +
  background_grid(major = 'xy') +
  scale_y_continuous(expand = c(0,0)) -> inf_explot
inf_explot


tibble(date = hosp_dates) %>% 
  count(date) %>% 
  ggplot(aes(date, n)) +
  geom_col() +
  scale_x_date(limits = c(mindate,maxdate)) +
  background_grid(major = 'xy') +
  labs(x = NULL, y = 'Hospital admissions') ->hosp_timing_explot
hosp_timing_explot

ex1_timing %>% 
  mutate(id = seq_along(zip)) %>% 
  unnest(inf_dates) %>% 
  group_by(id) %>% 
  mutate(cuminf = cumsum(new_infections)) %>% 
  ggplot(aes(date, cuminf, group = id)) + 
    geom_line(alpha = .1) +
    scale_x_date(limits = c(mindate,maxdate)) +
    scale_y_continuous(labels = scales::comma) +
    background_grid(major = 'xy') +
    labs(x = NULL, y = 'Cumulative estimated infections') -> inf_timing_explot
inf_timing_explot


plot_grid(plot_grid(ihr_explot,
                    inf_explot, 
                    nrow = 2,
                    labels = c('A', 'C'),
                    align='v'),
          plot_grid(hosp_timing_explot,
                    inf_timing_explot, 
                    nrow = 2,
                    labels = c('B', 'D'),
                    align='v'),
          labels = NA,
          rel_widths = c(1,2)
) -> ex_estimation_fig
ex_estimation_fig

save_plot('figs/ex_estimation_fig.png',
          ex_estimation_fig,
          base_height = 7,
          base_asp = 1.6,
          bg='white')
          



