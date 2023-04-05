##### Run analysis and make figures for manuscript

## Uncomment and run to make sure all packages are available
## Also need to install Python so reticulate function works...
# install.packages(c('tidyverse', 
#                    'readxl', 
#                    'writexl',
#                    'reticulate',
#                    'cowplot',
#                    'lubridate',
#                    'sf',
#                    'tidycensus',
#                    'R2jags',
#                    'rjags',
#                    'rmapzen',
#                    'lme4',
#                    'effects',
#                    'sjPlot',
#                    'merTools', 
#                    'padr'))


## Run to install python + needed libraries
# install_python()
# py_install("pandas")
# py_install("numpy")

## Only needed if you have the real data csv
# source('R/0_process-linelist-data.R')

## Gets Texas age-specific IHRs
source('R/1_tx-ihr-estimation.R')

## Runs the full estimation procedure and saves data
source('R/2_run-estimation.R')
source('R/2b_run-example.R')

## Generates figures and results for manuscript
source('R/3_fig1.R')
source('R/3b_sfigs-age-sero.R')
source('R/4_fig2.R')
source('R/5_fig3.R')
source('R/6_fig4.R')
source('R/6b_sfigs-svi-components.R')
source('R/7_supplemental-figs.R')

