## Disproportionate impacts of COVID-19 in a large US city
This repository contains the code that was used to estimate infection rates by age and ZIP code for the first 15 months of the pandemic in the City of Austin. The preprint for the work can be found here: https://www.medrxiv.org/content/10.1101/2022.11.04.22281855v1

## Data
Data were privately provided by Austin Public Health and cannot be shared. They can be made available upon reasonable requests to Austin Public Health. We include fake data so that others may run the analysis pipeline.

## Recreating the analysis
The script `R/00_run-all-analysis.R` can be used to run all analyses and produce all figures. All R packages can be installed at the beginning of the script, and it should produce the R outputs and figures from the manuscript (though not exactly if run with fake data). The main infection estimation code is found in `R/infection-estimation-fxns-jags.R`, and `R/2b_run-example.R` provides a simplified example to showcase each step of the methodology.

## Contact
With questions or comments please contact Spencer Fox via [email](mailto:sjfox@uga.edu) or [twitter](https://twitter.com/foxandtheflu).




