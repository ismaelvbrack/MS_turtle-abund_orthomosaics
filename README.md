# MS_turtle-abund_orthomosaics
This repository contains data and code scripts for the manuscript _"Estimating abundance of aggregated populations with drones while accounting for multiple sources of errors: a case study on the mass nesting of Giant South American River Turtles"_ (Brack et al. _not published_).

## Repository content

-   `data/`: original data
    - `mammal_info.txt`: information 
    - `mammal_records.txt`: observed dat
    - `plots_covariates.txt`: site-level 

-   `figs/`: figures generated from the data analysis
    -   `Fig.Lambda-dBurnAWB.png`: 
    -   `Fig.Psi_Coefs_Forest.png`: Slope 
    -   `Fig.Psi-Forest.png`: Relationship 

-   `ms/`: manuscript files used in the submission
  
-   `outputs/`: output objects from models fitted 

-   `R/`: R code scripts
    - `JAGS_models/`: folder containing JAGS models
      - `multisppDA_fatalNmix_covars3_6scales.txt`: selected model: psi(forest) lambda(deltaBurn+tanque) p(.)
      - `multisppDA_fatalNmix_covars3_6scales.txt`: full model: psi(forest+dist2wiw) lambda(deltaBurn+tanque) p(greenVeg)
    - `figures_relationship predictions.R`: code to make the figures
    - `results_coefficient estimates.R`: code to see model results
    - `run_multispp ZIP Nmix carcass.R`: code to import, arrange, and analyze data using JAGS
