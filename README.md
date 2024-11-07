# MS_turtle-abund_orthomosaics

This repository contains data and code scripts for the manuscript *"Estimating abundance of aggregated populations with drones while accounting for multiple sources of errors: a case study on the mass nesting of Giant South American River Turtles"* (Brack et al. *not published*).

## Repository content

-   `Appendix2_tutorial/`: files used to create Appendix S2, with the tutorial to simulate and analyze data

-   `data/`: original turtle data

    -   `encounter-history_counts.csv`: number of appearances for each individual in each occasion
    -   `encounter-history_states.csv`: observed state for each individual in each occasion
    -   `mark_occasions.csv`: marking occasion for each individual
    -   `mark_identification.txt`: number of marked individuals with marks identified and unidentified
    -   `PopulationCounts.txt`: Overall population counts

-   `ms/`: manuscript files used in the submission

-   `outputs/`: output objects from models fitted to the turtle data

    -   `MR5_postMCMC.txt`: posterior samples of the mark-resight model (step 1)
    -   `Counts5_phi1_postMCMC/`: folder containing the posterior samples for the population counts model (step 2) fitted for each posterior sample of step 1

-   `R/`: R code scripts to fit the turtle data

    -   `nimble_counts b(T)theta(.)phi(.)delta(.)omega(.).R`: nimble code for the population counts model (step 2)
    -   `nimble_MR5 theta(.)phi(occ)delta(.)omega(.).R`: nimble code for the mark-resight model (step 1)
    -   `see_results multi-samps.R`: code to see model results and create figures
    -   `turtle_fit step1 MR5.R`: code to fit the mark-resight model to the turtle data (step 1)
    -   `turtle_fit step2 counts5.R`: code to fit the population counts model to the turtle data (step 2)
