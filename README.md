# BayesPharmacogenetics

R data and code for the paper:
"A Bayesian Approach for Investigating the Pharmacogenetics of Combination Antiretroviral Therapy in People with HIV".

## Instructions for Use

This repository contains the R code generating a simulated dataset that has similar design as the real data in the paper "A Bayesian Approach for Investigating the Pharmacogenetics of Combination Antiretroviral Therapy in People with HIV" and the R code implementing the proposed method.

In the folder "Data_and_Code":

* The R data file "Treatment_History_Data.Rdata" contains the treatment history data for n=200 individuals randomly sampled from the Women's Interagency HIV Study (WIHS) dataset, which will be used to generated simulated dataset in the R script "Data_Generate.R"; 

* The R script “MCMC_R_Functions.R” provides R functions used for MCMC, the Rcpp script “MCMC_Rcpp_Functions.cpp” provides Rcpp functions used for MCMC, and the R script "MCMC_Update_Tree" provides R functions for tree-related MCMC steps;

* The R script "MCMC_Main.R" runs the MCMC algorithm for the proposed method.  
