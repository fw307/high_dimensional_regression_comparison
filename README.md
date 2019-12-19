##### 19/12/19: This repository will be updated over the next few days.
***

Code to produce results in "High-dimensional regression in practice: an empirical study of finite-sample prediction, variable selection and ranking" by Wang et al. 

Overview
R script to generate simulated data and perform linear regression using Lasso, Light elastic net, Heavy elastic net, Ridge regression, Dantzig Selector, SCAD and Stability Selection, and generate the following metrics (see paper for details):

Partial Area Under the Curve (pAUC), Root Mean Square Error (RMSE), True Positive Rate (TPR) and Positive Predictive Values (PPV).


User first chooses a setting from settings listed in "SETTINGS.txt" (or provide a setting in the same format which is applicable). Pass the setting to "data_generation" function to generate training and test data accordingly. Using the data generated, one can apply different penalised regression methods to obtain the metrics for ranking, prediction and seleciton, using the script "simulation.R".

The results in "results.7z" zip file are obtained as described above, named after the metrics and the settings (specifically, in the "metric_method_setting.Rdata" format); each individual file contains 64 metric values of a method in a simulation setting, over 64 replications

"plotting_code_main_simulation" folder contains scripts for generating Figures 1-10 for systematic simulation section in the paper, named after their figure numbers. 

Simulation studies were performed using R version 3.0.2. 

