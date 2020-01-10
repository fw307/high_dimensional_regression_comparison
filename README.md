# High-dimensional regression comparison
This repository contains R scripts and data files for the large-scale simulation study presented in:  
[Wang, Mukherjee, Richardson & Hill (2019) High-dimensional regression in practice: an empirical study of finite-sample prediction, variable selection and ranking. *Statistics and Computing*. Advance Online Publication. doi:10.1007/s11222-019-09914-9](https://doi.org/10.1007/s11222-019-09914-9).  
In particular, R scripts to generate data, run penalized regression methods, evaluate method performance and visualize results are provided, together with performance metric data.

### Simulation study overview
The main simulation study in Wang et al. (2019) included simulations based on synthetic data (simulated covariates and responses) and semisynthetic data (real covariates and simulated responses). Various data-generating factors were explored, including sample size, dimensionality, sparsity, signal strength and multicollinearity, resulting in over 2300 data-generating scenarios. Penalized linear regression methods included in the study were: Lasso, Adaptive Lasso, Elastic Net, Ridge Regression, SCAD, the Dantzig Selector and Stability Selection. Three goals were considered: variable ranking, variable selection and prediction. Metrics used to assess performance were: partial area under the ROC curve (pAUC) for ranking, true positive rate (TPR) and positive predictive value (PPV) for selection, and root-mean-squared error (RMSE) for prediction.  
See the paper for full details.

### Simulation scripts
The files [`syntheticDataSimulation.R`](syntheticDataSimulation.R) and [`semisyntheticDataSimulation.R`](semisyntheticDataSimulation.R) contain code to generate data, run the penalized regression methods and evaluate performance for synthetic data and semisynthetic data scenarios, respectively.

The scenarios included in the simulation study are contained in the text files [`syntheticDataScenarios.txt`](syntheticDataScenarios.txt) and [`semisyntheticDataScenarios.txt`](semisyntheticDataScenarios.txt). These are tables specifying the factors that define each scenario (rows are scenarios and columns are factors). The simulation scripts read in the corresponding scenario table (the scenario text file needs to be in the working directory) and a specific scenario (row) can then be selected by the user to use in a simulation. A new scenario, not included in the simulation study, could also be specified by the user.

Note that the synthetic data scenario file [`syntheticDataScenarios.txt`](syntheticDataScenarios.txt) also includes the Toeplitz correlation design simulation scenarios from the "Additional investigations" section of the paper (Section 4.1). Toeplitz design scenarios can also be simulated using the script [`syntheticDataSimulation.R`](syntheticDataSimulation.R).

##### Required installs:
```
install.packages('MASS', 'glmnet', 'parcor', 'ncvreg', 'monomvn'.
                    'flare', 'c060', 'Matrix', 'caret')
```

### Simulation results
The results in "results.7z" zip file are obtained as described above, named after the metrics and the settings (specifically, in the "metric_method_setting.Rdata" format); each individual file contains 64 metric values of a method in a simulation setting, over 64 replications

"plotting_code_main_simulation" folder contains scripts for generating Figures 1-10 for systematic simulation section in the paper, named after their figure numbers. 



