# High-dimensional regression comparison
This repository contains R scripts and data files for the large-scale simulation study presented in:  
[Wang, Mukherjee, Richardson & Hill (2019) High-dimensional regression in practice: an empirical study of finite-sample prediction, variable selection and ranking. *Statistics and Computing*. Advance Online Publication. doi:10.1007/s11222-019-09914-9](https://doi.org/10.1007/s11222-019-09914-9).  
R scripts to generate data, run penalized regression methods, evaluate method performance and visualize results are provided, together with performance metric data.

### Simulation study overview
The main simulation study in Wang et al. (2019) included simulations based on synthetic data (simulated covariates and responses) and semisynthetic data (real covariates and simulated responses). Various data-generating factors were explored, including sample size, dimensionality, sparsity, signal strength and multicollinearity, resulting in over 2300 data-generating scenarios. Penalized linear regression methods included in the study were: Lasso, Adaptive Lasso, Elastic Net, Ridge Regression, SCAD, the Dantzig Selector and Stability Selection. Three goals were considered: variable ranking, variable selection and prediction. Metrics used to assess performance were: partial area under the ROC curve (pAUC) for ranking, true positive rate (TPR) and positive predictive value (PPV) for selection, and root-mean-squared error (RMSE) for prediction.  
See the manuscript for full details.

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
Results (performance metric data) from the simulation study are provided in the [results](results) directory as RData files. For each simulation scenario, results were averaged across 64 simulated datasets. Both the non-averaged and averaged results are provided.

### Visualizing simulation results and reproducing manuscript figures
Functions to produce plots of the type presented in the manuscript are provided in the [plottingCode](plottingCode) directory in the R file [plottingFunctions.R](plottingFunctions.R). Each function is documented within this file and loads in results data from the [results](results) directory. These plotting functions can be used to further explore the results data.

All the manuscript figures pertaining to the main simulation study (including supplementary figures) can be reproduced as follows:
```
source("plottingCode/manuscriptFigsScript.R")
```
Figures are saved in the [manuscriptFigures](manuscriptFigures) directory.

##### Required installs:
```
install.packages('ggplot2', 'GGally', 'latex2exp', 'grid', 
                    'gridExtra',' tidyr', 'dplyr', 'RColorBrewer')
```
