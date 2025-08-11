# CLrenal_pred_22chem

This repository is constructed for the paper "An in vitro-in silico workflow for predicting renal clearance of environmental and pharmaceuticals compounds (TBD)" and contains scripts and data for executing modeling and analysis.

## Contents
1. 22che_mcmc_exe_code.R – R script used to execute MCMC simulations for the 22 chemicals tested using the Transwell™ assay.
2. "MCMC_outputs" folder – Contains the following results: (1) all iterations from four MCMC chains, (2) R-hat diagnostic table, (3) steps file, and (4) simulated concentration data for each chemical.
3. "MCMC_plots" folder – Includes MCMC trace plots, posterior distribution density plots, and comparison plots of measured vs. simulated concentrations for each chemical.
4. 22Che_2D_CL_calculations.xlsx – Spreadsheet detailing the renal CL calculation process for the 22 chemicals tested in the 2D assay.
5. CL_plot_code.R – R script for generating renal clearance plots.
6. "Results" folder - Contains the files used for plotting Figure S2 (P ratio) and Figure 6 (Observed vs Predicted renal clearance).
