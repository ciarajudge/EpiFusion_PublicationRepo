# Simulation Based Calibration

To further validate and examine model performance we characterised its ability to recapture epidemic trajectories and epidemiological parameters. To do this we randomly generated 500 unique epidemics and matching datasets with `beta`, `gamma`, `psi` and `phi` parameters drawn from distributions. These distributions were then used as the corresponding priors for EpiFusion models that tried to reconstruct the truth of the epidemics from the data, with absolutely no manual optimisation. We ran excess replicates, to allow the discarding of any non-convergent runs or runs with very small ESS, with the rationale that in a real scenario, where the truth was unknown, non-convergence or low ESS would indicate a poor model.

We have not added the raw model results to GitHub as the file sizes are exceedingly large, but the final output of this section (a selection of tables with summary statistics for each rep, trajectory fits, calibration tables etc) is in the folder and is used in the results plotting. However, if you want to recreated the interim raw model results (and simulated outbreaks), the code for this is in `simulation_based_calibration.R` and `simulation_based_calibration_results.R`. More detail on this code is below:

## Code
The code that automated this process is in `simulation_based_calibration.R`, as an R script rather than Rmd to make automation slightly easier. The steps taken in this script are walked through below. The automation process makes some assumptions: 

1. That you have java installed and it is accessible in the command line by typing 'java'
2. That you have `BEAST 2.7.6` and it is in `/Applications/BEAST\ 2.7.3/bin/beast`
3. That you have `ReMASTER` installed

### Data Simulation (`simulation_based_calibration.R`)
N draws from distributions for the values of `beta`, `gamma`, `psi` and `phi` are drawn, where N is the number of replicates. These draws are saved in a CSV file and used to automatically create a ReMASTER parameter file for each rep. The files associated with each rep are saved in folders automatically created in the format `rep_n`. ReMASTER is run, and the tree and incidence output is parsed and inserted into a template EpiFusion xml file for each rep.

### Running EpiFusion (`simulation_based_calibration.R`)
The EpiFusion param file for each rep is run with EpiFusion.

### Parsing EpiFusion Output (`simulation_based_calibration_results.R`)
The EpiFusion outputs, some information on timing, the corresponding true values, and other useful things for each rep are saved in a number of tables that are used for the final results. These are the tables currently in the repository, so if you rerun `simulation_based_calibration.R` and `simulation_based_calibration_results.R`, they will be overwritten.








