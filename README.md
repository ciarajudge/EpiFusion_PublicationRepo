# EpiFusion Publication Repository
This repository houses the data and code to fully reproduce the results of the paper 'EpiFusion: Joint inference of the effective reproduction number by integrating phylodynamic and epidemiological modelling with particle filtering'. By cloning this repository and following the instructions in this README and in the `replicate_results.Rmd` file, you will be able to fully reproduce the results in the paper, including the appendix.

## Main Sections of the Paper
### 1. Likelihood Comparison
To validate our discretised phylogenetic likelihood we compared the likelihood of ranges of parameters given the same dataset under a Phylo-Only EpiFusion model to a BDSky model. We repeated this on a selection of 7 total datasets with different true parameters for `beta`, `gamma` and `psi`, trying 3 different true values of each. This required the automated creation and running of 450 EpiFusion parameter files. Scripts to recreate this section, as well as further information, is in the `Likelihood_Comparison` folder.

### 2. Simulation Based Calibration
We also carried out a simulation based calibration of the model, where we simulated 500 unique epidemics and matching datasets with `beta`, `gamma`, `psi` and `phi` parameters drawn from distributions. These distributions were then used as the corresponding priors for EpiFusion models that tried to reconstruct the truth of the epidemics from the data. We ran excess replicates, to allow the discarding of any non-convergent runs or runs with very small ESS, with the rationale that in a real scenario, where the truth was unknown, non-convergence or low ESS would indicate a poor model. Scripts to recreate this section, as well as further information, is in the `Simulation_Based_Calibration` folder. The summary results of this section are housed in this section also in CSV files.

### 3. Scenario Testing 
This is the main section of the paper, where we tested EpiFusion on simulated data generated from a range of theoretical epidemiological scenarios, and compare the performance of the Combined EpiFusion approach against Phylo-Only and Epi-Only EpiFusion runs. This is the section we most recommend trying to reproduce from the code, as the number of EpiFusion models run in this section in the smallest (compared to Simulation Based Calibration, for example, where we run 500 models in a loop). The steps in reproducing the results of this section are below:

#### 3.1 Simulating Data
The information, XML files, and parsing code used to create each simulated dataset is in the `Scenario_Testing/Data_Simulation` folder, divided out by section and by simulation. The `Scenario_Testing/Data_Simulation` folder also contains a README with an index of every simulated dataset used in the Scenario Testing.

#### 3.2 Running EpiFusion Models (EpiFusion XMLs)
The `Scenario_Testing/EpiFusion_XMLs` folder contains the final EpiFusion XML parameter files (including any data simulated in step 1) for every model featured in the Scenario Testing section. Should you wish to rerun the models using the EpiFusion program, you can use these XML files, and the EpiFusion jar file used to run them `EpiFusion.jar`. The model XML files are preloaded with file names that will place them in a new folder, called `Scenario_Testing/Results2`, so as not to overwrite the actual paper results which are in `Scenario_Testing/Results`.

### 4. Noise Testing
To test model robustness to noise we simulated outbreaks with increasing levels of observation and transmission noise (10 levels, 20 outbreaks in total) and compared EpiFusion Epi-Only, Phylo-Only and Combined models on their ability to infer R_t. The code to create the data, form EpiFusion XML files and run the EpiFusion analyses is in `Noise_Testing`.

### 5. Benchmarking (against other programs)
This paper also features models run with other tools in the section where we benchmark EpiFusion performance. We do not include detailed information on how to rerun these other programs, but we do provide the parameter files in the folder `Benchmarking/Data`, and the results are already in the `Benchmarking/Results` folder.


## Reproducing plots and tables
This repository has all the results featured in the paper in the folders for their specific section, and you can use the .Rmd file `figures_and_tables.Rmd` to reproduce plots and tables. If you have rerun the scenario models (Step 1.2) and want to use your own output, be sure to change the `scenario_results` folder in the script to `Scenario_Testing/Results2` (this will be clearly marked in the script).


