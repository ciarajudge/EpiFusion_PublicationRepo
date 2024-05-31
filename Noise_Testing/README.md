# Noise Testing

For the noise testing section we simulated outbreaks with increasing observation and transmission noise (10 noise levels x 2 types of noise = 20 outbreaks), and examined the performance of EpiFusion phylo-only, epi-only and combined models on them. The interim results (raw simulated data and raw EpiFusion outputs) are not on the repo due to space but using the code below it's possible to reproduce them. The EpiFusion XMLs with the prepped data are hosted here on the repo, however, The summary of the analysis in csv form (`noisesummarytable.csv`) is what we use to make the plots in the report.

## Data Simulation
The code used to generate the datasets for the analysis is in `noise_testing_generate_datasets.R`. The automation process makes some assumptions: 

1. That you have java installed and it is accessible in the command line by typing 'java'
2. That you have `BEAST 2.7.6` and it is in `/Applications/BEAST\ 2.7.6/bin/beast`
3. That you have `ReMASTER` installed

The script generates datasets with increasing noise, and creates EpiFusion XMLs which are places in the `EpiFusion_XMLs` folder. Three XMLs are created per dataset that carry out combined, epi-only and phylo-only analyses.


## Parsing the results
The script `full_noise_results.Rmd` was used to parse the results of the EpiFusion runs, which are programmed by their XML files to place the results in the `Results` folder. The raw results are not uploaded to the github, but the folder is still there if you decide to rerun the XML files in `EpiFusion_XMLs`. The R markdown script inspects the runs, plots them, and creates the `noisesummarytable.csv` file which we use in the results plotting for the manuscript.
