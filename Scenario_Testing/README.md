# Scenario Testing

For the scenario testing section we simulated three outbreak scenarios: a baseline scenario, a scenario with a step-change in sampling, and a scenario with a step-change in transmission.

The code to generate these datasets, as well as the raw simulated data, is in `Data_Simulation`. The code is in the markdown file `Data_Simulation.Rmd`.

The EpiFusion XMLs (input files for EpiFusion) with the data and parameters are in `EpiFusion_XMLs`.

The EpiFusion output folders used to generate the results of this paper are in `Results`.

`Results2` is an empty folder that the EpiFusion XMLs in `EpiFusion_XMLs` have been programmed to write to if they are rerun. This is to allow you to rerun the XMLs without overwriting the paper results in `Results`.


