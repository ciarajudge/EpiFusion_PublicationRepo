# Likelihood Comparison

To validate our discretised phylogenetic likelihood we compared the likelihood of ranges of parameters given the same dataset under a Phylo-Only EpiFusion model to a BDSky model. This folder contains all of the code and results for this section.

## Code
The code that automated this process is in `Likelihood_Comparison.R`, as an R script rather than Rmd to make automation slightly easier. The steps taken in this script are walked through below. The automation process makes some assumptions: 

1. That you have java installed and it is accessible in the command line by typing 'java'
2. That you have `BEAST 2.7.6` and it is in `/Applications/BEAST\ 2.7.6/bin/beast`
3. That you have `ReMASTER` installed
4. That you have the `BEAST` package `FEAST` installed

---

## Method
### Dataset generation
We generated 7 total datasets with different true parameters for `beta`, `gamma` and `psi`, trying 3 different true values of each:

| beta | gamma | psi |
| :--- | :--- | :--- |
| 0.19 | 0.12 | 0.03 |
| 0.23 | 0.08 | 0.03 |
| 0.23 | 0.12 | 0.01 |
| 0.23 | 0.12 | 0.03 |
| 0.23 | 0.12 | 0.05 |
| 0.23 | 0.16 | 0.03 |
| 0.27 | 0.12 | 0.03 |


### Likelihood evaluation under different parameters
For the 7 datasets, we carried out 9 tests where two parameters were fixed to the true value, and we examined how the likelihood changed when the third parameter varied around the true value. The discrepancy between 7 datasets and 9 trials is due to the fact that for the `0.23, 0.12, 0.03` dataset, each of the parameters were varied in turn. This led to 9 'tests':

| truebeta | truegamma | truepsi | varyingparam | modelbeta | modelgamma | modelpsi |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0.19 | 0.12 | 0.03 | beta | 0.0475-0.327 | 0.12 | 0.03 |
| 0.23 | 0.08 | 0.03 | gamma | 0.23 | 0.02-0.14 | 0.03 |
| 0.23 | 0.12 | 0.01 | psi | 0.23 | 0.12 | 0.0025-0.0175 |
| 0.23 | 0.12 | 0.03 | beta | 0.0575-0.4025 | 0.12 | 0.03 |
| 0.23 | 0.12 | 0.03 | gamma | 0.23 | 0.03-0.21 | 0.03 |
| 0.23 | 0.12 | 0.03 | psi | 0.23 | 0.12 | 0.0075-0.0525 |
| 0.23 | 0.12 | 0.05 | psi | 0.23 | 0.12 | 0.0125-0.0875 |
| 0.23 | 0.16 | 0.03 | gamma | 0.23 | 0.04-0.28 | 0.03 |
| 0.27 | 0.12 | 0.03 | beta | 0.0675-0.4725 | 0.12 | 0.03 |


This required the creation of 450 EpiFusion parameter files (a range of 50 values around the true value for each of the 9 'tests'). These parameter files were automatically created from the `EpiFusion_template.xml` file, and run with a special version of an EpiFusion jar file, `EpiFusion_SaveAllLikelihoods.jar` that runs 1000 replicates of the density calculation and save each of them.

The same 450 combinations of parameters were evaluated on the data under a BDSky model using `densitymapper` from the `feast` package of `BEAST 2.7.6`.

The results of both sets of likelihood calculations were parsed, yielding a dataframe of 450 rows, and the columns:
`TrueBeta`, `TrueGamma`, `TruePsi`, `ModelBeta`, `ModelGamma`, `ModelPsi`, `AvgLikelihood` (Likelihood according to EpiFusion), `BDSkyDensity` (Likelihood according to BDSky). This data is what is used for the plots in the paper.
