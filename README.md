# EpiFusion Publication Repository
This is a repository to fully reproduce the results of the paper 'EpiFusion: Joint inference of the effective reproduction number by integrating phylodynamic and epidemiological modelling with particle filtering'. By cloning this repository and following the instructions in this README and in the `replicate_results.Rmd` file, you will be able to fully reproduce the results in the paper, including the appendix.

## Models in this paper and where they appear
### EpiFusion Models
1. **Simulated Baseline/Introduction Scenario, Combined:**
* The data for this model is pictured in Fig 2
* Table 2
* Figure 3
* Figure 4
* Figure 5
* Figure 6
* Table 4
2. **Simulated Baseline/Introduction Scenario, Case Incidence Only:** 
* The data for this model is pictured in Fig 2
* Table 2
* Figure 3
* Figure 4
* Figure 5
* Figure 6
* Table 4
3. **Simulated Baseline/Introduction Scenario, Phylogenetic Tree Only:** 
* The data for this model is pictured in Fig 2
* Table 2
* Figure 3
* Figure 4
* Figure 5
* Figure 6
* Table 4
4. **Simulated Sampling Step-Change Scenario, Combined:** 
* The data for this model is pictured in Supplementary Figure S1
* Table 2
* Figure 3
* Figure 4
* Figure 5
5. **Simulated Sampling Step-Change Scenario, Case Incidence Only:** 
* The data for this model is pictured in Supplementary Figure S1
* Table 2
* Figure 3
* Figure 4
* Figure 5
6. **Simulated Sampling Step-Change Scenario, Phylogenetic Tree Only:** 
* The data for this model is pictured in Supplementary Figure S1
* Table 2
* Figure 3
* Figure 4
* Figure 5
7. **Simulated Transmission Step-Change Scenario, Combined:** 
* The data for this model is pictured in Supplementary Figure S1
* Table 2
* Figure 3
* Figure 4
* Figure 5
8. **Simulated Transmission Step-Change Scenario, Case Incidence Only:** 
* The data for this model is pictured in Supplementary Figure S1
* Table 2
* Figure 3
* Figure 4
* Figure 5
9. **Simulated Transmission Step-Change Scenario, Phylogenetic Tree Only:** 
* The data for this model is pictured in Supplementary Figure S1
* Table 2
* Figure 3
* Figure 4
* Figure 5

## Rerunning whole models vs using paper results
### Paper results
The output results of the models run for this paper are in the `Results` folder. Each EpiFusion model has it's own sub-folder, and the benchmarking results (other programs) are in the subfolder `benchmarking`. If you work through the `replicate_results.Rmd` file in this repository, you will be able to fully recreate the graphs and tables in this paper using these output files

### Rerunning whole models
Should you wish to rerun the models quoted in the paper using the EpiFusion program, the XML files of the models are in `Data/`, and the EpiFusion jar file used to run them is included in the repository `EpiFusion.jar`. If you have cloned this repository and are in the working directory, you can use the below commands verbatim to rerun the models. The model XML files are preloaded with file names that will place them in a new folder, called `Results2`, so as not to overwrite the actual paper results. If you use this option of rerunning the models, be sure to change the results repository to `Results2` in the `replicate_results.Rmd` file when you are plotting the results. While EpiFusion runs quite quickly, there are X models run for the paper (with instructions on how to rerun each below), and this may take some time.

## Recreating all plots and table statistics

