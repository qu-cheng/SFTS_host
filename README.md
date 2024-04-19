# Contribution of different host species to the natural transmission of severe fever with thrombocytopenia syndrome virus in China

This repository contains data and code for reproducing the analyses in *Contribution of different host species to the natural transmission of severe fever with thrombocytopenia syndrome virus in China* by Cheng et al.

## Data
This foloder contains data collected from the systematic review and used for developing the mathematical models, including: 
* **StudyLocation.csv**: The names (*Location*), longitute (*Long*), latitude (*Lat*), and reference (*Ref*) for each Survey (*LocationID*).   
* **Species_N_Npos.csv**: Summary data for each species, including  
    *  *Species*, name of the species
    *  *Number of studies*, number of studies surveyed the species of interest
    *  *N*, number of animals tested for the species of interest
    *  *Npos*, number of species tested positive for the species of interest  
* **Species_N_Npos.csv**: Species-specific seroprevalence data by SurveyID, including:
    * *SurveyID*, survey ID
    * *References*, the data source
    * *Survey year*, year the survey was conducted
    * *Locations*, locations included in the survey
    * *Host species*, name of the host species
    * *n_positive*, number of animals tested positive
    * *n*, number of animals tested
    * *prev*, point estimate of the prevalence rate (%)
    * *N*, number of animals in the survey region collected from statistical yearbooks


## Code
This folder contains code for the analyses and figures:
* **0_Library_and_functions.R**: R scripts for loading libraries, running and calibrating the mathematical models
* **1_Descriptive_figures.R**: code for generating descriptive figures from the raw data, including Figure 1 and Figure S2
* **2_model_runs**: this folder contains code and results for calibarting the model, including
    * parameter values used in the mathematical model (*Animal_par.csv*)
    * code for calibrating the model (*2_1_Model_runs_R1.R*, *2_2_Model_runs_R2.R*, and *2_3_Model_runs_R3.R*)
    * code for processing the passes (*2_4_Process_passes.R*), and
    * the corresponding results (*Results* folder), including
       * the 1,000 passing parameter sets (*Passes_1000.csv*)
       * the 10,000 posterior samples obtained from the 1,000 passing parameter sets with the sampling-importance resampling approach (*Passes_1000_resampledID.csv*)
       * summary statistics for the species-level R0i by survey ID *Resampled_species_R0i.csvResampled_species_R0i.csv*)
* **3_Parameter permutation**: this folder contains code and results related to parameter perturbations, including:
    * code for generating perturbed parameter datasets (*3_1_Generate_perturbed_par_sets.R*, the generated dataset will be saved in a folder named *3_2_Perturbed_parameters*, which was not uploaded due to the large file size)
    * code for simulating with the perturbed parameter sets (*3_3_Sim_perturbed_par_sets.R*, the results will be saved in a folder named *3_4_Perturbed_results*, which was not uploaded due to the large file size)
    * code for fitting ranfom forest models to examine the importance of each parameter in determining the species-level R0i (*3_5_Fit_RF.R*)
    * Resulted permutation importance score from the random forest models (*3_6_RF_importance.csv*)
* **4_Optimal_intervention**: this folder contains code (*4_1_Estimate_R0_decrease.R*) and results (*4_2_R0_decrease.csv*) related to examining the effectiveness of each intervention
* **5_Missing_host**: this folder contains code (*5_1_Sim_Missing_Host_R1.R* and *5_2_Sim_Missing_Host_R2.R*) related to examining the impact of missing each individual species being surveyed in Survey 9
* **6_Model_result_figures.R** contains scripts for generating figures related to the model results
