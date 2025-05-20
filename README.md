## Supplemantary Code to manuscript ``Statistical-Analysis-of-Large-scale-Item-Response-Data-under-Measurement-Non-invariance``


### Abstract 

With the process of globalization and advancements in technology, International Large-scale Assessments in education (ILSAs), such as the Programme for International Student Assessment (PISA) and the Trends in International Mathematics and Science Study (TIMSS), have become increasingly important in educational research and policy-making. They collect valuable data on education quality and performance development across many education systems worldwide, allowing countries to share techniques, organizational structures, and policies that have proven efficient and successful. 
A key to analyzing ILSA data is an Item Response Theory (IRT) model, which is used to estimate the performance distributions of different groups (e.g., countries) and then produce a ranking. 
A major challenge in calibrating the IRT model is that some items suffer from Differential Item Functioning (DIF), i.e., different groups have different probabilities of correctly answering the items after controlling for individual proficiency levels. DIF is particularly common in ILSA due to the differences in test languages, cultural contexts, and curriculum designs across different groups. Ignoring or improperly accounting for DIF when calibrating the IRT model can lead to severe biases in the estimated performance distributions, which may further distort the ranking of the groups. Unfortunately,  existing methods cannot guarantee the statistically consistent recovery of the group ranking without unrealistic assumptions for ILSA, such as the existence and knowledge of reference groups and anchor items. To fill this gap, this paper proposes a new approach to DIF analysis across multiple groups. This approach is computationally efficient and statistically consistent, does not require assumptions about reference groups and anchor items, 
and provides uncertainty quantification for the group-specific parameters. 
The proposed method is applied to PISA 2022 data to analyze data from the mathematics, science, and reading domains, providing insights into their DIF structures and the performance rankings of countries. 

### Simulation 

#### Proposed Method
- `simulation-code/main_proposed_S1.R` to `simulation-code/main_proposed_S4.R`: Main scripts to run the proposed method under settings S1–S4.
- `simulation-code/functions.R`: Source file containing core functions for the proposed method.

Run `simulation-code/main_proposed_S1.R` to `simulation-code/main_proposed_S4.R` with random seed 80-89 (each seed corresponds to 10 independent replications) and sample size N = 10000 and 20000 produces `simulation-output/proposed_S1` to `simulation-output/proposed_S4`, respectively.

#### RMSD Method
- `simulation-code/main_RMSD_S1.R` to `simulation-code/main_RMSD_S4.R`: Main scripts to run the RMSD method under settings S1–S4.
- `simulation-code/functions_rmsd.R`: Source file containing core functions for the RMSD method.

Run `simulation-code/main_RMSD_S1.R` to `simulation-code/main_RMSD_S4.R` with random seed 80, 82, 84, 86, 88 (each seed corresponds to 20 independent replications) and sample size N = 10000 and 20000 produces `simulation-output/rmsd_S1` to `simulation-output/rmsd_S4`, respectively.

#### Baseline Method
- `simulation-code/main_baseline_S1.R` to `simulation-code/main_baseline_S4.R`: Main scripts to run the baseline method under settings S1–S4.
- `simulation-code/baseline_functions.R`: Source file containing core functions for the baseline method.

Run `simulation-code/main_proposed_S1.R` to `simulation-code/main_proposed_S4.R` with random seed 80-89 (each seed corresponds to 10 independent replications) and sample size N = 10000 and 20000 produces `simulation-output/proposed_S1` to `simulation-output/proposed_S4`, respectively.

#### Evaluation

Run `simulation-results-analysis/boxplot compare S1.R` to `simulation-results-analysis/boxplot compare S4.R` with N = 10000 and N = 20000 produces 8 plots in Figure 2 and Table 4.

Run `simulation-results-analysis/import_data_S1_S3.R` and then run `mu_gamma_MSE.R` and `a_d_sigma_MSE.R` to obtain S1 and S3 of Tables 2-3.

Run `simulation-results-analysis/import_data_S2_S4.R` and then run `mu_gamma_MSE.R` and `a_d_sigma_MSE.R` to obtain S2 and S4 of Tables 2-3.

