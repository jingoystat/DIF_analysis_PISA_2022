## Supplemantary Code to manuscript "Statistical Analysis of Large-scale Item Response Data under Measurement Non-invariance: A Statistically Consistent  Method and Its Application to PISA 2022"


### Paper Abstract 

With the process of globalization and advancements in technology, International Large-scale Assessments in education (ILSAs), such as the Programme for International Student Assessment (PISA) and the Trends in International Mathematics and Science Study (TIMSS), have become increasingly important in educational research and policy-making. They collect valuable data on education quality and performance development across many education systems worldwide, allowing countries to share techniques, organizational structures, and policies that have proven efficient and successful. 
A key to analyzing ILSA data is an Item Response Theory (IRT) model, which is used to estimate the performance distributions of different groups (e.g., countries) and then produce a ranking. 
A major challenge in calibrating the IRT model is that some items suffer from Differential Item Functioning (DIF), i.e., different groups have different probabilities of correctly answering the items after controlling for individual proficiency levels. DIF is particularly common in ILSA due to the differences in test languages, cultural contexts, and curriculum designs across different groups. Ignoring or improperly accounting for DIF when calibrating the IRT model can lead to severe biases in the estimated performance distributions, which may further distort the ranking of the groups. Unfortunately,  existing methods cannot guarantee the statistically consistent recovery of the group ranking without unrealistic assumptions for ILSA, such as the existence and knowledge of reference groups and anchor items. To fill this gap, this paper proposes a new approach to DIF analysis across multiple groups. This approach is computationally efficient and statistically consistent, does not require assumptions about reference groups and anchor items, 
and provides uncertainty quantification for the group-specific parameters. 
The proposed method is applied to PISA 2022 data to analyze data from the mathematics, science, and reading domains, providing insights into their DIF structures and the performance rankings of countries. 

---

### Simulation 

#### Proposed Method
- `simulation-code/main_proposed_S1.R` to `simulation-code/main_proposed_S4.R`: Main scripts to run the proposed method under settings S1–S4.
- `simulation-code/functions.R`: Source file containing core functions for the proposed method.

To run the proposed method, execute the scripts `main_proposed_S1.R` to `main_proposed_S4.R` with random seeds 80–89 (each seed corresponds to 10 independent replications) and sample sizes `N = 10000` and `N = 20000`. The results will be saved in `simulation-output/proposed_S1` to `simulation-output/proposed_S4`, respectively.

#### RMSD Method
- `simulation-code/main_RMSD_S1.R` to `simulation-code/main_RMSD_S4.R`: Main scripts to run the RMSD method under settings S1–S4.
- `simulation-code/functions_rmsd.R`: Source file containing core functions for the RMSD method.

To run the RMSD method, execute the scripts `main_RMSD_S1.R` to `main_RMSD_S4.R` with random seeds 80, 82, 84, 86, and 88 (each seed corresponds to 20 independent replications) and sample sizes `N = 10000` and `N = 20000`. The results will be saved in `simulation-output/rmsd_S1` to `simulation-output/rmsd_S4`, respectively.

#### Baseline Method
- `simulation-code/main_baseline_S1.R` to `simulation-code/main_baseline_S4.R`: Main scripts to run the baseline method under settings S1–S4.
- `simulation-code/baseline_functions.R`: Source file containing core functions for the baseline method.

To run the baseline method, execute the scripts `main_baseline_S1.R` to `main_baseline_S4.R` with random seeds 80–89 (each seed corresponds to 10 independent replications) and sample sizes `N = 10000` and `N = 20000`. The results will be saved in `simulation-output/baseline_S1` to `simulation-output/baseline_S4`, respectively.

#### Results Evaluation

- Run `simulation-results-analysis/boxplot compare S1.R` to `boxplot compare S4.R` with `N = 10000` and `N = 20000` to generate eight plots presented in Figure 2 and to generate Table 4.
- Run `simulation-results-analysis/import_data_S1_S3.R`, followed by `mu_gamma_MSE.R` and `a_d_sigma_MSE.R`, to produce results for settings S1 and S3 in Tables 2–3.
- Run `simulation-results-analysis/import_data_S2_S4.R`, followed by `mu_gamma_MSE.R` and `a_d_sigma_MSE.R`, to produce results for settings S2 and S4 in Tables 2–3.

---

### Real Data Application

#### Data Preprocessing

The raw datasets (`CY08MSP_STU_QQQ.SAS7BDAT` and `CY08MSP_STU_COG.SAS7BDAT`) can be downloaded from [OECD PISA 2022 Database](https://www.oecd.org/en/data/datasets/pisa-2022-database.html#data).

To preprocess the data:
- Run `real-data/math_cleaning.R`, `real-data/science_cleaning.R`, and `real-data/reading_cleaning.R`.  
- The cleaned datasets are saved in:
  - `real-data/math_cleaned_data.zip`
  - `real-data/science_cleaned_data.zip`
  - `real-data/reading_cleaned_data.zip`

#### Analysis Results

- Run `real-data-code/estimation_math.R` to obtain estimation results using the proposed method, saved in `real-data-output/math_proposed/`.
- Run `real-data-code/estimation_math_rmsd.R` for results using the RMSD method, saved in `real-data-output/math_rmsd/`.
- Run `real-data-code/estimation_math_baseline.R` for results using the baseline method, saved in `real-data-output/math_baseline/`.


Based on estimation output and cleaned data, run `real-data-analysis/math_results_analysis.R` to produce country ranking in Table 6, histogram in Figure 3, and MDS plot in Figure 4. 

Similar procedure can be applied to analysis on reading and science data.


