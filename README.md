# Statistical-Analysis-of-Large-scale-Item-Response-Data-under-Measurement-Non-invariance

## Abstract 

With the process of globalization and advancements in technology, International Large-scale Assessments in education (ILSAs), such as the Programme for International Student Assessment (PISA) and the Trends in International Mathematics and Science Study (TIMSS), have become increasingly important in educational research and policy-making. They collect valuable data on education quality and performance development across many education systems worldwide, allowing countries to share techniques, organizational structures, and policies that have proven efficient and successful. 
A key to analyzing ILSA data is an Item Response Theory (IRT) model, which is used to estimate the performance distributions of different groups (e.g., countries) and then produce a ranking. 
A major challenge in calibrating the IRT model is that some items suffer from Differential Item Functioning (DIF), i.e., different groups have different probabilities of correctly answering the items after controlling for individual proficiency levels. DIF is particularly common in ILSA due to the differences in test languages, cultural contexts, and curriculum designs across different groups. Ignoring or improperly accounting for DIF when calibrating the IRT model can lead to severe biases in the estimated performance distributions, which may further distort the ranking of the groups. Unfortunately,  existing methods cannot guarantee the statistically consistent recovery of the group ranking without unrealistic assumptions for ILSA, such as the existence and knowledge of reference groups and anchor items. To fill this gap, this paper proposes a new approach to DIF analysis across multiple groups. This approach is computationally efficient and statistically consistent, does not require assumptions about reference groups and anchor items, 
and provides uncertainty quantification for the group-specific parameters. 
The proposed method is applied to PISA 2022 data to analyze data from the mathematics, science, and reading domains, providing insights into their DIF structures and the performance rankings of countries. 

## Simulation 

#### Proposed Method
- `proposed/main_proposed_S1.R` to `proposed/main_proposed_S4.R`: Main scripts to run the proposed method under settings S1–S4.
- `proposed/functions.R`: Source file containing core functions for the proposed method.

#### RMSD Method
- `RMSD/main_RMSD_S1.R` to `RMSD/main_RMSD_S4.R`: Main scripts to run the RMSD method under settings S1–S4.
- `RMSD/functions_rmsd.R`: Source file containing core functions for the RMSD method.

#### Baseline Method
- `baseline/main_baseline_S1.R` to `baseline/main_baseline_S4.R`: Main scripts to run the baseline method under settings S1–S4.
- `baseline/baseline_functions.R`: Source file containing core functions for the baseline method.



