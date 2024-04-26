# bmi_genomics

This repository contains the code to reproduce the findings and figures presented in [Tang, et al., Obesity shapes selection for driver mutations in cancer, medRxiv (2024).](https://www.medrxiv.org/content/10.1101/2024.01.10.24301114v2)

Before running any scripts, please download the necessary data from Zenodo. The generated data structure is as follows:
  - g2use - list of genes tested
  - clin - clinical features of each patient
  - mutsns - binary matrix of all non-synonymous mutations
  - mutsonc -  binary matrix of all oncogenic mutations
  - mutssilent - binary matrix of all silent mutations
  - bmi_dfci - raw counts of KRAS and EGFR mutations in validation cohort
  - gaddy - BMI information to make supplemental figure 1
  - wl_lung - cachexia labels


In order to generate the main results reported in the paper, please run the script "BMI Analysis and Figures Script.R". This will run the logistic regression models and produce figure 1, supplemental figure 1, and supplemental figure 2.

