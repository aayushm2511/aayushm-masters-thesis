# aayushm-masters-thesis
Code for my master's thesis -  dealing with response diversity and ecosystem stability in the Western Ghats

Code files:

1.) imbalance_analysis.R - code for GAM fitting and imabalance calculations (ref. Polazzo et al., 2025) as per Section 2.2 of my thesis

2.) gam_evaluation_metrics - code to obtain model evaluation metrics, including but not limited to AUC and Deviance Explained

3.) sensitivity_analysis - code to re-run imbalance calculations using only species with significant CWD terms in the GAM (88/183).

4.) CWD_bin_geography_check - code to check within-group spatial spread of CWD groups (ref. Section 2.1).

5.) Temporal_variability.R - code to estimate temporal variability in Enhanced Vegetation Index following White et al., 2020, pertaining to section 2.3 of the thesis

6.) thesis_analysis.R - code to conduct statistical analysis amongs the variables of interest (CWD, imbalance, richness, variability).

In addition, the following files were used as primary data sources:

1.) geb13350-sup-0002-tree-species-occ.csv - field data taken from Dr. Meghna Krishnadas (ref. Krishnadas et al., 2021). 

2.) cleaned_evi.csv - EVI data from the MYD13Q1 product (Didan et al., 2021), extracted using AppEARS (https://appeears.earthdatacloud.nasa.gov/). Raw EVI data was cleaned by removing redundant/unneccesary columns
