This repository contains customized Matlab codes to reproduce the six figures in the paper, Liao et al., Compilation of longitudinal microbiota data and hospitalome from hematopoietic cell transplantation patients, Scientific Data (2021). The folders are structured as the following:

Figure 1, displaying a patient timeline: example1_display_patient_timeline/main.m

Figure 2, visualizing the entire dataset of microbiota compositions: example2_visualize_compositional_states/main.m

Figure 3, relative frequency of antibiotics administered to HCT patients at MSK by oral and intravenous routes: example3_drug_administration_stats/main.m

Figure 4, inferring the impact of antibiotics on the microbiota using regularized regression: example4_impacts_of_antibiotics/main.m

Figure 5,6, intestinal domination increases a patient’s risk of bloodstream infection: example5_survival_analysis/main.m

To run our codes, you do not need to download any data other than the github folder. A copy of the "hospitalome" data, which is also available in Figshare (https://figshare.com/collections/Compilation_of_longitudinal_microbiota_data_and_hospitalome_from_hematopoietic_cell_transplantation_patients/5271128), was pre-downloaded in the folder "deidentified_data_tables". The only difference of the local database compared to the online version is that the local version has two taxonomy tables that classify amplicon sequence variants (ASVs) using Silva 132 (tblASVtaxonomy_silva132_v4v5_filter.csv) and 138 (tblASVtaxonomy_silva138_v4v5_filter.csv). The online version only includes tblASVtaxonomy_silva132_v4v5_filter.csv, which is the default file used to reproduce all figures in our Scientific Data paper.
