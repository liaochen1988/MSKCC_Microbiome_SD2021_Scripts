This repository contains customized Matlab codes to reproduce the six figures in the paper (https://www.nature.com/articles/s41597-021-00860-8), Liao et al., Compilation of longitudinal microbiota data and hospitalome from hematopoietic cell transplantation patients, Scientific Data (2021). The folders are structured as the following:

Figure 1, displaying a patient timeline: example1_display_patient_timeline/main.m

Figure 2, visualizing the entire dataset of microbiota compositions: example2_visualize_compositional_states/main.m

Figure 3, relative frequency of antibiotics administered to HCT patients at MSK by oral and intravenous routes: example3_drug_administration_stats/main.m

Figure 4, inferring the impact of antibiotics on the microbiota using regularized regression: example4_impacts_of_antibiotics/main.m

Figure 5,6, intestinal domination increases a patient’s risk of bloodstream infection: example5_survival_analysis/main.m

To run our codes, you need to create a new folder "deidentified_data_tables" and download the "hospitalome" data from Figshare (https://figshare.com/collections/Compilation_of_longitudinal_microbiota_data_and_hospitalome_from_hematopoietic_cell_transplantation_patients/5271128) into the folder. The "deidentified_data_tables" folder should be structured as follows:

Subfolder counts: tblcounts_asv_melt.csv, tblcounts_asv_wide.csv, tblcounts_class_wide.csv, tblccounts_family_wide.csv, tblcounts_genus_wide.csv, tblcounts_order_wide.csv, tblcounts_phylum_wide.csv, tblqpcr.csv

Subfolder meta_data: tblbc.csv, tbldrug.csv, tblhctmeta.csv, tblInfectionsCidPapers.csv, tbltemperature.csv, tblVanA.csv

Subfolder samples: tblASVsamples.csv

Subfolder taxonomy: tblASVtaxonomy_silva132_v4v5_filter.csv
