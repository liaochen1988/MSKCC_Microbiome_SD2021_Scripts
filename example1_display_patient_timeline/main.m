%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script plots time series of microbiota composition and selected meta
% data (white blood cells, infection, antibiotic administration, temperature)
% of MSK patients
% last updated: April. 23, 2020

addpath('../utils');

%% choose a patient id to plot
patientID2plot = '1511';

%% path to data
data_path = '../deidentified_data_tables/';

%% load samples table
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'));
tblsamples.PatientID = categorical(tblsamples.PatientID);

%% load counts table
tblcounts = readtable(strcat(data_path, 'counts/tblASVcounts_human_filter.csv'));

%% load taxonomy table
tbltaxonomy = readtable(strcat(data_path,'taxonomy/tblASVtaxonomy_silva_v4v5_filter.csv'));

%% load HCT information table
tblhctmeta = readtable(strcat(data_path, 'meta_data/tblhctmeta.csv'));
tblhctmeta.PatientID = categorical(tblhctmeta.PatientID);

%% load blood cell counts table
tblbc = readtable(strcat(data_path, 'meta_data/tblbc.csv'));
tblbc.PatientID = categorical(tblbc.PatientID);

%% load drug administration table
tbldrug = readtable(strcat(data_path, 'meta_data/tbldrug.csv'));
tbldrug.PatientID = categorical(tbldrug.PatientID);

%% load temperature table
tbltemp = readtable(strcat(data_path, 'meta_data/tbltemperature.csv'));
tbltemp.PatientID = categorical(tbltemp.PatientID);

%% load bacterial infection table
tblinfection = readtable(strcat(data_path, 'meta_data/tblInfectionsCidPapers.csv'));
tblinfection.PatientID = categorical(tblinfection.PatientID);

%% define a time window during which period data will be plotted
% the range should be expressed relative to HCT date
% if HCT date is missing for the given patient, the middle day of all
% samples is used as the reference point
window2plot = [-5,21];

%% define what blood cell data to plot (use column name in the table)
bc2plot = 'Neutrophils_unit_K_per_uL'; % neutrophils

%% plot time series of microbiome composition, white blood cell count, anti-infective drug administration, temperature
plot_community(patientID2plot, window2plot, bc2plot, tblsamples, tblcounts, tbltaxonomy, tblhctmeta, tblbc, tbldrug, tbltemp, tblinfection);
