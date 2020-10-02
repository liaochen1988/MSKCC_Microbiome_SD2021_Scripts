%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script plots time series of microbiota composition and selected meta
% data (white blood cells, infection, antibiotic administration, temperature)
% of MSK patients
% last updated: Oct. 1, 2020

addpath('../utils');
data_path = '../deidentified_data_tables/'; % path to data

%% choose a patient id to plot
PatientID2Plot = '1511';

%% define a time window during which period data will be plotted
% the range should be expressed relative to HCT date
% this exemplified code assumes this range is relative to the first HCT day
% if multiple HCT exists
RelativeTimePeriod2Plot = [-5,21];

%% define what blood cell data to plot (use column name in the table)
BloodCellType2Plot = 'Neutrophils'; % neutrophils

%% load HCT information
% this exemplified code only works for patients that have documented bone
% marrow transplant date (time point). Quit otherwise.
tblhctmeta = readtable(strcat(data_path, 'meta_data/tblhctmeta.csv'));
tblhctmeta.PatientID = categorical(tblhctmeta.PatientID);
tblhctmeta = tblhctmeta(tblhctmeta.PatientID==PatientID2Plot, :);
if (height(tblhctmeta)==0)
    warning("No HCT information for patient %s. Quit.", PatientID2Plot);
    return;
end
if (isnan(tblhctmeta{1,'TimepointOfTransplant'}))
    warning("Transplant day of patient %s is not available. Quit.", PatientID2Plot);
    return;
end

%% load samples table
% at least one stool sample exists for any patient with HCT meta data
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'));
tblsamples.PatientID = categorical(tblsamples.PatientID);
tblsamples = tblsamples((tblsamples.PatientID==PatientID2Plot) & (tblsamples.DayRelativeToNearestHCT >= RelativeTimePeriod2Plot(1)) & (tblsamples.DayRelativeToNearestHCT <= RelativeTimePeriod2Plot(2)), :);
tblsamples = sortrows(tblsamples, 'Timepoint'); % sort rows by time point of samples

%% find corresponding timepoint of TimePeriodRelativeToHCT
curr_relative_day = tblsamples{1,'DayRelativeToNearestHCT'};
for i=2:height(tblsamples)
    if (tblsamples{i,'DayRelativeToNearestHCT'}<curr_relative_day)
        tblsamples = tblsamples(1:i-1,:);
        break;
    else
        curr_relative_day = tblsamples{i,'DayRelativeToNearestHCT'};
    end
end
AbsoluteTimePeriod = [RelativeTimePeriod2Plot(1)-tblsamples{1,'DayRelativeToNearestHCT'}+tblsamples{1,'Timepoint'},...
                      RelativeTimePeriod2Plot(2)-tblsamples{1,'DayRelativeToNearestHCT'}+tblsamples{1,'Timepoint'}];

%% load counts table
tblcounts = readtable(strcat(data_path, 'counts/tblASVcounts_human_filter.csv'));
tblcounts = tblcounts(contains(tblcounts.SampleID, tblsamples.SampleID), :);

%% unstack counts table and normalize ASV counts to relative abundance
tblcounts = unstack(tblcounts, 'Count', 'ASV');
counts_matrix = tblcounts{:, 2:end}; % the first column is "SampleID"
counts_matrix(isnan(counts_matrix)) = 0; % missing count value is filled with 0
counts_matrix = counts_matrix ./ sum(counts_matrix, 2); % convert to relative abundance
tblcounts{:, 2:end} = counts_matrix;
tblcounts = innerjoin(tblsamples(:, {'SampleID', 'Timepoint', 'DayRelativeToNearestHCT'}), tblcounts);
tblcounts = sortrows(tblcounts, 'Timepoint'); % sort rows by time point of samples

%% load taxonomy table
tbltaxonomy = readtable(strcat(data_path,'taxonomy/tblASVtaxonomy_silva_v4v5_filter.csv'));
tbltaxonomy = tbltaxonomy(ismember(tbltaxonomy.ASV,tblcounts.Properties.VariableNames(4:end)), :);

%% load blood cell counts table
opts = detectImportOptions(strcat(data_path, 'meta_data/tblbc.csv'));
tblbc = readtable(strcat(data_path, 'meta_data/tblbc.csv'), opts);
tblbc.PatientID = categorical(tblbc.PatientID);
tblbc.BloodCellType = categorical(tblbc.BloodCellType);
tblbc = tblbc((tblbc.PatientID == PatientID2Plot) & (tblbc.BloodCellType==BloodCellType2Plot) & (tblbc.Timepoint >= AbsoluteTimePeriod(1)) & (tblbc.Timepoint <= AbsoluteTimePeriod(2)), :);
tblbc = sortrows(tblbc, 'Timepoint');
if (isempty(tblbc))
    warning("No %s data for patient %s.", BloodCellType2Plot, PatientID2Plot);
else
    tblbc = tblbc(:,{'DayRelativeToNearestHCT','Value'});
end

%% load drug administration table
opts = detectImportOptions(strcat(data_path, 'meta_data/tbldrug.csv'));
tbldrug = readtable(strcat(data_path, 'meta_data/tbldrug.csv'), opts);
tbldrug.PatientID = categorical(tbldrug.PatientID);
tbldrug = tbldrug(strcmp(tbldrug.AntiInfective,'True'),:); % select for only anti-infectives
tbldrug = tbldrug((tbldrug.PatientID == PatientID2Plot) & (tbldrug.StartTimepoint <= AbsoluteTimePeriod(2)) & (tbldrug.StopTimepoint >= AbsoluteTimePeriod(1)), :);
if (isempty(tbldrug))
    warning("No anti-infective drug records for patient %s.", PatientID2Plot);
end

%% load temperature table
tbltemp = readtable(strcat(data_path, 'meta_data/tbltemperature.csv'));
tbltemp.PatientID = categorical(tbltemp.PatientID);
tbltemp = tbltemp((tbltemp.PatientID == PatientID2Plot) & (tbltemp.Timepoint >= AbsoluteTimePeriod(1)) & (tbltemp.Timepoint <= AbsoluteTimePeriod(2)), :);
if (isempty(tbltemp))
    warning("No temperature data for patient %s.", PatientID2Plot);
end

%% load bacterial infection table
tblinfection = readtable(strcat(data_path, 'meta_data/tblInfectionsCidPapers.csv'));
tblinfection.PatientID = categorical(tblinfection.PatientID);
tblinfection = tblinfection((tblinfection.PatientID == PatientID2Plot) & (tblinfection.Timepoint >= AbsoluteTimePeriod(1)) & (tblinfection.Timepoint <= AbsoluteTimePeriod(2)), :);
if (isempty(tblinfection))
    warning("No bacterial infection data for patient %s.", PatientID2Plot);
end

%% plot time series of microbiome composition, blood cell count, anti-infective drug administration, temperature, and indicate infections if possible
plot_community(PatientID2Plot, RelativeTimePeriod2Plot, BloodCellType2Plot, tblcounts, tbltaxonomy, tblbc, tbldrug, tbltemp, tblinfection);