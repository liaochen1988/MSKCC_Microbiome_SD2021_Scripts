function [] = abundanceTrajectories(tblAsvs, tblCounts, tblSamples, tblInfections,...
    GENUS, INFECTION, MINDAY, MAXDAY, EXCLUDEFMT, colorValue)

%% Joao Xavier, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script makes a plot comparing the trajectories of abundances of
% GENUS in patients that develop infections by INFECTIO
%
% last updated: Mar. 22, 2020 by Joao Xavier


%% load ASV table that correspont to ASVs belonging to genus
tblFocusAsvs = tblAsvs;

%% keep only the GENUS of interest
tblFocusAsvs(~contains(tblFocusAsvs.Genus, GENUS), :) = [];
fprintf('found %d ASVs belonging to %s\n', height(tblFocusAsvs), GENUS);

%% sum the number of total counts per sample
tblSampleTotalCounts = grpstats(tblCounts, 'SampleID', 'sum', 'DataVars', 'Count');

%% make table with the ASVs of interest
tblCounts = innerjoin(tblCounts, tblFocusAsvs(:, {'ASV' 'Genus'}));

%%
tblCounts = innerjoin(tblCounts, tblSampleTotalCounts(:, {'SampleID' 'sum_Count'}));

%% compute fraction of each ASV of the GENUS in each sample
tblCounts.fraction = tblCounts.Count ./tblCounts.sum_Count;

%% Add abundance of all ASVs of GENUS in a sample
tblCounts = grpstats(tblCounts, 'SampleID', 'sum', 'DataVars', 'fraction');

%% add the PatientID and the day of sample collection
tblCounts = outerjoin(tblSamples(:, {'SampleID' 'PatientID' 'DayRelativeToNearestHCT'}), tblCounts, 'Type', 'left',...
    'MergeKeys', true);
tblCounts.sum_fraction(isnan(tblCounts.sum_fraction)) = 0;

%% keep only samples between MINDAY and MAXDAY, exlude NaN and exclue FMT patients
tblCounts(tblCounts.DayRelativeToNearestHCT<MINDAY | tblCounts.DayRelativeToNearestHCT>MAXDAY | isnan(tblCounts.DayRelativeToNearestHCT), :) = [];
if EXCLUDEFMT
    tblCounts(contains(tblCounts.SampleID, 'FMT'), :) = [];
end

%% find is the patient has an GENUS infection
tblInfections(~contains(tblInfections.InfectiousAgent, INFECTION), :) = [];
tblInfections(tblInfections.DayRelativeToNearestHCT<MINDAY | tblInfections.DayRelativeToNearestHCT>MAXDAY | isnan(tblInfections.DayRelativeToNearestHCT), :) = [];
if EXCLUDEFMT
    tblInfections(contains(tblInfections.PatientID, 'FMT'), :) = []; 
end
% keep only patients present in the samples table
tblInfections(~ismember(tblInfections.PatientID, tblCounts.PatientID), :) = [];
% keep only the first infection
tblInfections = grpstats(tblInfections, 'PatientID', 'min', 'DataVars', 'DayRelativeToNearestHCT');
tblInfections = tblInfections(:, {'PatientID' 'min_DayRelativeToNearestHCT'});
tblInfections.Properties.VariableNames{2} = 'infectionDay';

%% Plot all patients that develop invection
tblCounts = outerjoin(tblCounts, tblInfections, 'Type', 'left',...
    'MergeKeys', true);


%% make list of patients that says if the patient had infection
tblPatient = unique(tblCounts(:, {'PatientID' 'infectionDay'}));


%% unstack
% plot samples of infected patients
m = fitglm(tblCounts, 'sum_fraction ~ DayRelativeToNearestHCT - 1',...
    'CategoricalVars', 'DayRelativeToNearestHCT',...
    'Exclude', isnan(tblCounts.infectionDay));
day = cellfun(@(x) str2num(strrep(x, 'DayRelativeToNearestHCT_', '')), m.CoefficientNames);
plot(day, m.Coefficients.Estimate, '-', 'LineWidth', 2, 'Color', colorValue)
hold on;
% plot the uninfected patients
m = fitglm(tblCounts, 'sum_fraction ~ DayRelativeToNearestHCT - 1',...
    'CategoricalVars', 'DayRelativeToNearestHCT',...
    'Exclude', ~isnan(tblCounts.infectionDay));
day = cellfun(@(x) str2num(strrep(x, 'DayRelativeToNearestHCT_', '')), m.CoefficientNames);
plot(day, m.Coefficients.Estimate, 'k', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
hold off
set(gca, 'YLim', [0 1])

fprintf('%s: %d patients infected, %d patients uninfected\n',...
    INFECTION,...
    length(unique(tblCounts.PatientID(~isnan(tblCounts.infectionDay)))),...
    length(unique(tblCounts.PatientID(isnan(tblCounts.infectionDay)))))

