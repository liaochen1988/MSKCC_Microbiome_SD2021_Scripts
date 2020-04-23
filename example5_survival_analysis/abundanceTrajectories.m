function [] = abundanceTrajectories(tblAsvs, tblSample, tblSamplePatientDay, tblInfections,...
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
tblSampleTotalCounts = grpstats(tblSample, 'SampleID', 'sum', 'DataVars', 'Count');

%% make table with the ASVs of interest
tblSample = innerjoin(tblSample, tblFocusAsvs(:, {'ASV' 'Genus'}));

%%
tblSample = innerjoin(tblSample, tblSampleTotalCounts(:, {'SampleID' 'sum_Count'}));

%% compute fraction of each ASV of the GENUS in each sample
tblSample.fraction = tblSample.Count ./tblSample.sum_Count;

%% Add abundance of all ASVs of GENUS in a sample
tblSample = grpstats(tblSample, 'SampleID', 'sum', 'DataVars', 'fraction');

%% add the PatientID and the day of sample collection
tblSample = outerjoin(tblSamplePatientDay(:, {'SampleID' 'PatientID' 'Day'}), tblSample, 'Type', 'left',...
    'MergeKeys', true);
tblSample.sum_fraction(isnan(tblSample.sum_fraction)) = 0;

%% keep only samples between MINDAY and MAXDAY, exlude NaN and exclue FMT patients
tblSample(tblSample.Day<MINDAY | tblSample.Day>MAXDAY | isnan(tblSample.Day), :) = [];
if EXCLUDEFMT
    tblSample(contains(tblSample.SampleID, 'FMT'), :) = [];
end

%% find is the patient has an GENUS infection
tblInfections(~contains(tblInfections.InfectiousAgent, INFECTION), :) = [];
tblInfections(tblInfections.Day<MINDAY | tblInfections.Day>MAXDAY | isnan(tblInfections.Day), :) = [];
if EXCLUDEFMT
    tblInfections(contains(tblInfections.PatientID, 'FMT'), :) = []; 
end
% keep only patients present in the samples table
tblInfections(~ismember(tblInfections.PatientID, tblSample.PatientID), :) = [];
% keep only the first infection
tblInfections = grpstats(tblInfections, 'PatientID', 'min', 'DataVars', 'Day');
tblInfections = tblInfections(:, {'PatientID' 'min_Day'});
tblInfections.Properties.VariableNames{2} = 'infectionDay';

%% Plot all patients that develop invection
tblSample = outerjoin(tblSample, tblInfections, 'Type', 'left',...
    'MergeKeys', true);


%% make list of patients that says if the patient had infection
tblPatient = unique(tblSample(:, {'PatientID' 'infectionDay'}));


%% unstack
% plot samples of infected patients
m = fitglm(tblSample, 'sum_fraction ~ Day - 1',...
    'CategoricalVars', 'Day',...
    'Exclude', isnan(tblSample.infectionDay));
day = cellfun(@(x) str2num(strrep(x, 'Day_', '')), m.CoefficientNames);
plot(day, m.Coefficients.Estimate, '-', 'LineWidth', 2, 'Color', colorValue)
hold on;
% plot the uninfected patients
m = fitglm(tblSample, 'sum_fraction ~ Day - 1',...
    'CategoricalVars', 'Day',...
    'Exclude', ~isnan(tblSample.infectionDay));
day = cellfun(@(x) str2num(strrep(x, 'Day_', '')), m.CoefficientNames);
plot(day, m.Coefficients.Estimate, 'k', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
hold off
set(gca, 'YLim', [0 1])

fprintf('%s: %d patients infected, %d patients uninfected\n',...
    INFECTION,...
    length(unique(tblSample.PatientID(~isnan(tblSample.infectionDay)))),...
    length(unique(tblSample.PatientID(isnan(tblSample.infectionDay)))))

