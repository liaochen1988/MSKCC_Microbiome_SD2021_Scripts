function [b,logl,H,stats] = survivalAnalysis(tblAsvs, tblSample, tblSamplePatientDay, tblInfections,...
    GENUS, INFECTION, MINDAY, MAXDAY, DOMINATION, EXCLUDEFMT)

%% Joao Xavier, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script makes a survival analysis of the risk of VRE bacteremia
% due to intestinal domination by6 Enterococcus bacteria.
% The analysis follows the appraoch in Taur et al 2012 Clin. Infec. Dis.
% where domination is defined as 30% or more of intestinal composition
% and is used as a time dependent covariate in a Cox proportional hazards
% model. Analysis uses only the patients with at least one microbiota
% sample between day MINDAY and day MAXDAY relative to the transplant.
%
% last updated: Mar. 14, 2020 by Joao Xavier


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
tblSample.fraction = tblSample.Count ./ tblSample.sum_Count;

%% find dominatations (where 1 ASV is > DOMINATION)
tblSample.dominated = tblSample.fraction >= DOMINATION;
fprintf('found %d dominated samples out of %d (%0.1f%% samples)\n',...
    sum(tblSample.dominated), height(tblSample),...
    sum(tblSample.dominated)/height(tblSample) * 100);

%% get the PatientID and the day of sample collection
tblSample = outerjoin(tblSamplePatientDay(:, {'SampleID' 'PatientID' 'Day'}), tblSample, 'Type', 'left',...
    'MergeKeys', true);
tblSample.dominated(isnan(tblSample.dominated)) = false;

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


%% merge infection table with domination table
tblSample = outerjoin(tblSample, tblInfections, 'Type', 'left',...
    'MergeKeys', true);

%% if the patient never gets infection set infection day to infinite
tblSample.infectionDay(isnan(tblSample.infectionDay)) = Inf;
% remove samples after infection
tblSample(tblSample.Day > tblSample.infectionDay, :) = [];


%% get the first samplpe, last sample, and domination dates for each patient
tblLastSample = grpstats(tblSample, 'PatientID', {'min' 'max'}, 'DataVars', 'Day');
tblLastSample.GroupCount = [];
tblLastSample.Properties.VariableNames{end-1} = 'firstSample';
tblLastSample.Properties.VariableNames{end} = 'lastSample';

tblDomination = tblSample(tblSample.dominated == 1, {'PatientID' 'Day'});
tblDomination.Properties.VariableNames{end} = 'donimationDay';
tblDomination = sortrows(tblDomination, 'donimationDay');
[~, idx] = unique(tblDomination.PatientID);
tblDomination = tblDomination(idx, :);

tblPatientData = outerjoin(tblLastSample, tblDomination, 'Type', 'left', 'MergeKeys', true);
tblPatientData = outerjoin(tblPatientData, tblInfections, 'Type', 'left', 'MergeKeys', true);
tblPatientData = sortrows(tblPatientData, 'PatientID');

%% build table for survival analysis with time-depdendent covariates
nRows = height(tblPatientData) + sum(tblPatientData.donimationDay > tblPatientData.firstSample);
tblSurvival =...
    table('Size',[nRows 5],...
    'VariableTypes', {'string' 'double' 'double' 'double' 'logical'},...
    'VariableNames', {'PatientID' 'startT' 'stopT' 'dominated' 'censored'});
hInternal = waitbar(0, 'Please wait...');

c = 1;
for i = 1:height(tblPatientData)
%for i = 1:2
    waitbar(i/height(tblPatientData),hInternal)
    first = tblPatientData.firstSample(i);
    last = tblPatientData.lastSample(i);
    domination = tblPatientData.donimationDay(i);
    infection = tblPatientData.infectionDay(i);
    tblSurvival.PatientID(c) = tblPatientData.PatientID(i);
    tblSurvival.startT(c) = first;
    tblSurvival.stopT(c) = nanmax([last nanmin([infection MAXDAY])]) + 23/24;
    % the first row is only dominated if it domination occurred at MINDAY
    tblSurvival.dominated(c) = (domination == first);
    tblSurvival.censored(c) = ~(infection == floor(tblSurvival.stopT(c)));
    c = c+1;
    % if patient has domination after firstSample add another row
    if domination > first
        % correct the stop date of pre-domination
        tblSurvival.stopT(c-1) = domination-1/24;
        % add post domination period
        tblSurvival.PatientID(c) = tblPatientData.PatientID(i);
        tblSurvival.startT(c) = domination;
        tblSurvival.stopT(c) = nanmax([last nanmin([infection MAXDAY])]) + 23/24;
        tblSurvival.dominated(c) = 1;
        tblSurvival.censored(c) = ~(infection == floor(tblSurvival.stopT(c)));
        c = c+1;
    end
end
close(hInternal)

%% Cox proportional hazzards with time dependent covariate
[b,logl,H,stats] = coxphfit(tblSurvival.dominated,...
    [tblSurvival.startT tblSurvival.stopT],...
    'Censoring',tblSurvival.censored, 'Baseline', 0);

fprintf('Domination by %s (defined >%0.1e) increases %s infection by %0.1fX [%0.1fX to %0.1fX 95%% confidence]\n',...
    GENUS, DOMINATION, INFECTION, exp(b), exp(b-stats.se*1.96), exp(b+stats.se*1.96))
fprintf('p-value = %0.2e\n', stats.p)

fprintf('%d patients dominated of %d patients\n',...
    length(unique(tblSample.PatientID(tblSample.dominated))),...
    length(unique(tblSample.PatientID)));

