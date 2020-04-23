%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script plots tSNE scores of microbiota composition data of all MSK patients in a reduced 2-dimensional space
% and overlays compositional change of microbiome of a selected patient
% last updated: April. 23, 2020

addpath('../utils');

%% path to data
data_path = '../deidentified_data_tables/';

%% choose a patient id to show imeline trajectory on top of tSNE plot
patientID2plot = '1511';

%% find microbiome samples of the patient and sort by sample collection date
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'));
tblsamples.PatientID = categorical(tblsamples.PatientID);
tblsamples = tblsamples(tblsamples.PatientID==patientID2plot, :);

%% This should not happen as repeated samples for the same patient on the same day were removed from our dataset
% In case this happens, pick the first sample
tblsamples = sortrows(tblsamples, 'Day');
[~,idx2keep] = unique(tblsamples,'first');
tblsamples = tblsamples(idx2keep,:);

%% load and unstack the counts table
tblcounts = readtable(strcat(data_path, 'counts/tblASVcounts_human_filter.csv'));
tblcounts = unstack(tblcounts, 'Count', 'ASV');

%% load the taxonomy table
tbltaxonomy = readtable(strcat(data_path, 'taxonomy/tblASVtaxonomy_silva_v4v5_filter.csv'));

%% keep only ASVs in counts table that also appears in the taxonomy table
ASV_start_idx = find(contains(tblcounts.Properties.VariableNames, 'ASV'), 1); % index of first ASV in column names
[sharedASVs, idx, ~] = intersect(tbltaxonomy.ASV, tblcounts.Properties.VariableNames(ASV_start_idx:end), 'stable');
tblcounts = tblcounts(:, [tblcounts.Properties.VariableNames(1:ASV_start_idx-1),sharedASVs{:}]);
tbltaxonomy = tbltaxonomy(idx,:);

%% normalize ASV counts to relative abundance
matcounts = tblcounts{:, ASV_start_idx:end};
matcounts(isnan(matcounts)) = 0;
matcounts = matcounts ./ sum(matcounts, 2);
tblcounts{:, ASV_start_idx:end} = matcounts;

%% IMPORTANT! Run tSNE on the entire dataset takes ~ 30 minutes on my computer (MacBook Pro, 2018, 6 cpu parallelization)
% To reduce runninng time, you can set "test_mode" to true, which randomly
% choose ~1000 samples (including all samples of the chosen patient)
test_mode = false;
if (test_mode)
    samples2test = tblcounts.SampleID(randperm(height(tblcounts), 1000));
    for i=1:height(tblsamples)
        if(~strcmp(samples2test, tblsamples.SampleID(i)))
            samples2test(end+1) = tblsamples.SampleID(i);
        end
    end
    [~, samples2test_idx] = ismember(samples2test,tblcounts.SampleID);
    tblcounts = tblcounts(samples2test_idx, :);
end

%% get the dominant taxa of each sample and its rgb color
[~, dominant_ASV_idx] = max(tblcounts{:, ASV_start_idx:end}, [], 2);
ASV_color = hex2rgb(tbltaxonomy(dominant_ASV_idx, :).HexColor);

%% tSNE using functions downloaded from https://lvdmaaten.github.io/tsne/code/bh_tsne.tar.gz
% note: tSNE coefficients vary from run to run
original_dir = pwd;
software_path = '../utils/';
fastTSNE_path = strcat(software_path, 'bh_tsne');
cd(fastTSNE_path);
system('g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2');
fprintf('running tSNE...\n');
scoreLin = fast_tsne(tblcounts{:, ASV_start_idx:end}, 2, 20, 30, 0.5);
cd(original_dir);
fprintf('tSNE done.\n');

%% plot tSNE scores in 2 a reduced 2-dimensional space
figure();
hold on;

scatter(scoreLin(:, 1), scoreLin(:, 2), 100, ASV_color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
axis square;
box on;
set(gca,'Color',[176,196,222]/255); % set background color
xlabel('tSNE 1');
ylabel('tSNE 2');

xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;

axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
set(gca,'XTick',xtick_lb:10:xtick_ub);
set(gca,'YTick',ytick_lb:10:ytick_ub);

%% plot microbiome trajectory of the chosen patient on top of tSNE

% define a time window during which period data will be plotted
% the range should be expressed relative to HCT date
% if HCT date is missing for the given patient, the middle day of all
% samples is used as the reference point
window2plot = [-5,21];

% load HCT information table
tblhctmeta = readtable(strcat(data_path, 'meta_data/tblhctmeta.csv'));
tblhctmeta.PatientID = categorical(tblhctmeta.PatientID);

% find reference date to which all sample dates are relative
hct_idx = find(tblhctmeta.PatientID == patientID2plot);
if (isempty(hct_idx))
    warning('No hematopoietic cell transplant date found for patient %s.', patientID2plot);
    reference_day = tblsamples.Day(ceil(end/2));
else
    reference_day = tblhctmeta(hct_idx, :).DayOfTransplant;
    % in case the patient has received multiple HCT days, give warning and use the most recent one
    if (length(reference_day)>1)
        warning('More than one hematopoietic cell transplant dates found for patient %s. Use the most recent one.', patientID2plot);
        reference_day = max(reference_day);
    end
end
tblsamples.Day = tblsamples.Day-reference_day;

% keep data within the time window for plot
idx2keep = find(tblsamples.Day >= window2plot(1) & tblsamples.Day <= window2plot(2));
if (isempty(idx2keep))
    warning('No samples for patient %s within [%d, %d] days of the reference date.\n', patientID2plot, window2plot(1), window2plot(2));
else
    tblsamples = tblsamples(idx2keep, :);
end

% find the row index of the chosen patient's samples in the counts table
[~,idx2keep] = ismember(tblsamples.SampleID, tblcounts.SampleID);

% plot lines that connect consecutive samples and draw arrows at the middle of each line connection
x = scoreLin(idx2keep(1:end-1), 1);
y = scoreLin(idx2keep(1:end-1), 2);
u = (scoreLin(idx2keep(2:end), 1)-x)/2;
v = (scoreLin(idx2keep(2:end), 2)-y)/2;
for i=1:length(idx2keep)-1
    h = annotation('arrow');
    set(h,'parent', gca, ...
        'position', [x(i),y(i),u(i),v(i)], ...
        'HeadLength', 10, 'HeadWidth', 10, 'HeadStyle', 'cback1', ...
        'Color', 'w', 'LineWidth', 2, 'LineStyle', 'none');
end
plot(scoreLin(idx2keep, 1), scoreLin(idx2keep, 2), 'w-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
set(gca,'FontSize',14);