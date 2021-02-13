%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script plots tSNE scores of microbiota composition data of all MSK patient samples in a reduced 2-dimensional space
% and overlays compositional change of microbiome of a selected patient
% last updated: Oct. 1, 2020

addpath('../utils');

%% path to data
data_path = '../deidentified_data_tables/';

%% load and unstack the counts table
tblcounts = readtable(strcat(data_path, 'counts/tblcounts_asv_melt.csv'));
tblcounts = unstack(tblcounts, 'Count', 'ASV');

%% normalize ASV counts to relative abundance
matcounts = tblcounts{:, 2:end}; % the first column is 'SampleID'
matcounts(isnan(matcounts)) = 0; % replace missing taxa count with 0
matcounts = matcounts ./ sum(matcounts, 2); % convert to relative abundance
tblcounts{:, 2:end} = matcounts;

%% load the taxonomy table
tbltaxonomy = readtable(strcat(data_path, 'taxonomy/tblASVtaxonomy_silva132_v4v5_filter.csv'));

%% get the dominant taxa of each sample and its rgb color
% each sample in the tSNE plot will be colored by its dominant taxa
[~, dominant_ASV_idx] = max(tblcounts{:, 2:end}, [], 2);
ASV_color = hex2rgb(tbltaxonomy(dominant_ASV_idx, :).HexColor);

%% tSNE using functions downloaded from https://lvdmaaten.github.io/tsne/code/bh_tsne.tar.gz
% IMPORTANT! Run tSNE on the entire dataset takes ~ 30 minutes on my PC (MacBook Pro, 2018, 6 cpu parallelization)
% note: tSNE coefficients vary from run to run
original_dir = pwd;
software_path = '../utils/';
fastTSNE_path = strcat(software_path, 'bh_tsne');
cd(fastTSNE_path);
system('g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2');
fprintf('running tSNE...\n');
scoreLin = fast_tsne(tblcounts{:, 2:end}, 2, 20, 30, 0.5);
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

% choose a patient id to show timeline trajectory on top of tSNE plot
patientID2plot = '1511';

% load sample table and find microbiome samples of the patient
opts = detectImportOptions(strcat(data_path, 'samples/tblASVsamples.csv'));
opts = setvartype(opts,{'PatientID'},'categorical');
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'), opts);
tblsamples = tblsamples(tblsamples.PatientID==patientID2plot, :);

% define a time window during which period data will be plotted
% the range should be expressed relative to HCT date
window2plot = [-5,21];

% load HCT information table
tblhctmeta = readtable(strcat(data_path, 'meta_data/tblhctmeta.csv'));
tblhctmeta.PatientID = categorical(tblhctmeta.PatientID);

% keep data within the time window for plot
tblsamples = tblsamples((tblsamples.DayRelativeToNearestHCT >= window2plot(1)) & (tblsamples.DayRelativeToNearestHCT <= window2plot(2)), :);

% find the row index of the chosen patient's samples in the count table
[~,count_table_index] = ismember(tblsamples.SampleID, tblcounts.SampleID);

% plot lines that connect consecutive samples and draw arrows at the middle of each line connection
x = scoreLin(count_table_index(1:end-1), 1);
y = scoreLin(count_table_index(1:end-1), 2);
u = (scoreLin(count_table_index(2:end), 1)-x)/2;
v = (scoreLin(count_table_index(2:end), 2)-y)/2;
for i=1:length(count_table_index)-1
    h = annotation('arrow');
    set(h,'parent', gca, ...
        'position', [x(i),y(i),u(i),v(i)], ...
        'HeadLength', 10, 'HeadWidth', 10, 'HeadStyle', 'cback1', ...
        'Color', 'w', 'LineWidth', 2, 'LineStyle', 'none');
end
plot(scoreLin(count_table_index, 1), scoreLin(count_table_index, 2), 'w-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
set(gca,'FontSize',14);