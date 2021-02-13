% costumizable parameters
dataDir = '../deidentified_data_tables/';

% load sample table
tblSamples = readtable([dataDir 'samples/tblASVsamples.csv'],...
    'Format','%s%s%f%s%s%s%f');

% load infection table
tblInfections = readtable([dataDir 'meta_data/tblInfectionsCidPapers.csv'],...
    'Format','%s%f%s%f');

% load count table
tblCounts = readtable([dataDir 'counts/tblcounts_asv_melt.csv']);

% load taxonomy table
tblAsvs = readtable([dataDir 'taxonomy/tblASVtaxonomy_silva132_v4v5_filter.csv']);

%% plot the incidences of bacterial infections relative to day of transplant
MINDAY = -15;
MAXDAY = 35;

tblInfectionsUnStacked = tblInfections;
tblInfectionsUnStacked(tblInfectionsUnStacked.DayRelativeToNearestHCT<MINDAY | tblInfectionsUnStacked.DayRelativeToNearestHCT > MAXDAY, :) = [];
tblInfectionsUnStacked.PatientID = [];
tblInfectionsUnStacked.Timepoint = [];
tblInfectionsUnStacked.count = ones(height(tblInfectionsUnStacked), 1);
tblInfectionsUnStacked = unstack(tblInfectionsUnStacked,...
    'count', 'InfectiousAgent');
m = tblInfectionsUnStacked{:, 2:end};
m(isnan(m)) = 0;
tblInfectionsUnStacked{:, 2:end} = m;

%% plot the distribution of positive blood culters relative to HCT
% make color map
cmap = gray(size(m, 2));
organismDetected = tblInfectionsUnStacked.Properties.VariableNames(2:end);

figure(1)
set(gcf, 'Position', [44         407        1397         398]);
ba = bar(tblInfectionsUnStacked.DayRelativeToNearestHCT, tblInfectionsUnStacked{:, 2:end}, 'stacked');

green = 1;
red = 1;
for i = 1:length(ba)
    if contains(organismDetected{i}, 'enterococcus', 'IgnoreCase',true)
        ba(i).FaceColor = [0.2 0.4 0.2]*green ;
        green = green*0.9;
    elseif contains(organismDetected{i}, 'escherichia', 'IgnoreCase',true)
        ba(i).FaceColor = [0.9 0 0]*red;
        red = red*0.9;
    else
        ba(i).FaceColor = cmap(i, :);
    end
end
xlabel('Day relative to transplant')
ylabel('Number of positive blood cultures')
h = legend(organismDetected, 'Location', 'eastoutside');
set(h,'Interpreter', 'none');

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, 'fig5.eps', '-depsc');

%% plot the tajectories of enterococcus abundances in patients that develop
% infections
figure(2)
set(gcf, 'Position', [274    97   954   708]);
subplot(2, 2, 1)
GENUS = 'Enterococcus';
INFECTION = 'Enterococcus';
EXCLUDEFMT = true; % exclude any patients enroled in the FMT trial
abundanceTrajectories(tblAsvs, tblCounts, tblSamples, tblInfections,...
    GENUS, INFECTION, MINDAY, MAXDAY, EXCLUDEFMT, [0.2 0.4 0.2])
xlabel('Day')
ylabel('Average Enterococcus in stool')
set(gca, 'XLim', [-15 35])
set(gca, 'YLim', [0 1])
legend({'Infected patients' 'Uninfected patients'},...
    'Location', 'northwest')
grid on

subplot(2, 2, 2)
GENUS = 'Escherichia';
INFECTION = 'Escherichia';
EXCLUDEFMT = true; % exclude any patients enroled in the FMT trial
abundanceTrajectories(tblAsvs, tblCounts, tblSamples, tblInfections,...
    GENUS, INFECTION, MINDAY, MAXDAY, EXCLUDEFMT, [0.9 0 0])
xlabel('Day')
ylabel('Average Escherichia in stool')
set(gca, 'XLim', [-15 35])
set(gca, 'YLim', [0 1])
legend({'Infected patients' 'Uninfected patients'},...
    'Location', 'northwest')
grid on

%%
GENUS = 'Enterococcus';
INFECTION = 'Enterococcus';
EXCLUDEFMT = true; % exclude any patients enroled in the FMT trial
dominationEnterococcus = [0.001 0.01 0.1 0.3];

h = waitbar(0, 'Please wait...');
for i = 1:length(dominationEnterococcus)
    waitbar(i/length(dominationEnterococcus),h)
    [bEnterococcus(i),logl,H,stats] = survivalAnalysis(tblAsvs, tblCounts, tblSamples, tblInfections,...
        GENUS, INFECTION, MINDAY, MAXDAY, dominationEnterococcus(i), EXCLUDEFMT);
    pEnterococcus(i) = stats.p;
    errEnterococcus(i) = stats.se*1.96;
end
close(h);

figure(2),
subplot(2, 2, 3)
h = bar(exp(bEnterococcus(end:-1:1)));
h(1).FaceColor = [0.2 0.4 0.2];
set(gca, 'XTickLabel', {'30%' '10%' '1%' '0.1%'});
xlabel('Stool Enterococcus threshold')
ylabel('HR of Enterococcus positive blood culture')
hold on
plot(get(gca, 'XLim'), [1 1], 'k--')
hold off
set(gca, 'YLim', [0 10])
set(gca, 'YTick', 0:10)

%%
GENUS = 'Escherichia';
INFECTION = 'Escherichia';
EXCLUDEFMT = true; % exclude any patients enroled in the FMT trial
dominationEscherichia = [0.001 0.01 0.1 0.3];

h = waitbar(0, 'Please wait...');
for i = 1:length(dominationEscherichia)
    waitbar(i/length(dominationEscherichia),h)
    [bEscherichia(i),logl,H,stats] = survivalAnalysis(tblAsvs, tblCounts, tblSamples, tblInfections,...
        GENUS, INFECTION, MINDAY, MAXDAY, dominationEscherichia(i), EXCLUDEFMT);
    pEscherichia(i) = stats.p;
    errEscherichia(i) = stats.se*1.96;
end
close(h);


figure(2),
subplot(2, 2, 4)
h = bar(exp(bEscherichia(end:-1:1)));
h(1).FaceColor = [0.9 0 0];
set(gca, 'XTickLabel', {'30%' '10%' '1%' '0.1%'});
xlabel('Stool Escherichia threshold')
ylabel('HR of Escherichia positive blood culture')
hold on
plot(get(gca, 'XLim'), [1 1], 'k--')
hold off
set(gca, 'YLim', [0 10])
set(gca, 'YTick', 0:10)

%
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, 'fig6.eps', '-depsc');