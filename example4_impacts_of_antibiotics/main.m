%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script runs linear regression to estimate effects of antibiotics on microbial ASVs
% last updated: April. 23, 2020

addpath('../utils');

%% data to path
data_path = '../deidentified_data_tables/'; % path to data

%% set up simulation parameters

% do we use pseudo fractions? No
use_pseudo_fraction = 0;

% do we inlude antibiotic effect administered at the last day of the period
% between two consecutive samples? Yes
%
% Explanantion from Jonas Schluter:
% Consider an antibiotic that rapidly kills it's targets. Say this is administered on Wednesday.
% You got poop on Tuesday,  Wednesday and Thursday. 
% If you consider this antibiotic to be "on" only during the Wednesday to Thursday interval then the following can happen:
% the drug was given a few hours before stool was collected. All the killing happened. 
% You observe a flat line between wednesday and Thursday, inferring "no effect" of the drug. 
% Therefore, we ought to "backpropagate" any administration, i.e. consider any drug "on" when it was given at either end of the interval. 
% Thus in the example, the drug would be on during the tuesday-wednesday and Wednesday-thursday interval.
% This might sometimes underestimate the true slope (by stretching it out too much) but does not risk missing the very strongest effects entirely.
include_lastday_drug = 1;

% include one or more of the following:
% 'oral', 'gastric or jejunum tube', 'injection', 'intravenous', 'topical', 'unknown', 'urethral/vaginal/irrigation'
drug_admin_route = {'intravenous'};

max_interval = 3; % maximum sample intervals considered in regression

num_ASVs = 20; % maximum number of ASVs considered in regression

%% create X and Y matrices for regression (Y ~ X)
[tblX, tblY, patientID_of_rows] = create_regression_matrices(data_path, use_pseudo_fraction, include_lastday_drug, drug_admin_route);

%% cut X and Y matrices by limiting sample intervals
logic2keep = tblX.growth <= max_interval;
tblX = tblX(logic2keep,:);
tblY = tblY(logic2keep,:);
patientID_of_rows = patientID_of_rows(logic2keep,:); % patientID

%% cut table Y by keeping only top ASVs with most finite data
[~,sorted_idx] = sort(sum(isfinite(tblY{:,:})), 'descend');
tblY = tblY(:,sorted_idx(1:num_ASVs));
ASVs = tblY.Properties.VariableNames;

%% remove ASVs with all non-finite values (NaN, Inf/-Inf)
tblY = tblY(:,any(isfinite(tblY{:,:})));

%% remove X matrix columns with all zero values
tblX = tblX(:,~all(tblX{:,:}==0));

%% get predictors (variable names in the X matrix)
predictors = tblX.Properties.VariableNames;

%% cross validation using lasso
fprintf('running cross-validated fitting ...\n');
tic;

figure('Name','Cross validation error');
hold on;

opt_coefs = ones(length(predictors),num_ASVs) * NaN; % optimal coefficients
observed_logdiff = [];  % observed log-difference of absolute abundance
predicted_logdiff = []; % predicted log-difference of absolute abundance

kfold = 3;  % k-fold cross validation
nrep  = 10; % number of repetitions
lambda = 10.^[-2:0.1:2]; % coefficients of panelty
nrows2plot = 4;  % number of rows in plot
ncols2plot = 5;  % number of columns in plot

for i=1:num_ASVs

    % remove nonfinite number in Y matrix
    idx2keepY = isfinite(tblY{:,i});
    tblX_ASVi = tblX(idx2keepY,:);
    tblY_ASVi = tblY(idx2keepY,i);
    patientID_of_rows_ASVi = patientID_of_rows(idx2keepY,1);
    unique_patientID_ASVi = unique([patientID_of_rows_ASVi{:}]);
    
    % remove all-0 columns in X matrix
    idx2keepX = ~all(tblX_ASVi{:,:}==0);
    tblX_ASVi = tblX_ASVi(:,idx2keepX);

    % partition patient patientID to training and test datasets    
    rng default;
    c = cvpartition(length(unique_patientID_ASVi),'KFold',kfold);
    crossval_mse = zeros(nrep,length(lambda));
    for j=1:nrep
        cnew = c.repartition;
        for r=1:kfold
            idxTrain = ismember([patientID_of_rows_ASVi{:}], unique_patientID_ASVi(cnew.training(r)));
            idxTest = ~idxTrain;
            XTrain = tblX_ASVi{idxTrain,:};
            yTrain = tblY_ASVi{idxTrain,1};
            XTest = tblX_ASVi{idxTest,:};
            yTest = tblY_ASVi{idxTest,1};
            % alpha~=0 for ridge regression
            [opt_beta,~] = lasso(XTrain,yTrain,'Lambda',lambda,'Alpha',1e-4,'CV','resubstitution','Standardize',false,'RelTol',1e-6);     
            crossval_mse(j,:) = crossval_mse(j,:)+sum((XTest*opt_beta-repmat(yTest,1,length(lambda))).^2,1)/length(yTest);
        end
        crossval_mse(j,:) = crossval_mse(j,:)/kfold;
    end

    % find lambda that gives minimal cross validation error
    mean_crossval_mse = mean(crossval_mse,1);
    [~,idx_opt_lambda] = min(mean_crossval_mse);
    opt_lambda = lambda(idx_opt_lambda);
    
    % do lasso again using optimal lambda and all data
    [opt_beta,~] = lasso(tblX_ASVi{:,:},tblY_ASVi{:,1},'Lambda',opt_lambda,'RelTol',1e-6,'Alpha',1e-4,'Standardize',false);
    opt_coefs(idx2keepX,i) = opt_beta;
    
    % observed vs predicted log difference values
    observed_logdiff = [observed_logdiff;tblY_ASVi{:,1}];
    predicted_logdiff = [predicted_logdiff;tblX_ASVi{:,:}*opt_beta];
    
    % plot cross-validated fits
    subplot(nrows2plot,ncols2plot,i);
    plot(lambda, mean_crossval_mse, 'k.-');
    set(gca,'XScale','log');
    xlim([1e-2,1e2]);
    axis square;
    box on;
    set(gca,'XTick', 10.^[-2:2:2]);
    set(gca,'XTicklabel', {'10^{-2}'; '10^{0}'; '10^{2}'});
    title(strrep(ASVs{i}, '_', ' ')); 
    if (mod(i,ncols2plot)==1)
        ylabel('CV error');
    end
    if (ceil(i/ncols2plot)==nrows2plot)
        xlabel('lambda');
    end
end
toc;
fprintf('cross validation done.\n');

%% plot model vs. data
figure('Name','Model vs data');
hold on;
plot(predicted_logdiff, observed_logdiff, 'k.','MarkerSize',10);
axis square;
axis([-5,5,-20,20]);
set(gca,'XTick',[-5:5:5]);
set(gca,'YTick',[-20:10:20]);
xlabel('Predicted log-difference','FontSize',16);
ylabel('Observed log-difference','FontSize',16);
ax = gca;
ax.FontSize = 12;
box on;

% calculate Pearson correlation coefficients
pearson_corr = corr(predicted_logdiff, observed_logdiff);
title(strcat('Pearson correlation (', sprintf('%2.2f',pearson_corr), ')'));

%% turn opt_coefs into a table
tbl_opt_coefs = array2table(opt_coefs','VariableNames',predictors,'RowNames',ASVs);

%% visualize data using heatmap

% rename ASV to format 'taxonomy (ASV)' (taxonomy is the lowest taxa that is not unclassified)
tbltaxonomy = readtable(strcat(data_path, 'taxonomy/tblASVtaxonomy_silva_v4v5_filter.csv'));
ASVs_w_taxa = cell(length(ASVs),1);
for i=1:length(ASVs)
    taxonomy_i = tbltaxonomy(strcmp(tbltaxonomy.ASV, ASVs{i}), {'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'});
    lowest_available_taxonomy_classification_idx = find(~strcmp(taxonomy_i{:,:}, '<not present>'), 1, 'last');
    if(~isempty(lowest_available_taxonomy_classification_idx))
        ASVs_w_taxa{i} = sprintf('%s (%s)', ASVs{i}, taxonomy_i{1, lowest_available_taxonomy_classification_idx}{1});
    else
        ASVs_w_taxa{i} = sprintf('%s (Unknown)', ASVs{i});
    end
end

% make first letter upper case
for i=1:length(ASVs)
    ASVs_w_taxa{i} = strcat(upper(extractBefore(ASVs_w_taxa{i},2)), extractAfter(ASVs_w_taxa{i},1)); % first letter upper case
end
for i=1:length(predictors)
    predictors{i} = strcat(upper(extractBefore(predictors{i},2)), extractAfter(predictors{i},1)); % first letter upper case
end

% replace '_' with white space
ASVs_w_taxa = strrep(ASVs_w_taxa, '_', ' ');
predictors = strrep(predictors, '_', ' ');

% use redblue map
my_map = brewermap(1000, '*RdBu');

% replace NaN with zero
tbl_opt_coefs_replace_NaN_w_zero = tbl_opt_coefs;
if (sum(sum(isnan(tbl_opt_coefs_replace_NaN_w_zero{:,:})))>0)
    tbl_opt_coefs_replace_NaN_w_zero = fillmissing(tbl_opt_coefs_replace_NaN_w_zero,'constant',0);
end
    
% using heatmap
figure('Name','Heatmap');
h = heatmap(tbl_opt_coefs{:,2:end}, 'Colormap', my_map);
h.CellLabelFormat = '%.2f';
ax = gca;
ax.XData = predictors(2:end);
ax.YData = ASVs_w_taxa;
ax.FontSize=10;
set(h.NodeChildren(3), 'XTickLabelRotation', 90);

% edit colorbar
axs = struct(ax);
cbh = axs.Colorbar;
ylabel(cbh, 'Susceptibilities (1/day)');
caxis([-0.3, 0.3]);