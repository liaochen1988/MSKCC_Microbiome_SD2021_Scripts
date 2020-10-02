%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script plots frequency of drugs administered for MSK patients
% (administration routes separated for oral and intravenous)
% last updated: April. 23, 2020

addpath('../utils');

%% path to data
data_path = '../deidentified_data_tables/meta_data/';

%% load drug table
opts = detectImportOptions(strcat(data_path, 'tbldrug.csv'));
tbldrug = readtable(strcat(data_path, 'tbldrug.csv'),opts);

%% keep only anti-bacterial drug classes
class2include = {'aminoglycosides', 'carbapenems', 'cephalosporins' ...
    'glycopeptide antibiotics', 'glycylcyclines', 'leprostatics', ...
    'lincomycin derivatives', 'macrolide derivatives', 'miscellaneous antibiotics',...
    'oxazolidinone antibiotics', 'penicillins', 'quinolones',...
    'sulfonamides', 'tetracyclines'};
tbldrug = tbldrug(ismember(tbldrug.Category, class2include), :);

%% remove atovaquone
tbldrug = tbldrug(~strcmp(tbldrug.Factor, 'atovaquone'), :);

%% replace ' ' and '/' with '_' in drug name and category
tbldrug.Factor = strrep(tbldrug.Factor, ' ', '_');
tbldrug.Factor = strrep(tbldrug.Factor, '/', '_');
tbldrug.Category = strrep(tbldrug.Category, ' ', '_');
tbldrug.Category = strrep(tbldrug.Category, '-', '_');

%%  add columns All_Routes_Administration, Intravenous_Administration, Oral_Administration
tbldrug = tbldrug(:,{'PatientID';'Factor';'Category';'Route'});
tbldrug{:,'All_Routes_Administration'} = true;
tbldrug{:,'Intravenous_Administration'} = strcmp(tbldrug.Route,'intravenous');
tbldrug{:,'Oral_Administration'} = strcmp(tbldrug.Route,'oral');
tbldrug = unique(tbldrug);

%% plot frequence of antibiotic use
keywords = {'Oral_Administration';'Intravenous_Administration';'All_Routes_Administration'};
for k=1:length(keywords)
    
    % unstack table
    tbldrug2 = tbldrug(:,{'PatientID';'Factor';keywords{k}});
    tbldrug_unstacked  = unstack(tbldrug2, keywords{k}, 'Factor');
    tbldrug_unstacked_matrix = tbldrug_unstacked{:,2:end};
    tbldrug_unstacked_matrix(isnan(tbldrug_unstacked_matrix)) = false;
    fraction_of_patients = sum(tbldrug_unstacked_matrix,1)/height(tbldrug_unstacked);
    
    % calculate total counts and relative frequency of each drug
    tblfreq = unique(tbldrug(:,{'Factor';'Category'}));
    [~,loc] = ismember(tbldrug_unstacked.Properties.VariableNames(2:end),tblfreq.Factor);
    tblfreq{loc,'Frequency'} = fraction_of_patients';
    
    % sort rows by category
    tblfreq = sortrows(tblfreq, 'Category');
    
    % barplot
    figure()
    hold on;
    
    unique_cat = unique(tblfreq.Category, 'stable');
    my_colors = distinguishable_colors(length(unique_cat));
    hbar = bar([1:height(tblfreq)], tblfreq.Frequency, 'BarWidth', 1);
    hbar.FaceColor = 'flat';
    
    for i=1:length(unique_cat)
        
        % assign the same color for all bars of the same category
        factor_idx_class_i = find(strcmp(tblfreq.Category, unique_cat{i}));
        hbar.CData(factor_idx_class_i,:) = repmat(my_colors(i,:), length(factor_idx_class_i), 1);
        
        % patch
        xdata = [factor_idx_class_i(1)-0.5,factor_idx_class_i(end)+0.5, factor_idx_class_i(end)+0.5, factor_idx_class_i(1)-0.5];
        ydata = [0,0,1,1];
        hpi = patch(xdata,ydata,my_colors(i,:),'FaceAlpha',0.2);
        uistack(hpi, 'bottom');
        
        % add text
        hti=text(mean(factor_idx_class_i), 0.6, strrep(unique_cat{i}, '_', ' '), 'FontSize', 14, 'HorizontalAlignment', 'center', 'Rotation', 90, 'Interpreter', 'None');
        uistack(hti, 'top');
    end
    set(gca,'xtick',[1:height(tblfreq)]);
    set(gca,'xticklabel', strrep(tblfreq.Factor, '_', ' '));
    xtickangle(90);
    set(gca,'ytick',[0:0.2:1.0]);
    box on;
    set(gca,'ticklength',[0,0]);
    grid on;
    set(gca,'GridLineStyle','-');
    set(gca,'GridAlpha',0.2);
    axis([0.5,height(tblfreq)+0.5, 0, 1]);
    ylabel('Fraction of patients receiving antibiotic');
    title(keywords{k}, 'interpreter', 'None');
end