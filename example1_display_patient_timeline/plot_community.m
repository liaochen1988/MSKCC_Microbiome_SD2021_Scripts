function plot_community(patientID2plot, window2plot, bc2plot, tblsamples, tblcounts, tbltaxonomy, tblhctmeta, tblbc, tbldrug, tbltemp, tblinfection)

%% find samples of the patient with patientID "patientID2plot"
tblsamples = tblsamples(tblsamples.PatientID==patientID2plot, :);
if (isempty(tblsamples))
    warning("No samples data for patient %s.", patientID2plot);
end

%% This should not happen as repeated samples for the same patient on the same day were removed from our dataset
% In case this happens, pick the first sample
tblsamples = sortrows(tblsamples,'Day');
[~,idx2keep] = unique(tblsamples,'first');
tblsamples = tblsamples(idx2keep,:);

%% find ASV counts of these samples
tblcounts = tblcounts(contains(tblcounts.SampleID, tblsamples.SampleID), :);
if (isempty(tblcounts))
    warning("No ASV count data for patient %s.", patientID2plot);
end

%% unstack counts table and normalize ASV counts to relative abundance
tblcounts = unstack(tblcounts, 'Count', 'ASV');
ASV_start_idx = find(contains(tblcounts.Properties.VariableNames, 'ASV'), 1);
counts_matrix = tblcounts{:, ASV_start_idx:end};
counts_matrix(isnan(counts_matrix)) = 0;
counts_matrix = counts_matrix ./ sum(counts_matrix, 2);
tblcounts{:, ASV_start_idx:end} = counts_matrix;

%% join samples and counts by SampleID
tbljoined = innerjoin(tblsamples(:, {'SampleID', 'Day'}), tblcounts);
if (isempty(tbljoined))
    warning('No data after joining the sample and count table.');
end

%% sort rows by day of samples
tbljoined = sortrows(tbljoined, 'Day');

%% find taxonomy for ASVs in the joined table
ASV_start_idx = find(contains(tbljoined.Properties.VariableNames, 'ASV'), 1);
[sharedASVs, idx, ~] = intersect(tbltaxonomy.ASV, tbljoined.Properties.VariableNames(ASV_start_idx:end), 'stable');
tbljoined = tbljoined(:, {'SampleID','Day',sharedASVs{:}});
tbltaxonomy = tbltaxonomy(idx,:);

%% convert sample days to relative scale
% if HCT date is available, make samples days reltive to HCT date
% if HCT date is not available, make sample days relative to the median day of all samples
if (~isempty(tbljoined))
    hct_idx = find(tblhctmeta.PatientID == patientID2plot);
    if (isempty(hct_idx))
        warning('No hematopoietic cell transplant date found for patient %s.', patientID2plot);
        reference_day = tbljoined.Day(ceil(end/2));
    else
        reference_day = tblhctmeta(hct_idx, :).DayOfTransplant;
        % in case the patient has received multiple HCT days, give warning and use the most recent one
        if (length(reference_day)>1)
            warning('More than one hematopoietic cell transplant dates found for patient %s. Use the most recent one.', patientID2plot);
            reference_day = max(reference_day);
        end
    end
    tbljoined.Day = tbljoined.Day - reference_day;
else
    reference_day = 0;
end

%% keep data within the time window for plot
if (~isempty(tbljoined) && ~isempty(window2plot))
    idx2keep = find(tbljoined.Day >= window2plot(1) & tbljoined.Day <= window2plot(2));
    if (isempty(idx2keep))
        warning('No ASV count data for patient %s within [%d, %d] days of the reference date.', patientID2plot, window2plot(1), window2plot(2));
        tbljoined = [];
    else
        tbljoined = tbljoined(idx2keep, :);
    end
end

%% get blood cell counts for the patient with patientID "patientID2plot"
tblbc = tblbc(tblbc.PatientID == patientID2plot, :);
if (isempty(tblbc))
    warning("No white blood cell data for patient %s.", patientID2plot);
else
    tblbc.Day = tblbc.Day - reference_day;
    if (~isempty(window2plot))
        idx2keep = find(tblbc.Day >= window2plot(1) & tblbc.Day <= window2plot(2));
        if (isempty(idx2keep))
            warning('No blood cell data for patient %s within [%d, %d] days of the reference point.', patientID2plot, window2plot(1), window2plot(2));
            tblbc = [];
        else
            tblbc = tblbc(idx2keep, :);
        end
    end
    
    % impute missing values
    if (~isempty(tblbc))
        tblbc = fillmissing(tblbc,'spline','DataVariables',{bc2plot});
        tblbc{tblbc.(bc2plot)<0, bc2plot} = 0;
        
        % This should not happen as repeated measurement of blood cell counts have been removed
        % In case blood cell counts have duplicates; use their average
        [Day,~,ic] = unique(tblbc.Day,'stable');
        bc2plot_ave = accumarray(ic,tblbc.(bc2plot),[],@mean);
        tblbc = table(Day,bc2plot_ave,'VariableNames',{'Day',bc2plot});
    end
end

%% get anti-infective drug adminisitration records for the patient with patientID "patientID2plot"
tbldrug = tbldrug(strcmp(tbldrug.AntiInfective,'True'),:); % select for only anti-infectives
tbldrug = tbldrug(tbldrug.PatientID == patientID2plot, :);
if (isempty(tbldrug))
    warning("No anti-infective drug adminisitration data for patient %s.", patientID2plot);
else
    tbldrug.StartDay = tbldrug.StartDay - reference_day;
    tbldrug.StopDay = tbldrug.StopDay - reference_day;
    if (~isempty(window2plot))
        idx2keep = find(tbldrug.StopDay >= window2plot(1) & tbldrug.StartDay <= window2plot(2));
        if (isempty(idx2keep))
            warning('No anti-infective drug adminisitration data for patient %s within [%d, %d] days of the reference date.\n', patientID2plot, window2plot(1), window2plot(2));
            tbldrug = [];
        else
            tbldrug = tbldrug(idx2keep, :);
        end
    end
end

%% get temperature data for the patient with patientID "patientID2plot"
tbltemp = tbltemp(tbltemp.PatientID == patientID2plot, :);
if (isempty(tbltemp))
    warning("No temperature data for patient %s.", patientID2plot);
else
    tbltemp.Day = tbltemp.Day - reference_day;
    if (~isempty(window2plot))
        idx2keep = find(tbltemp.Day >= window2plot(1) & tbltemp.Day <= window2plot(2));
        if (isempty(idx2keep))
            warning('No temperature data for patient %s within [%d, %d] days of the reference date.', patientID2plot, window2plot(1), window2plot(2));
            tbltemp = [];
        else
            tbltemp = tbltemp(idx2keep, :);
        end
    end
end

%% get bacterial infection data for the patient with patientID "patientID2plot"
tblinfection = tblinfection(tblinfection.PatientID == patientID2plot, :);
if (isempty(tblinfection))
    warning("No infection data for patient %s.", patientID2plot);
else
    tblinfection.Day = tblinfection.Day - reference_day;
    if (~isempty(window2plot))
        idx2keep = find(tblinfection.Day >= window2plot(1) & tblinfection.Day <= window2plot(2));
        if (isempty(idx2keep))
            warning('No infection data for patient %s within [%d, %d] days of the reference date.', patientID2plot, window2plot(1), window2plot(2));
            tblinfection = [];
        else
            tblinfection = tblinfection(idx2keep, :);
        end
    end
end

%% plot time series of microbiome composition (indicate infections if available)
figure();
title(['patientID: ', patientID2plot]);
set(gca, 'FontName', 'Arial');

subplot(5,1,[4,5]);
hold on;
box on;
xtick_lb = round(window2plot(1)/5)*5;      % lower bound of xtick
xtick_ub = (round(window2plot(2)/5)+1)*5;  % upper bound of xtick
    
if (~isempty(tbljoined))
    
    % 1st and 2nd columns of tbljoined should be SampleID and day
    abundance_matrix = tbljoined{:, ASV_start_idx:end}; 

    % calculate the cumulative sum of taxa with same color_order
    % unique_color_order should be automatically sorted
    [unique_color_order,ia,ic] = unique(tbltaxonomy.ColorOrder);
    uni_color_hex = tbltaxonomy.HexColor(ia);
    
    color_grouped_abundance = zeros(size(abundance_matrix,1), length(unique_color_order));
    for k = 1:length(unique_color_order)
        currsum = sum(abundance_matrix(:,ic==k),2);
        color_grouped_abundance(:,k) = currsum;
    end
    
    % plot stacked bars
    if (height(tbljoined) == 1)
        h=bar([tbljoined.Day,nan], [color_grouped_abundance; nan(1, size(color_grouped_abundance,2))], 'stacked');
    else
        h=bar(tbljoined.Day, color_grouped_abundance, 'stacked');
    end
    
    % add text to indicate dominant ASV with more than 15% of relative abundance
    taxonomy_vars = tbltaxonomy.Properties.VariableNames;
    [dominant_ASV_row, dominant_ASV_col] = find(abundance_matrix > 0.15);
    for k = 1:length(dominant_ASV_row)
        curr_row = dominant_ASV_row(k); % sample
        curr_col = dominant_ASV_col(k); % asv
        dominant_tax = tbltaxonomy.(taxonomy_vars{strcmp(taxonomy_vars,'Genus')}){curr_col};
        if (strcmp(dominant_tax, '<not present>'))
            try_higher_level_classification = {'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'ASV'};
            for q=1:length(try_higher_level_classification)
                dominant_tax = tbltaxonomy.(taxonomy_vars{strcmp(taxonomy_vars,try_higher_level_classification{q})}){curr_col};     
                if(~strcmp(dominant_tax, '<not present>'))
                    break;
                end
            end
            dominant_tax = strcat('Unknown', {' '}, dominant_tax);
        end
        curr_color = tbltaxonomy.HexColor(curr_col); % taxonomy and joined_samples_counts share the same order of ASVs
        curr_box = find(strcmp(uni_color_hex,curr_color));
        cum_sum_cmap = cumsum(color_grouped_abundance(curr_row,:));
        heightbefore = [0 cum_sum_cmap(1:end-1)];
        heightafter = cum_sum_cmap;
        text(tbljoined.Day(curr_row), (heightafter(curr_box)+heightbefore(curr_box))/2, dominant_tax, ...
            'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 12);
    end
    
    % set barplot color
    ASV_cmap = hex2rgb(uni_color_hex);
    for cc = 1:size(color_grouped_abundance,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap(cc,:);
    end
    colormap(gca, ASV_cmap);
    
    ylim([0 1]);
    set(gca,'Ytick',[0:0.2:1]);
else
    set(gca,'YTick', []);
    set(gca,'YTicklabel', []);
end

% add indicators for infection
if (~isempty(tblinfection))
    for k=1:height(tblinfection)
        text(tblinfection.Day(k), 1, '*', 'HorizontalAlignment', 'center', 'FontSize', 50);
    end
end

if isempty(hct_idx)
    xlabel('Days relative to date of middle sample');
else
    xlabel('Days relative to HCT date');
end
xlim(window2plot);
set(gca,'Xtick',[xtick_lb:5:xtick_ub]);
ylabel('16S relative abundance');
set(gca,'FontSize',12);

%% plot time series of blood cell counts of type "bc2plot"
subplot(5,1,3);
hold on;
box on;

if (~isempty(tblbc))
    yyaxis left;
    plot(tblbc.Day, tblbc.(bc2plot), 'k.-', 'MarkerFaceColor', 'w', 'MarkerSize', 20, 'LineWidth', 1);
    if(contains(bc2plot, 'Neutrophils'))
        idx_neutropenia = find(tblbc.(bc2plot) <= 0.5); % threshold for neutropenia
        plot(tblbc.Day(idx_neutropenia), tblbc.(bc2plot)(idx_neutropenia), 'r.', 'MarkerSize', 20);
    end
    set(gca,'Xtick',xtick_lb:5:xtick_ub);
    
    ytick_lb = floor(min(tblbc.(bc2plot)));
    ytick_ub = ceil(max(tblbc.(bc2plot)));
    if (mod(ytick_lb,2)~=0)
        ytick_lb = ytick_lb-3;
    else
        ytick_lb = ytick_lb-2;
    end
    if (ytick_lb<0)
        ytick_lb=0;
    end
    if (mod(ytick_ub,2)~=0)
        ytick_ub = ytick_ub+3;
    else
        ytick_ub = ytick_ub+2;
    end
    ylim([ytick_lb ytick_ub]);
    set(gca,'Ytick',[ytick_lb,ytick_lb+(ytick_ub-ytick_lb)/2,ytick_ub]);
    
    % plot reference line for neutropenia
    yyaxis right;
    ylim([ytick_lb ytick_ub]);
    plot([window2plot(1), window2plot(2)], [0.5, 0.5], 'k--');
    set(gca,'Ytick',0.5);
    
    % set black axis
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    yyaxis left;
else
    set(gca,'YTick', []);
    set(gca,'YTicklabel', []);
    set(gca,'XTick', []);
end

xlabel('');
xlim(window2plot);
set(gca, 'XTicklabel', []);
if(contains(bc2plot,'RBCtotal'))
    ylabel({strrep(bc2plot,'_unit_M_per_uL',''), 'count (M/\muL)'});
else
    ylabel({strrep(bc2plot,'_unit_K_per_uL',''), 'count (K/\muL)'});
end
set(gca,'FontSize',12);

%% plot time series for anti-infective drug administration
subplot(5,1,2);
hold on;
box on;

boxwidth = 0.4;
boxheight = 0.4;
color2patch = [192,192,192]/255;

if (~isempty(tbldrug))
    unique_factors = unique(tbldrug.Factor);
    for i=1:length(unique_factors)
        factor_i_idx = find(strcmp(tbldrug.Factor, unique_factors(i)));
        if (isempty(factor_i_idx))
            continue
        else
            tbldrug_factor_i = tbldrug(factor_i_idx, :);
            administered_days =  []; % days that factor i was administered
            for j=1:height(tbldrug_factor_i)
                administered_days = [administered_days,tbldrug_factor_i.StartDay(j):1:tbldrug_factor_i.StopDay(j)];              
                % plot a single period from first day to last day administered for factor i
                patch([tbldrug_factor_i.StartDay(j)-boxwidth, tbldrug_factor_i.StopDay(j)+boxwidth, tbldrug_factor_i.StopDay(j)+boxwidth, tbldrug_factor_i.StartDay(j)-boxwidth], ...
                      [i-boxheight, i-boxheight, i+boxheight, i+boxheight], color2patch, 'EdgeColor', 'none');
            end
            % add text
            center_x = (min(administered_days) + max(administered_days))/2;
            text(center_x, i, unique_factors(i), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'None');
        end
    end

    ytick_lb = 0;
    ytick_ub = length(unique_factors);
    if (mod(ytick_lb,2)~=0)
        ytick_ub = ytick_ub+1;
    else
        ytick_ub = ytick_ub+2;
    end
    ylim([0 ytick_ub]); 
    
    set(gca,'Xtick',xtick_lb:5:xtick_ub);
else
    set(gca,'XTick', []);
end

xlabel('');
xlim(window2plot);
set(gca, 'XTicklabel', []);
ylabel('Medications');
set(gca,'Ytick',[]); 
set(gca,'Yticklabel',[]); 
set(gca,'FontSize',12);

%% plot time series for temperature
subplot(6,1,1);
hold on;
box on;

if (~isempty(tbltemp))
    yyaxis left;
    plot(tbltemp.Day, tbltemp.MaxTemperature, 'k.-', 'MarkerFaceColor', 'w', 'MarkerSize', 20, 'LineWidth', 1); 
    idx_fever = find(tbltemp.MaxTemperature>=100.4);
    plot(tbltemp.Day(idx_fever), tbltemp.MaxTemperature(idx_fever), 'r.', 'MarkerSize', 20);
    set(gca,'Xtick',xtick_lb:5:xtick_ub);
    
    ytick_lb = floor(min(tbltemp.MaxTemperature));
    ytick_ub = ceil(max(tbltemp.MaxTemperature));
    if (mod(ytick_lb,2)~=0)
        ytick_lb = ytick_lb-1;
    end
    if (mod(ytick_ub,2)~=0)
        ytick_ub = ytick_ub+1;
    end
    ylim([ytick_lb ytick_ub]);
    set(gca,'Ytick',[ytick_lb,ytick_lb+(ytick_ub-ytick_lb)/2,ytick_ub]);
    
    % plot reference line for fever
    yyaxis right;
    plot([window2plot(1), window2plot(2)], [100.4, 100.4], 'k--');
    ylim([ytick_lb ytick_ub]);
    set(gca,'Ytick',100.4);
    
    % set black axis
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    yyaxis left;
else
    set(gca,'YTick', []);
    set(gca, 'YTicklabel', []);
    set(gca,'XTick', []);
end

xlabel('');
xlim(window2plot);
set(gca, 'XTicklabel', []);
ylabel('Temperature (F)');
set(gca,'FontSize',12);

end

