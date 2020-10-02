function plot_community(PatientID2Plot, RelativeTimePeriod2Plot, BloodCellType2Plot, tblcounts, tbltaxonomy, tblbc, tbldrug, tbltemp, tblinfection)

figure();
title(['patientID: ', PatientID2Plot]);
set(gca, 'FontName', 'Arial');

subplot(5,1,[4,5]);
hold on;
box on;
xtick_lb = round(RelativeTimePeriod2Plot(1)/5)*5;      % lower bound of xtick
xtick_ub = (round(RelativeTimePeriod2Plot(2)/5)+1)*5;  % upper bound of xtick
    
%% plot time series of microbiome composition and indicate infections if possible

if (~isempty(tblcounts))    
    % the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT 
    abundance_matrix = tblcounts{:, 4:end}; 

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
    if (height(tblcounts) == 1)
        h=bar([tblcounts.DayRelativeToNearestHCT,nan], [color_grouped_abundance; nan(1, size(color_grouped_abundance,2))], 'stacked', 'BarWidth', 0.1);
    else
        h=bar(tblcounts.DayRelativeToNearestHCT, color_grouped_abundance, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
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
        text(tblcounts.DayRelativeToNearestHCT(curr_row), (heightafter(curr_box)+heightbefore(curr_box))/2, dominant_tax, ...
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
        text(tblinfection.DayRelativeToNearestHCT(k), 1, '*', 'HorizontalAlignment', 'center', 'FontSize', 50);
    end
end
xlabel('Days relative to HCT date');
xlim(RelativeTimePeriod2Plot);
set(gca,'Xtick',[xtick_lb:5:xtick_ub]);
ylabel('16S relative abundance');
set(gca,'FontSize',12);

%% plot time series of blood cell counts of type "bc2plot"
subplot(5,1,3);
hold on;
box on;

if (~isempty(tblbc))
    
	% impute missing values in tblbc
    rows2insert = {};
    for i=RelativeTimePeriod2Plot(1):1:RelativeTimePeriod2Plot(2)
        if (sum(tblbc.DayRelativeToNearestHCT==i)==0)
            rows2insert = [rows2insert;{i,nan}];
        end
    end
    tblbc = [tblbc;rows2insert];
    tblbc = sortrows(tblbc, 'DayRelativeToNearestHCT');
    tblbc = fillmissing(tblbc,'linear');
    tblbc{tblbc.Value<0, 'Value'} = 0;
        
    yyaxis left;
    plot(tblbc.DayRelativeToNearestHCT, tblbc.Value, 'k.-', 'MarkerFaceColor', 'w', 'MarkerSize', 20, 'LineWidth', 1);
    if(strcmp(BloodCellType2Plot, 'Neutrophils'))
        idx_neutropenia = find(tblbc.Value <= 0.5); % threshold for neutropenia
        plot(tblbc.DayRelativeToNearestHCT(idx_neutropenia), tblbc.Value(idx_neutropenia), 'r.', 'MarkerSize', 20);
    end
    set(gca,'Xtick',xtick_lb:5:xtick_ub);
    
    ytick_lb = floor(min(tblbc.Value));
    ytick_ub = ceil(max(tblbc.Value));
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
    plot([RelativeTimePeriod2Plot(1), RelativeTimePeriod2Plot(2)], [0.5, 0.5], 'k--');
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
xlim(RelativeTimePeriod2Plot);
set(gca, 'XTicklabel', []);
if(strcmp(BloodCellType2Plot,'RBCtotal'))
    ylabel({BloodCellType2Plot, 'count (M/\muL)'});
else
    ylabel({BloodCellType2Plot, 'count (K/\muL)'});
end
set(gca,'FontSize',12);

%% plot time series of anti-infective drug administration
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
                administered_days = [administered_days,tbldrug_factor_i.StartDayRelativeToNearestHCT(j):1:tbldrug_factor_i.StopDayRelativeToNearestHCT(j)];              
                % plot a single period from first day to last day administered for factor i
                patch([tbldrug_factor_i.StartDayRelativeToNearestHCT(j)-boxwidth, tbldrug_factor_i.StopDayRelativeToNearestHCT(j)+boxwidth, tbldrug_factor_i.StopDayRelativeToNearestHCT(j)+boxwidth, tbldrug_factor_i.StartDayRelativeToNearestHCT(j)-boxwidth], ...
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
xlim(RelativeTimePeriod2Plot);
set(gca, 'XTicklabel', []);
ylabel('Medications');
set(gca,'Ytick',[]); 
set(gca,'Yticklabel',[]); 
set(gca,'FontSize',12);

%% plot time series of temperature
subplot(6,1,1);
hold on;
box on;

if (~isempty(tbltemp))
    yyaxis left;
    plot(tbltemp.DayRelativeToNearestHCT, tbltemp.MaxTemperature, 'k.-', 'MarkerFaceColor', 'w', 'MarkerSize', 20, 'LineWidth', 1); 
    idx_fever = find(tbltemp.MaxTemperature>=100.4);
    plot(tbltemp.DayRelativeToNearestHCT(idx_fever), tbltemp.MaxTemperature(idx_fever), 'r.', 'MarkerSize', 20);
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
    plot([RelativeTimePeriod2Plot(1), RelativeTimePeriod2Plot(2)], [100.4, 100.4], 'k--');
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
xlim(RelativeTimePeriod2Plot);
set(gca, 'XTicklabel', []);
ylabel('Temperature (F)');
set(gca,'FontSize',12);

end

