%% Chen Liao, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script creates X (independent variables) and Y (response variables) matrices used for linear regression (Y~X)
% last updated: April. 23, 2020

function [tblX, tblY, patientID_of_rows] = create_regression_matrices(data_path, use_pseudo_fraction, include_lastday_drug, drug_admin_route)

%% load sample table
opts = detectImportOptions(strcat(data_path, 'samples/tblASVsamples.csv'));
opts = setvartype(opts,{'PatientID'},'categorical');
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'),opts);
tblsamples.PatientID = categorical(tblsamples.PatientID);

%% load counts table
tblcounts = readtable(strcat(data_path, 'counts/tblASVcounts_human_filter.csv'));

%% load drug table and select for administration route
opts = detectImportOptions(strcat(data_path, 'meta_data/tbldrug.csv'));
opts = setvartype(opts,{'PatientID'},'categorical');
tbldrug = readtable(strcat(data_path, 'meta_data/tbldrug.csv'), opts);
tbldrug = tbldrug(ismember(tbldrug.Route,drug_admin_route), :);
if (height(tbldrug)==0)
    error('tbldrug is empty. check drug_admin_route.');
end

%% treat metronidazole, aztreonam as independent cateogry
special_factor_idx = find(ismember(tbldrug.Factor, {'metronidazole', 'aztreonam'}));
tbldrug.Category(special_factor_idx) = tbldrug.Factor(special_factor_idx);

%% replace ' ' and '-' with '_' in drug category
tbldrug.Category = strrep(tbldrug.Category, ' ', '_');
tbldrug.Category = strrep(tbldrug.Category, '-', '_');
tbldrug.Category = strrep(tbldrug.Category, '/', '_');

%% include only anti-bacterial drugs
drugcat2include = {'aminoglycosides', 'carbapenems', 'cephalosporins' ...
                   'glycopeptide_antibiotics', 'glycylcyclines', 'leprostatics', ...
                   'lincomycin_derivatives', 'macrolide_derivatives',...
                   'oxazolidinone_antibiotics', 'penicillins', 'quinolones',...
                   'sulfonamides', 'tetracyclines', 'metronidazole', 'aztreonam'};
tbldrug = tbldrug(ismember(tbldrug.Category, drugcat2include), :);

%% load qPCR table
tblqpcr = readtable(strcat(data_path, 'counts/tblqpcr.csv'));
tblqpcr = tblqpcr(tblqpcr.qPCR16S>0, :); % remove samples if its qPCR value = 0

%% join tblsamples, tblqpcr, and tblcounts
tbljoined = innerjoin(tblsamples, tblqpcr);
tbljoined = innerjoin(tbljoined, tblcounts);
tbljoined = tbljoined(:, {'PatientID','SampleID','Timepoint','qPCR16S','ASV','Count'});

%% unstack and calculate total count for each samples
tbljoined = unstack(tbljoined, 'Count', 'ASV');
ASV_start_idx = find(contains(tbljoined.Properties.VariableNames, 'ASV'), 1);
tbljoined.totalCount = nansum(tbljoined{:, ASV_start_idx:end},2); % ignore nan in sum
tbljoined = movevars(tbljoined, 'totalCount', 'After', 'qPCR16S');

%% normalize ASV counts to relative abundance
ASV_start_idx = find(contains(tbljoined.Properties.VariableNames, 'ASV'), 1);
counts_matrix = tbljoined{:, ASV_start_idx:end};
counts_matrix(isnan(counts_matrix)) = 0;
counts_matrix = counts_matrix ./ sum(counts_matrix, 2);
tbljoined{:, ASV_start_idx:end} = counts_matrix;

%% remove ASVs that not not present in all samples
tbljoined = tbljoined(:, [true(1,ASV_start_idx-1), sum(tbljoined{:, ASV_start_idx:end}, 1) > 0]);

%% get unique patientIDs, ASVs, and antibiotic categories
unique_patientIDs = unique(tbljoined.PatientID);
unique_ASVs = tbljoined.Properties.VariableNames(ASV_start_idx:end);
unique_drugCategories = unique(tbldrug.Category);

%% construct regression matrices
total_consecutive_sample_pairs = height(tbljoined) - length(unique_patientIDs);
matX = zeros(total_consecutive_sample_pairs, length(unique_drugCategories)+1);
matY = zeros(total_consecutive_sample_pairs, length(unique_ASVs));
patientID_of_rows = cell(total_consecutive_sample_pairs, 1); % patientID associated with each row

local_iter=1;
for i=1:length(unique_patientIDs)
    
    patientID_i = unique_patientIDs(i); % get patient id
    patientID_i_samples = tbljoined(tbljoined.PatientID==patientID_i, :); % get all samples of the patient
    patientID_i_samples = sortrows(patientID_i_samples, 'Timepoint'); % sort samples based on sample collection date
    patientID_i_drug  = tbldrug(tbldrug.PatientID==patientID_i,:); % get antibiotic administration records for the patient
    
    for j=2:height(patientID_i_samples)
        day_1 = patientID_i_samples.Timepoint(j-1); % day of the first sample
        day_2 = patientID_i_samples.Timepoint(j); % day of the second sample
        if (day_1==day_2)
            patientID_i_samples(:,1:5)
            error('sample interval cannot be zero.');
        end
        matX(local_iter,1) = day_2-day_1;
        
        % add pseudo fraction to relative abundance data
        relative_abundance_1 = patientID_i_samples{j-1, ASV_start_idx:end}; % taxonomic composition of the first sample
        relative_abundance_2 = patientID_i_samples{j, ASV_start_idx:end}; % taxonomic composition of the second sample
        if (use_pseudo_fraction)
            for k=1:length(relative_abundance_2)
                if (relative_abundance_1(k)==0 && relative_abundance_2(k) ==0)
                    % do not add pseudo fraction if both are zero
                else
                    pseudo_fraction_1 = 1.0 / patientID_i_samples.totalCount(j-1);
                    pseudo_fraction_2 = 1.0 / patientID_i_samples.totalCount(j);
                    % do not add pseudo fraction for a taxon if the fraction is higher than the relative abundance of the same taxon in another sample
                    if (relative_abundance_1(k)==0 && (pseudo_fraction_1 < relative_abundance_2(k)))
                        relative_abundance_1(k) = pseudo_fraction_1;
                    end
                    if (relative_abundance_2(k)==0 && (pseudo_fraction_2 < relative_abundance_1(k)))
                        relative_abundance_2(k) = pseudo_fraction_2;
                    end
                end
            end
        end
        
        % convert relative abundance to absolute abundance
        absolute_abundance_1 = patientID_i_samples.qPCR16S(j-1)*relative_abundance_1;
        absolute_abundance_2 = patientID_i_samples.qPCR16S(j)*relative_abundance_2;
        
        % compute log-difference of two samples
        % do not worry if absolute_abundance_1 or absolute_abundance_2 is
        % zero; we will filter out non-finite number later
        matY(local_iter, :) = log(absolute_abundance_2)-log(absolute_abundance_1);
        
        % limit drug administrations to those occuring between the sample interval
        if (include_lastday_drug)
            patientID_i_int_j_drug = patientID_i_drug(patientID_i_drug.StartTimepoint<=day_2 & patientID_i_drug.StopTimepoint>=day_1,:);
            patientID_i_int_j_drug{patientID_i_int_j_drug.StartTimepoint<day_1,'StartTimepoint'} =  day_1;
            patientID_i_int_j_drug{patientID_i_int_j_drug.StopTimepoint>day_2,'StopTimepoint'} =  day_2;
        else
            patientID_i_int_j_drug = patientID_i_drug(patientID_i_drug.StartTimepoint<day_2 & patientID_i_drug.StopTimepoint>=day_1,:);
            patientID_i_int_j_drug{patientID_i_int_j_drug.StartTimepoint<day_1,'StartTimepoint'} =  day_1;
            patientID_i_int_j_drug{patientID_i_int_j_drug.StopTimepoint>=day_2,'StopTimepoint'} =  day_2-1;
        end

        % compute cumulative exposure (in unit days) to specific abtibiotic categories
        for k=1:length(unique_drugCategories)
            patientID_i_int_j_drug_cat_k = patientID_i_int_j_drug(strcmp(patientID_i_int_j_drug.Category, unique_drugCategories{k}),:);
            if (height(patientID_i_int_j_drug_cat_k)>0)
                assert(sum(patientID_i_int_j_drug_cat_k.StopTimepoint<patientID_i_int_j_drug_cat_k.StartTimepoint)==0);
                matX(local_iter, k+1) = sum(patientID_i_int_j_drug_cat_k.StopTimepoint-patientID_i_int_j_drug_cat_k.StartTimepoint + 1);
            else
                matX(local_iter, k+1) = 0;
            end
        end
        
        patientID_of_rows{local_iter} = patientID_i;
        local_iter = local_iter+1;
    end
end

%% convert X and Y matrices to tables
tblX = array2table(matX, 'VariableNames', {'growth', unique_drugCategories{:}});
tblY = array2table(matY, 'VariableNames', unique_ASVs);

end