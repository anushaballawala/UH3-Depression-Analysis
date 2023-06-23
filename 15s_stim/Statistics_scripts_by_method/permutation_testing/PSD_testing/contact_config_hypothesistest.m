%% contact_config_hypothesistest.m 

%%% FXN: Hypothesis testing to determine significant differences between
%%% contact configurations within a DBS target (SCC or VCVS) using ANOVA 

%%%% 
%% Load file that contains data for each trial 

clc 
close all 
clear
PatientID = 'DBSTRD002'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
%data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 
%Load channels labels that will be used later. 
load(sprintf('%s_goodch.mat',PatientID)) 

%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_elecconfig_conds_%s_new_leftright.mat',PatientID));
comp_labels = comp_labels.comps;

%% z-score data 

% calculate zscore for all values except for stim1,2,3

%get data that excludes stim (artifact) 
for i = 1:num_trials
    
    is_stim = contains(tbl.Stim_Labels(i), 'stim1') == 1 || contains(tbl.Stim_Labels(i), 'stim2') == 1 ||contains(tbl.Stim_Labels(i), 'stim3') == 1;  
    if is_stim == 1
        nostim_idx(i) = 0;
    elseif is_stim == 0
        nostim_idx(i) = 1;
    end
end

data_nostim = data(logical(nostim_idx'),:,:); 


%% relabel data after stim trials thrown out from zscoring 

comp_labels_old = comp_labels; 

for i = 1:numel(comp_labels)
    
    comp_labels{1,i}.label_vector = comp_labels{1,i}.label_vector(logical(nostim_idx'));  
end 

labels = {}; 
for i = 1:numel(comp_labels)
    labels{i} = comp_labels{i}.label_vector;
end 

tbl_new = tbl(logical(nostim_idx'),{'Stim_Labels','DBS_target','Elec_label','Frequency'}); 
%% Calculate z-score 

data_zscore = zeros(size(data_nostim)); 
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - mean(data_nostim(:,j,k)))./std(data_nostim(:,j,k)); 
        end 
    end 
end 

%% Define which type of data, Channel or ROI 

data_type_list = {'ROI', 'Ch'}; 
[idx_data_list,tf_data] = listdlg('ListString', data_type_list, 'ListSize',[150,250]); 
data_list_config = data_type_list{idx_data_list}; 
disp(data_list_config)
data_dim = data_list_config; 
 
if strcmp(data_dim,'ROI')==1 
    % get data that will be input into permutations fxn. 
    [matrix_ROI,ROI_labels] = generate_ROI_from_ch(data_zscore,PatientID,2);
end 

%%

lSCC_labels = comp_labels{1, 1}.label_vector  ;
rSCC_labels = comp_labels{1, 2}.label_vector  ;

lVCVS_labels = comp_labels{1, 3}.label_vector  ;
rVCVS_labels = comp_labels{1, 4}.label_vector  ;

alldata_lSCC = matrix_ROI(lSCC_labels>0,:,:); 
alldata_lVCVS = matrix_ROI(lVCVS_labels>0,:,:); 

alldata_rSCC = matrix_ROI(rSCC_labels>0,:,:); 
alldata_rVCVS = matrix_ROI(rVCVS_labels>0,:,:); 

alldata_baseline = matrix_ROI(lSCC_labels==0,:,:); 

idx_lSCC_labels = lSCC_labels>0; 
idx_lVCVS_labels = lVCVS_labels>0; 
idx_rSCC_labels = rSCC_labels>0; 
idx_rVCVS_labels = rVCVS_labels>0;


diff_lSCC_labels = lSCC_labels(idx_lSCC_labels); 
diff_lVCVS_labels = lVCVS_labels(idx_lVCVS_labels);
diff_lSCC = alldata_lSCC - mean(alldata_baseline,1); 
diff_lVCVS = alldata_lVCVS - mean(alldata_baseline,1); 

diff_rSCC_labels = rSCC_labels(idx_rSCC_labels); 
diff_rVCVS_labels = rVCVS_labels(idx_rVCVS_labels);
diff_rSCC = alldata_rSCC - mean(alldata_baseline,1); 
diff_rVCVS = alldata_rVCVS - mean(alldata_baseline,1); 


%% Run ANOVA for each DBS target LSCC

%First SCC on zscored data subtracted by baseline 
 
    for i = 1:size(diff_lSCC,2)
        for j = 1:size(diff_lSCC,3)          
        y_diff_lSCC = diff_lSCC(:,i,j); 
        group_diff_lSCC = diff_lSCC_labels; 
        [p_diff_lSCC(i,j),tbl_lSCC{i,j}] = anova1(y_diff_lSCC,group_diff_lSCC);
        close
        close all hidden
        end
    end 

 %%    rSCC
    
    for i = 1:size(diff_rSCC,2)
        for j = 1:size(diff_rSCC,3)          
        y_diff_rSCC = diff_rSCC(:,i,j); 
        group_diff_rSCC = diff_rSCC_labels; 
        [p_diff_rSCC(i,j),tbl_rSCC{i,j}] = anova1(y_diff_rSCC,group_diff_rSCC);
        close
        close all hidden
        end
    end 
    
%% Run for lVCVS 
 
    for i = 1:size(matrix_ROI,2)
        for j = 1:size(matrix_ROI,3)          
        y_diff_lVCVS = diff_lVCVS(:,i,j); 
        groupdiff_lVCVS = diff_lVCVS_labels; 
        [p_diff_lVCVS(i,j),tbl_lVCVS{i,j}] = anova1(y_diff_lVCVS,groupdiff_lVCVS);
        close
        close all hidden
        end
    end 

 %% rVCVS
 
    for i = 1:size(matrix_ROI,2)
        for j = 1:size(matrix_ROI,3)          
        y_diff_rVCVS = diff_rVCVS(:,i,j); 
        groupdiff_rVCVS = diff_rVCVS_labels; 
        [p_diff_rVCVS(i,j),tbl_rVCVS{i,j}] = anova1(y_diff_rVCVS,groupdiff_rVCVS);
        close
        close all hidden
        end
    end  
    
    %% save stat results 
    
filename = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/%s_anova_results_%s',PatientID,date); 
save(filename,'p_diff_lSCC','p_diff_rSCC','p_diff_lVCVS','p_diff_rVCVS',...
    'tbl_lSCC','tbl_rSCC','tbl_lVCVS','tbl_rVCVS','ROI_labels'); 


%% Run permutations and shuffle labels 

num_perm = 1000;                                                                                                                                                                                 
% t_stat structure is num_perm+1 X data_dim X num_freq. 
labels = {}; 
labels{1} = SCC_labels; 
labels{2} = VCVS_labels; 

%shuffled_labels = cell(numel(comp_labels),1); 

%permutations done separately for each condition being compared, just put
%in one cell. 

[p_val_anova_SCC, shuffled_labels_SCC] = permutations(alldata_SCC,diff_SCC_labels,num_perm,@compute_f_statistic_anova); 

%% Run permutation for VCVS 

[p_val_anova_VCVS, shuffled_labels_VCVS] = permutations(alldata_VCVS,diff_VCVS_labels,num_perm,@compute_f_statistic_anova); 
% 
%% 
% calculate p_value by looking at comparison of t-stat values and dividing
% by total number of permuations 

% Taking absolute value to check against both extreme ends of sample
% population. 

p_value_perm_VCVS = squeeze(mean(-log(p_val_anova_VCVS)>=-log(p_val_anova_VCVS(1,:,:)))); 

p_value_perm_SCC = squeeze(mean(-log(p_val_anova_SCC)>=-log(p_val_anova_SCC(1,:,:)))); 

%% Make heatmap of SCC p values that is ROI X FOI 

FOI_labels = {'delta', 'theta', 'alpha', 'beta', 'lowgamma','highgamma'}; 
figure('Position',[500 300 800 600]) 
h = heatmap(p_value_perm_VCVS); 
colormap redbluecmap(4)
co.YDisplayLabels = 'p values';
h.XDisplayLabels = FOI_labels; 
h.YDisplayLabels = ROI_labels; 
h.ColorLimits = [0.04 0.08]; 
title('VCVS poststim') 
% filename = 'lSCC_stim.png'; 
% saveas(gcf, filename)
% %% for 002 - match up the order for 001 
% p_value_perm_VCVS_new(1,:) = p_value_perm_VCVS(4,:); 
% p_value_perm_VCVS_new(2,:) = p_value_perm_VCVS(5,:); 
% p_value_perm_VCVS_new(3,:) = p_value_perm_VCVS(6,:); 
% p_value_perm_VCVS_new(4,:) = p_value_perm_VCVS(2,:); 
% p_value_perm_VCVS_new(5,:) = p_value_perm_VCVS(1,:); 
% p_value_perm_VCVS_new(6,:) = p_value_perm_VCVS(7,:); 
% p_value_perm_VCVS_new(7,:) = p_value_perm_VCVS(8,:); 
% p_value_perm_VCVS_new(8,:) = p_value_perm_VCVS(3,:); 
% p_value_perm_VCVS_new(9,:) = p_value_perm_VCVS(9,:); 
% p_value_perm_VCVS_new(10,:) = p_value_perm_VCVS(10,:); 
% %% 
% p_value_perm_SCC_new(1,:) = p_value_perm_SCC(4,:); 
% p_value_perm_SCC_new(2,:) = p_value_perm_SCC(5,:); 
% p_value_perm_SCC_new(3,:) = p_value_perm_SCC(6,:); 
% p_value_perm_SCC_new(4,:) = p_value_perm_SCC(2,:); 
% p_value_perm_SCC_new(5,:) = p_value_perm_SCC(1,:); 
% p_value_perm_SCC_new(6,:) = p_value_perm_SCC(7,:); 
% p_value_perm_SCC_new(7,:) = p_value_perm_SCC(8,:); 
% p_value_perm_SCC_new(8,:) = p_value_perm_SCC(3,:); 
% p_value_perm_SCC_new(9,:) = p_value_perm_SCC(9,:); 
% p_value_perm_SCC_new(10,:) = p_value_perm_SCC(10,:); 
% 
% %% ROI-labels match up 002 to 001 
% 
% ROI_labels_new = {'ACC','AMY','DPF','LOF','MOF','STG','SFG','VPF','PHG','EC'}; 
