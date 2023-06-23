%% DBStarget_hypothesistest.m 

%%% FXN: Permutation and hypothesis testing for all combinations of 
%%%      15-s experiments in 001 and 002. 

%%% ***Bug-fixes remaining: 
%%% 1. check whether individual current configurations can be tested 
%%% 2. check if lead  hemisphere data can be combined 
%%% 3. test poststim 
%%% 4. test stim on 
%%% 5. Save list of channels that are significant 
%%% 6. save experiment tested in metadata. 

%% Load file that contains data for each trial

clc 
close all 
clear
PatientID = 'DBSTRD003'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
%data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 

%Load channels labels that will be used later. 
ROI_labels = load(sprintf('%s_ROI_labels.mat',PatientID));
ch_labels = ROI_labels.good_ch_labels;
label_tbl = ROI_labels.summary_ch_info;
data(:,ROI_labels.goodch_idx_remove,:) = [];

%% **** Fix this in the preprocessing code, temporary code for now ***** 
% switch PatientID 
%     case 'DBSTRD001'
%         disp('change nothing')                                                                                                                                               
%     case 'DBSTRD002'
%         tmp_idx = [1:57,59:66,68:101,103:129];
%         data = data(:,tmp_idx,:); 
%     case 'DBSTRD003'
%         keep_stim_idx = [1:5,7:28,30:48,50,52:66];
%         good_ch_labels = good_ch_labels(keep_stim_idx); 
%         good_ch_idx = good_ch_idx(keep_stim_idx); 
% end 
% 
% 
% %58, 67(LOF,rLOF), 102 (deleted from rAMY)
%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
comp_labels = comp_labels.comps;

switch PatientID
    case 'DBSTRD001'
        stim_idx_lSCC = find(comp_labels{1, 9}.label_vector==1);
        stim_idx_rSCC = find(comp_labels{1, 17}.label_vector == 1);
        stim_idx_lVCVS = find(comp_labels{1, 31}.label_vector == 1);
        stim_idx_rVCVS = find(comp_labels{1,37}.label_vector == 1);
        baseline_idx = find(comp_labels{1, 1}.label_vector==0);        
        % Get vector of labels
        
        labels = {};
        labels{1} = comp_labels{9}.label_vector;
        labels{2} = comp_labels{17}.label_vector;
        labels{3} = comp_labels{31}.label_vector;
        labels{4} = comp_labels{37}.label_vector;
    case 'DBSTRD002'
        
        stim_idx_lSCC = find(comp_labels{1, 9}.label_vector==1);
        stim_idx_rSCC = find(comp_labels{1, 17}.label_vector == 1);
        stim_idx_lVCVS = find(comp_labels{1, 33}.label_vector == 1);
        stim_idx_rVCVS = find(comp_labels{1,41}.label_vector == 1); 
        baseline_idx = find(comp_labels{1,9}.label_vector == 0); 
        labels = {};
        labels{1} = comp_labels{9}.label_vector;
        labels{2} = comp_labels{17}.label_vector;
        labels{3} = comp_labels{33}.label_vector;
        labels{4} = comp_labels{41}.label_vector;
        
        
    
    case 'DBSTRD003'
        
        stim_idx_lSCC = find(comp_labels{1, 4}.label_vector==1);
        stim_idx_rSCC = find(comp_labels{1, 7}.label_vector == 1);
        stim_idx_lVCVS = find(comp_labels{1, 13}.label_vector == 1);
        stim_idx_rVCVS = find(comp_labels{1, 16}.label_vector == 1); 
        baseline_idx = find(comp_labels{1, 4}.label_vector == 0); 
        labels = {};
        labels{1} = comp_labels{4}.label_vector;
        labels{2} = comp_labels{7}.label_vector;
        labels{3} = comp_labels{13}.label_vector;
        labels{4} = comp_labels{16}.label_vector;
        
end
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

idx = logical(nostim_idx'); 

new_tbl = tbl(idx,:); 
    
%% relabel data after stim trials thrown out from zscoring 

    for i = 1:numel(labels)
        labels{i} = labels{i}(logical(nostim_idx'));  
    end 


%% Calculate z-score 

data_zscore = zeros(size(data_nostim)); 
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - nanmean(data_nostim(:,j,k)))./nanstd(data_nostim(:,j,k)); 
        end 
    end 
end 


%% Define which type of data, Channel or ROI 

data_type_list = {'ROI', 'Ch'}; 
[idx_data_list,tf_data] = listdlg('ListString', data_type_list, 'ListSize',[150,250]); 
data_list_config = data_type_list{idx_data_list}; 
disp(data_list_config)
data_dim = data_list_config; 
 
if strcmp(PatientID,'DBSTRD003')==1
    load('DBSTRD003_referencedsig_SCC_f130_234.mat','metadata')
end 
if strcmp(data_dim,'ROI')==1 
    % get data that will be input into permutations fxn. 
    [matrix_ROI,ROI_labels1] = generate_ROI_from_ch(data_zscore,PatientID,2);
end 

%% Subtract from baseline. 


baseline_data = matrix_ROI(labels{1}==0,:,:); 
baseline_cond = new_tbl(labels{1}==0,:); 
mean_bl = mean(baseline_data,1); 

% SCC. 
lSCC_data = matrix_ROI(labels{1}==1,:,:); 
rSCC_data = matrix_ROI(labels{2}==1,:,:); 
%SCC conds labels. 
lSCC_cond = new_tbl(labels{1}==1,:); 
rSCC_cond = new_tbl(labels{2}==1,:); 

% VCVS. 
lVCVS_data = matrix_ROI(labels{3}==1,:,:); 
rVCVS_data = matrix_ROI(labels{4}==1,:,:); 
% VCVS conds labels. 
lVCVS_cond = new_tbl(labels{3} ==1,:); 
rVCVS_cond = new_tbl(labels{4} ==1,:); 

%% 

diff_lSCC = lSCC_data - mean_bl; 
diff_rSCC = rSCC_data - mean_bl;

diff_lVCVS = lVCVS_data - mean_bl;
diff_rVCVS = rVCVS_data - mean_bl; 

%% New matrix with just SCC and VCVS (for comparison) . 

%Compare lVCVS and lSCC 
num_lSCC_trials = size(diff_lSCC,1); 
num_lVCVS_trials = size(diff_lVCVS,1); 
labels_lSCC = ones(num_lSCC_trials,1);
labels_lVCVS = zeros(num_lVCVS_trials,1); 
matrix_ROI_left = vertcat(diff_lSCC,diff_lVCVS); 
labels_left = vertcat(labels_lSCC,labels_lVCVS); 
assert(size(matrix_ROI_left,1)==length(labels_left)); 

%Compare rSCC and rVCVS 
num_rSCC_trials = size(diff_rSCC,1); 
num_rVCVS_trials = size(diff_rVCVS,1); 
labels_rSCC = ones(num_rSCC_trials,1);
labels_rVCVS = zeros(num_rVCVS_trials,1); 
matrix_ROI_right = vertcat(diff_rSCC,diff_rVCVS); 
labels_right = vertcat(labels_rSCC,labels_rVCVS); 
assert(size(matrix_ROI_right,1)==length(labels_right)); 


%% Concatenate condition tables and SCC/VCVS data. 

left_cond = vertcat(lSCC_cond,lVCVS_cond);
right_cond = vertcat(rSCC_cond,rVCVS_cond); 


%% Save data for python viz of current steering. 

currsteer_fname_l = sprintf('/Users/anushaallawala/Python_projects/%s_currsteer_data_viz_%s.mat',PatientID,date); 
save(currsteer_fname, 'matrix_ROI_left'); 

%% 

tmp = matrix_ROI_left(:,1,1); 

left_cond{:,5} = tmp; 
left_cond.Properties.VariableNames{5} = 'PSD';  
%% 

tmp = matrix_ROI_left(1:num_lSCC_trials,1,1); 
lSCC_cond{:,5} = tmp; 
lSCC_cond.Properties.VariableNames{5} = 'PSD';
writetable(lSCC_cond,'/Users/anushaallawala/Python_projects/test_currsteer_lSCC.csv'); 

%% 
currsteer_fname_r = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/%s_currsteer_PSD_%s.mat',PatientID,date); 
save(currsteer_fname, 'lSCC_cond','lVCVS_cond','rSCC_cond','rVCVS_cond',...
    'left_cond','right_cond','matrix_ROI_left','matrix_ROI_right','ROI_labels',...
    'ROI_labels1')


%% Permutation and shuffle labels 

num_perm = 1000;                                                                                                                                                                                 
% t_stat structure is num_perm+1 X data_dim X num_freq. 

%shuffled_labels = cell(numel(comp_labels),1); 

%permutations done separately for each condition being compared, just put
%in one cell. 
%for i = 1:numel(comp_labels)
    switch data_dim 
        case 'Ch'
            [t_stat,shuffled_labels] = permutations(data_zscore,labels,num_perm,@Tfunc); 
        case 'ROI' 
            [t_stat_left, shuffled_labels_left] = permutations(matrix_ROI_left,labels_left,num_perm,@Tfunc); 
            [t_stat_right, shuffled_labels_right] = permutations(matrix_ROI_right,labels_right,num_perm,@Tfunc); 

    end 

%end 




%% 
p_value_left = (squeeze(mean(abs(t_stat_left)>=abs(t_stat_left(1,:,:))))); 
p_value_right = (squeeze(mean(abs(t_stat_right)>=abs(t_stat_right(1,:,:))))); 

%% Save data 

filename = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/%s_PSD_tstats_%s.mat',PatientID,date); 
save(filename,'ch_labels','comp_labels','data_zscore','label_tbl','labels','labels_left',...
    'labels_right','labels_lSCC','labels_lVCVS','labels_right','labels_rSCC',...
    'labels_rVCVS','lSCC_data','lVCVS_data','matrix_ROI','matrix_ROI_left',...
    'matrix_ROI_right','p_value_left','p_value_right','ROI_labels','ROI_labels1','shuffled_labels_left',...
    'shuffled_labels_right','t_stat_left','t_stat_right','tbl'); 


% %% Save data for python viz. 
% diff_data_left = cat(1,diff_lSCC,diff_lVCVS); 
% py_filename_left = sprintf('/Users/anushaallawala/Python_projects/%s_zscore_bl_poststim_left_SCC_VCVS_052522.mat',PatientID); 
% save(py_filename_left,'diff_data_left'); 
% 
% diff_data_right = cat(1,diff_rSCC,diff_rVCVS); 
% py_filename_right = sprintf('/Users/anushaallawala/Python_projects/%s_zscore_bl_poststim_right_SCC_VCVS_052522.mat',PatientID); 
% save(py_filename_right,'diff_data_right'); 







