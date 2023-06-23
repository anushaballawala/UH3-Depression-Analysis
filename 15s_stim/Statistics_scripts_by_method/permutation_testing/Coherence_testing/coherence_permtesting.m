%% permandhypothesistest.m 

%%% FXN: Permutation and hypothesis testing for all combinations of 
%%%      15-s experiments in 001 and 002. 

%% Load file that contains data for each trial 
%clear 

clc 
close all 
clear
PatientID = 'DBSTRD001'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
data_file = load(sprintf('%s_autocorrbl_data_for_coherence_stats.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 

ROI_labels = load(sprintf('%s_ROI_labels.mat',PatientID)); 
    ch_labels = ROI_labels.good_ch_labels; 
    label_tbl = ROI_labels.summary_ch_info; 
%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

% comp_labels = load(sprintf('coherence_labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
comp_labels = load(sprintf('coherence_labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
comp_labels = comp_labels.comps;

%% Get labels for poststim for std. 

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
        labels{4} = comp_labels{4}.label_vector;
        
end


%% Get data for each condition 

data_lSCC_diff = squeeze(mean(data(labels{1}==1,:,:,:),1)) - squeeze(mean(data(labels{1}==0,:,:,:),1)); 
data_rSCC_diff = squeeze(mean(data(labels{2}==1,:,:,:),1)) - squeeze(mean(data(labels{2}==0,:,:,:),1)); 
data_lVCVS_diff = squeeze(mean(data(labels{3}==1,:,:,:),1)) - squeeze(mean(data(labels{3}==0,:,:,:),1)); 
data_rVCVS_diff = squeeze(mean(data(labels{4}==1,:,:,:),1)) - squeeze(mean(data(labels{4}==0,:,:,:),1)); 

%% 
data_baseline = mean(data(labels{1}==0,:,:,:),1); 
data_lSCC = mean(data(labels{1}==1,:,:,:),1); 
data_rSCC = mean(data(labels{2}==1,:,:,:),1);
data_lVCVS = mean(data(labels{3}==1,:,:,:),1);
data_rVCVS  = mean(data(labels{4}==1,:,:,:),1);

%% Not z-scoring -- existing code might be incorrect. 

% data_zscore = zeros(size(data)); 
% for i = 1:size(data,1)
%     for j = 1:size(data,2)        
%         for k = 1:size(data,4)
%             data_zscore(i,j,j,k) = (data(i,j,j,k) - mean(data(:,j,j,k)))./std(data(:,j,j,k)); 
%         end 
%     end 
% end 

%% Define which type of data, Channel or ROI. 

data_type_list = {'ROI', 'Ch'}; 
[idx_data_list,tf_data] = listdlg('ListString', data_type_list, 'ListSize',[150,250]); 
data_list_config = data_type_list{idx_data_list}; 
disp(data_list_config)
data_dim = data_list_config; 
 
if strcmp(data_dim,'ROI')==1 
    % get data that will be input into permutations fxn. 
    [matrix_ROI,ROI_labels] = generate_ROI_from_ch(data,PatientID,4);
elseif strcmp(data_dim,'Ch') ==1
    ROI_labels = load(sprintf('%s_ROI_labels.mat',PatientID)); 
    matrix_ROI = data; 
    ch_labels = ROI_labels.good_ch_labels; 
    tbl = ROI_labels.summary_ch_info; 
end 

%% Run permutations and shuffle labels 

num_perm = 1000; 
num_ch = size(matrix_ROI,2); 
num_freq = size(matrix_ROI,4); 
% t_stat structure is num_perm+1 X data_dim X num_freq for each cell 
t_stat = cell(1,numel(labels)); 
shuffled_labels = cell(1,numel(labels)); 

parfor i = 1:numel(labels)
    [t_stat{i},shuffled_labels{i}] = permutations_coh(data,labels{i},num_perm,@Tfunc_coh);    
end

% ----------------------------------------------

%% PERFORM MULTIPLE COMPARISONS TESTING FIRST. 

%restructure t-stat matrix to run MT for ALL FOI and ROI
num_ch = size(data_lSCC_diff,1); 
num_foi = 6; 


t_stat_for_HT_all = permute(t_stat{1}, [2,3,4,1]); 

%----SEND TO MH 
tmp = reshape(t_stat_for_HT_all,[128*128*6,1001]); 

%reshape(t_stat, [48, 1001]); 

[pvals_all,~,~] = exchangeHT(abs(tmp),0.05,true); 

p_vals_all_FOI = reshape(pvals_all(:,4),[num_ch,num_ch,num_foi]); %grab p values from k == 3

%----- SEND TO MH 
%% Check that p-value in first k column is equivalent to what I had obtained previously after perm testing. 

%p-value from perm testing: 
p_value; 


% ----------------------------------------------
%% Convert t_stat into matrix of conditions 

% concatenate and make into matrix that is perm X ROI X freq X SCC conds X
% VCVS conds. 

t_stat_lSCC = t_stat{1}; 
t_stat_rSCC = t_stat{2}; 
t_stat_lVCVS = t_stat{3}; 
t_stat_rVCVS = t_stat{4}; 


%% Get p value 

p_value_lSCC = (squeeze(mean(abs(t_stat_lSCC)>=abs(t_stat_lSCC(1,:,:,:))))); 
p_value_rSCC = (squeeze(mean(abs(t_stat_rSCC)>=abs(t_stat_rSCC(1,:,:,:))))); 
p_value_lVCVS = (squeeze(mean(abs(t_stat_lVCVS)>=abs(t_stat_lVCVS(1,:,:,:))))); 
p_value_rVCVS = (squeeze(mean(abs(t_stat_rVCVS)>=abs(t_stat_rVCVS(1,:,:,:))))); 
%% Save output from permutations 

metadata.stats.labels = labels; 
metadata.stats.num_perm = num_perm; 
metadata.stats.tbl = tbl;
metadata.stats.num_trials = num_trials; 
metadata.stats.DBS_labels = DBS_labels; 

metadata.PatientID = PatientID; 
metadata.good_ch_idx = good_ch_idx; 
metadata.good_ch_labels = good_ch_labels; 
%metadata.good_ch_note = note; 
metadata.ROI.ROI_labels = ROI_labels;

%add hard drive to path. 
addpath(genpath('/Users/anushaallawala/Data/DBSTRD/Stats_data')); 
filepath = '/Users/anushaallawala/Data/DBSTRD/Stats_data'; 
filename = sprintf('%s/%s_t_stat_vals_%s_zscore.mat',filepath,PatientID,date); 
save(filename,'data','data_zscore','data_nostim','t_stat','shuffled_labels'); 

disp('t stat saved')
%% Convert t_stat into matrix of conditions 

% concatenate and make into matrix that is perm X ROI X freq X SCC conds X
% VCVS conds. 

t_stat_SCC_cell = {}; 
for i = 1:SCC_idx
t_stat_SCC_cell{i} = t_stat{i}; 
end 

t_stat_VCVS_cell = {};
for j = VCVS_idx:numel(comp_labels)
    t_stat_VCVS_cell{j} = t_stat{j}; 
end 
t_stat_VCVS_cell = t_stat_VCVS_cell(~cellfun('isempty',t_stat_VCVS_cell)); 

%% convert to matrix and concatenate

switch PatientID
    case 'DBSTRD001'
        t_stat_SCC = cat(4,t_stat_SCC_cell{:}); 
        t_stat_VCVS = cat(4,t_stat_VCVS_cell{:}); 
        tmp = zeros(size(t_stat_SCC)); 
        tmp(:,:,:,1:18) = t_stat_VCVS; 
        t_stat_VCVS = tmp; 
        %pad VCVS so that it can be concatenated 
        t_stat_mat = cat(5,t_stat_SCC,t_stat_VCVS); 
               
    case 'DBSTRD002'
        t_stat_SCC = cat(4,t_stat_SCC_cell{:}); 
        t_stat_VCVS = cat(4,t_stat_VCVS_cell{:});
        t_stat_mat = cat(5,t_stat_SCC,t_stat_VCVS); 
end 

% compute p values X ROI X freq X SCC conds X VCVS conds  


%% Hypothesis testing - get p-value 
% p values X ROI X freq X SCC conds X VCVS conds
% calculate p_value by looking at comparison of t-stat values and dividing
% by total number of permuations 

% Taking absolute value to check against both extreme ends of sample
% population. 
p_value = (squeeze(mean(abs(t_stat_mat)>=abs(t_stat_mat(1,:,:,:,:))))); 

%% save data for jupyter 

alldata_SCC_1 = alldata_SCC{1}; 
alldata_VCVS_1 = alldata_VCVS{1}; 

diff_SCC = alldata_SCC_1 - mean(alldata_baseline,1); 
diff_VCVS = alldata_VCVS_1 - mean(alldata_baseline,1); 
%% 
diff_data = cat(1,diff_SCC,diff_VCVS); 

py_filename = sprintf('/Users/anushaallawala/Python_projects/%s_zscore_bl_poststim_SCC_VCVS.mat',PatientID); 
save(py_filename,'diff_data'); 




