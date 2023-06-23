%% Load file that contains data for each trial 
%clear 

%% Load file that contains data for each trial 
%clear 

clc 
close all 
clear
PatientID = 'DBSTRD001'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 
%Load channels labels that will be used later. 
load(sprintf('%s_goodch.mat',PatientID)) 

%% Load ROI info. 

ROI_info = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/Electrodes Info/ROI_labels_mat/%s_ROI_labels.mat',...
    PatientID)); 

switch PatientID 
    % for 001
    case 'DBSTRD001'
        bipolar_remove_idx = [8,22,25,32,38,42,46,52,65,73,82,96,108,115,116,123];
        table_remove_idx = [59,60,61,132:136]; 
        del_idx = horzcat(bipolar_remove_idx,table_remove_idx); 
        % for 002 
    case 'DBSTRD002' 
        bipolar_remove_idx = [9,19,34,40,46,67,80,94,99,104,119]; 
        table_remove_idx = [57,58,102]; 
        del_idx = horzcat(bipolar_remove_idx,table_remove_idx); 
end 
        

%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
%comp_labels = load(sprintf('labels_for_poststim_elecconfig_conds_%s_new',PatientID));
comp_labels = comp_labels.comps;

%% Get labels for poststim for std. 

stim_idx_lSCC = find(comp_labels{1,9}.label_vector==1);
stim_idx_rSCC = find(comp_labels{1,17}.label_vector==1); 


switch PatientID 
    case 'DBSTRD001' 
       stim_idx_lVCVS = find(comp_labels{1,31}.label_vector==1);
       stim_idx_rVCVS = find(comp_labels{1,37}.label_vector==1); 
    case 'DBSTRD002' 
        stim_idx_lVCVS = find(comp_labels{1,33}.label_vector==1);
        stim_idx_rVCVS = find(comp_labels{1,41}.label_vector==1); 
end 

baseline_idx = find(comp_labels{1,1}.label_vector==0);  

lSCC = data(stim_idx_lSCC, :,:);  % lSCC elec all, 
rSCC = data(stim_idx_rSCC, :,:);  % rSCC elec all, 

lVCVS = data(stim_idx_lVCVS, :, :); % lVCVS
rVCVS = data(stim_idx_rVCVS, :,:); % rVCVS

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
idx_nostim = logical(nostim_idx'); 
stim_tbl = tbl(idx_nostim,:); 

%% relabel data after stim trials thrown out from zscoring 

comp_labels_old = comp_labels; 

for i = 1:numel(comp_labels)
    
    comp_labels{1,i}.label_vector = comp_labels{1,i}.label_vector(logical(nostim_idx'));  
    
end 

% Get new labels with poststim and BL only. 
switch PatientID 
    case 'DBSTRD001'                 
        labels = {}; 
        labels{1} = comp_labels{1,9}.label_vector ; 
        labels{2} = comp_labels{1,17}.label_vector ; 
        labels{3} = comp_labels{1,31}.label_vector ; 
        labels{4} = comp_labels{1,37}.label_vector ; 
    case 'DBSTRD002' 
        labels = {}; 
        labels{1} = comp_labels{1,9}.label_vector ; 
        labels{2} = comp_labels{1,17}.label_vector ; 
        labels{3} = comp_labels{1,33}.label_vector ; 
        labels{4} = comp_labels{1,41}.label_vector ;     
end 

%% Calculate z-score 

data_zscore = zeros(size(data_nostim)); 
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - mean(data_nostim(:,j,k)))./std(data_nostim(:,j,k)); 
        end 
    end 
end 

%% Define channels that will be analyzed. 

ch_data = data_zscore; 
ch_data(:,del_idx,:) = []; 

%% Exclude ROIs that are wonky from bipolar re-ref. 

ch_labels = ROI_info.summary_ch_info; 

switch PatientID 
    case 'DBSTRD001' 
        ch_labels(bipolar_remove_idx,:) = [];
    case 'DBSTRD002' 
        ch_labels(bipolar_remove_idx,:) = []; 
end 
 

%% Combine and get ROIs of interest for analysis. 


% remove channels with gray matter and get indices for ROIs of interest. 
is_grey = contains(ch_labels.Matter, 'Grey');
is_grey_roi = @(roi) find(contains(ch_labels.ROI, roi) & is_grey); 
switch PatientID 
    case 'DBSTRD001' 
        acc_idx = is_grey_roi('ACC') ; 
        amy_idx = is_grey_roi( 'AMY');

    case 'DBSTRD002' 
        acc_idx = is_grey_roi( 'aCC') ; 
        amy_idx = is_grey_roi( 'Amy');

end 
        
vpf_idx = is_grey_roi( 'VPF') ; 
mof_idx = is_grey_roi( 'MOF') ; 
lof_idx = is_grey_roi( 'LOF'); 
dpf_idx = is_grey_roi( 'DPF') ;
%%

% get detailed info for each ROI from bigger table. 
acc_info = ch_labels(acc_idx,:); 
vpf_info = ch_labels(vpf_idx,:); 
mof_info = ch_labels(mof_idx,:); 
lof_info = ch_labels(lof_idx,:); 
amy_info = ch_labels(amy_idx,:); 

% indices of all ROIs of interst for paper. 

%% combine all rois. 

acc  = mean(ch_data(:,acc_idx,:),2); 
amy = mean(ch_data(:,amy_idx,:),2); 
mof = mean(ch_data(:,mof_idx,:),2); 
lof = mean(ch_data(:,lof_idx,:),2);
vpf = mean(ch_data(:,vpf_idx,:),2); 

matrix_roi = [acc,amy,mof,lof,vpf]; 
ROI_labels = {'ACC','AMY','MOF','LOF','VPF'}; 

%% get t-statistic between off and on for SCC and VCVS respectively. 

num_perm = 1000; % change to a 1000 after verified.                                                                                                                                                                                 
% t_stat structure is num_perm+1 X data_dim X num_freq. 

%shuffled_labels = cell(numel(comp_labels),1); 
shuffled_labels = cell(numel(labels),1);

%permutations done separately for each condition being compared, just put
%in one cell. 
for i = 1:numel(labels)
    [t_stat{i},shuffled_labels{i}] = permutations(matrix_roi,labels{i},num_perm,@Tfunc);    
end


%% Get num of observations. 

%num_ROI = size(t_stat_left,2); % number of regions of interest 
num_ROI = 5; 
num_FOI = 6; % number of spectral features of interest 
num_obs = num_ROI*num_FOI; 

%% restructure t-stat matrix to run exchangeHT 
% testing data from stimulation experiments in left hemisphere 

% lSCC
k = 4; 
t_stat_lSCC = permute(t_stat{1}, [2,3,1]); 
tmp_t_lSCC = reshape(t_stat_lSCC,[num_obs,1001]); 
[pvals_lSCC,~,~] = exchangeHT(abs(tmp_t_lSCC),0.05,false); 
p_vals_all_lSCC = reshape(pvals_lSCC(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

% rSCC 
t_stat_rSCC = permute(t_stat{2}, [2,3,1]); 
tmp_t_rSCC = reshape(t_stat_rSCC,[num_obs,1001]); 
[pvals_rSCC,~,~] = exchangeHT(abs(tmp_t_rSCC),0.05,false); 
p_vals_all_rSCC = reshape(pvals_rSCC(:,k),[num_ROI,num_FOI]); %grab p values from k == 4


%lVCVS 
t_stat_lVCVS = permute(t_stat{3}, [2,3,1]); 
tmp_t_lVCVS = reshape(t_stat_lVCVS,[num_obs,1001]); 
[pvals_lVCVS,~,~] = exchangeHT(abs(tmp_t_lVCVS),0.05,false); 
p_vals_all_lVCVS = reshape(pvals_lVCVS(:,k),[num_ROI,num_FOI]); %grab p values from k == 4


%rVCVS
t_stat_rVCVS = permute(t_stat{4}, [2,3,1]); 
tmp_t_rVCVS = reshape(t_stat_rVCVS,[num_obs,1001]); 
[pvals_rVCVS,~,~] = exchangeHT(abs(tmp_t_rVCVS),0.05,false); 
p_vals_all_rVCVS = reshape(pvals_rVCVS(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

%% 



%% 












