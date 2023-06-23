%% permandhypothesistest.m 

%% Load file that contains data for each trial 
%clear 

clc 
close all 
clear
PatientID = 'DBSTRD002'; 
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
labels = {}; 
labels{1} = comp_labels{1,9}.label_vector ; 
labels{2} = comp_labels{1,17}.label_vector ; 
labels{3} = comp_labels{1,31}.label_vector ; 
labels{4} = comp_labels{1,37}.label_vector ; 
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
switch PatientID 
    case 'DBSTRD001' 
        acc_idx = find(contains(ch_labels.ROI, 'ACC') == 1 & contains(ch_labels.Matter, 'Grey') == 1); 
        amy_idx = find(contains(ch_labels.ROI, 'AMY') == 1 & contains(ch_labels.Matter, 'Grey') == 1);

    case 'DBSTRD002' 
        acc_idx = find(contains(ch_labels.ROI, 'aCC') == 1 & contains(ch_labels.Matter, 'Grey') == 1); 
        amy_idx = find(contains(ch_labels.ROI, 'Amy') == 1 & contains(ch_labels.Matter, 'Grey') == 1);

end 
        
vpf_idx = find(contains(ch_labels.ROI, 'VPF') == 1 & contains(ch_labels.Matter, 'Grey') == 1); 
mof_idx = find(contains(ch_labels.ROI, 'MOF') == 1 & contains(ch_labels.Matter, 'Grey') == 1); 
lof_idx = find(contains(ch_labels.ROI, 'LOF') == 1 & contains(ch_labels.Matter, 'Grey') == 1); 
dpf_idx = find(contains(ch_labels.ROI, 'DPF') == 1 & contains(ch_labels.Matter, 'Grey') == 1);

% get detailed info for each ROI from bigger table. 
acc_info = ch_labels(acc_idx,:); 
vpf_info = ch_labels(vpf_idx,:); 
mof_info = ch_labels(mof_idx,:); 
lof_info = ch_labels(lof_idx,:); 
dpf_info = ch_labels(dpf_idx,:);
amy_info = ch_labels(amy_idx,:); 

% indices of all ROIs of interst for paper. 

%all_rois_idx = sort([acc_idx;vpf_idx;mof_idx;lof_idx],'ascend'); 
all_rois_idx = sort([acc_idx;vpf_idx;mof_idx;lof_idx;dpf_idx;amy_idx],'ascend'); 

ch_roi = ch_data(:,all_rois_idx,:); 

%% get t-statistic between off and on for SCC and VCVS respectively. 

num_perm = 1000; % change to a 1000 after verified.                                                                                                                                                                                 
% t_stat structure is num_perm+1 X data_dim X num_freq. 

%shuffled_labels = cell(numel(comp_labels),1); 
shuffled_labels = cell(numel(labels),1);

%permutations done separately for each condition being compared, just put
%in one cell. 
for i = 1:numel(labels)
    [t_stat{i},shuffled_labels{i}] = permutations(ch_roi,labels{i},num_perm,@Tfunc);    
end

t_stat_lSCC = squeeze(t_stat{1}(1,:,:)); 
t_stat_rSCC = squeeze(t_stat{2}(1,:,:)); 
t_stat_lVCVS = squeeze(t_stat{3}(1,:,:)); 
t_stat_rVCVS = squeeze(t_stat{4}(1,:,:)); 


%% viz. tstat values between off and on for VCVS and SCC. 

% figure() 
% histogram(t_stat_lSCC(:,6)); 
% hold on 
% histogram(t_stat_lVCVS(:,6)); 
% 
% 


%% Save data for python viz. 

tbl_delta = table(t_stat_lSCC(:,1),t_stat_lVCVS(:,1),'VariableNames',{'lSCC','lVCVS'}); 
tbl_theta = table(t_stat_lSCC(:,2),t_stat_lVCVS(:,2),'VariableNames',{'lSCC','lVCVS'}); 
tbl_alpha = table(t_stat_lSCC(:,3),t_stat_lVCVS(:,3),'VariableNames',{'lSCC','lVCVS'}); 
tbl_beta = table(t_stat_lSCC(:,4),t_stat_lVCVS(:,4),'VariableNames',{'lSCC','lVCVS'}); 
tbl_lowg = table(t_stat_lSCC(:,5),t_stat_lVCVS(:,5),'VariableNames',{'lSCC','lVCVS'}); 
tbl_highg = table(t_stat_lSCC(:,6),t_stat_lVCVS(:,6),'VariableNames',{'lSCC','lVCVS'}); 

%delta
py_filename_left_delta = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_left_delta.csv',PatientID,date); 
writetable(tbl_delta,py_filename_left_delta); 
%theta 
py_filename_left_theta = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_left_theta.csv',PatientID,date); 
writetable(tbl_theta,py_filename_left_theta); 
%alpha 
py_filename_left_alpha = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_left_alpha.csv',PatientID,date); 
writetable(tbl_alpha,py_filename_left_alpha); 
%beta 
py_filename_left_beta = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_left_beta.csv',PatientID,date); 
writetable(tbl_beta,py_filename_left_beta); 
%lowgamma
py_filename_left_lowg = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_left_lowg.csv',PatientID,date); 
writetable(tbl_lowg,py_filename_left_lowg); 
%highgamma
py_filename_left_highg = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_left_highg.csv',PatientID,date); 
writetable(tbl_highg,py_filename_left_highg); 

%% repeat for right 

tbl_delta = table(t_stat_rSCC(:,1),t_stat_rVCVS(:,1),'VariableNames',{'rSCC','rVCVS'}); 
tbl_theta = table(t_stat_rSCC(:,2),t_stat_rVCVS(:,2),'VariableNames',{'rSCC','rVCVS'}); 
tbl_alpha = table(t_stat_rSCC(:,3),t_stat_rVCVS(:,3),'VariableNames',{'rSCC','rVCVS'}); 
tbl_beta = table(t_stat_rSCC(:,4),t_stat_rVCVS(:,4),'VariableNames',{'rSCC','rVCVS'}); 
tbl_lowg = table(t_stat_rSCC(:,5),t_stat_rVCVS(:,5),'VariableNames',{'rSCC','rVCVS'}); 
tbl_highg = table(t_stat_rSCC(:,6),t_stat_rVCVS(:,6),'VariableNames',{'rSCC','rVCVS'}); 

%delta
py_filename_right_delta = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_right_delta.csv',PatientID,date); 
writetable(tbl_delta,py_filename_right_delta); 
%theta 
py_filename_right_theta = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_right_theta.csv',PatientID,date); 
writetable(tbl_theta,py_filename_right_theta); 
%alpha 
py_filename_right_alpha = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_right_alpha.csv',PatientID,date); 
writetable(tbl_alpha,py_filename_right_alpha); 
%beta 
py_filename_right_beta = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_right_beta.csv',PatientID,date); 
writetable(tbl_beta,py_filename_right_beta); 
%lowgamma
py_filename_right_lowg = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_right_lowg.csv',PatientID,date); 
writetable(tbl_lowg,py_filename_right_lowg); 
%highgamma
py_filename_right_highg = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_allroi_tstat_ch_right_highg.csv',PatientID,date); 
writetable(tbl_highg,py_filename_right_highg); 

%% save all vars. 

filename = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/%s_%s_allroi_tstat_ch.mat',PatientID,date); 
save(filename);
%% Save output from permutations 


% metadata.stats.labels = labels; 
% metadata.stats.num_perm = num_perm; 
% metadata.stats.tbl = tbl;
% metadata.stats.num_trials = num_trials; 
% metadata.stats.DBS_labels = DBS_labels; 
% 
% metadata.PatientID = PatientID; 
% % metadata.good_ch_idx = good_ch_idx; 
% % metadata.good_ch_labels = good_ch_labels; 
% %metadata.good_ch_note = note; 
% metadata.ROI.ROI_labels = ROI_labels;
% 
% %add hard drive to path. 
% addpath(genpath('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data')); 
% filepath = '/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data'; 
% filename = sprintf('%s/%s_t_stat_vals_ch_offon_%s_zscore.mat',filepath,PatientID,date); 
% save(filename,'data','data_zscore','data_nostim','t_stat','shuffled_labels', 'ROI_labels1','labels','shuffled_labels','comp_labels'); 
% 
% disp('t stat saved')
