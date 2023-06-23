
clear all 
clc 
PatientID = 'DBSTRD002'; 
load(sprintf('%s_PSD_tstats_25-May-2022.mat',PatientID));
load(sprintf('%s_zscore_bl_poststim_left_SCC_VCVS_052522.mat',PatientID)); 
mult_compare = load(sprintf('%s_mult_compare_data.mat',PatientID));

%% 
num_lSCC_trials = length(labels_lSCC);
num_lVCVS_trials = length(labels_lVCVS); 
total_trials = num_lSCC_trials + num_lVCVS_trials; 

% ACC alpha
% l_ACC_SCC_alpha= diff_data_left(1:num_lSCC_trials,4,3); 
% l_ACC_VCVS_alpha = diff_data_left(num_lSCC_trials+1:total_trials,4,3); 
% ACC_label_SCC_alpha = string(repmat('ACC_alpha',num_lSCC_trials,1)); 
% ACC_label_VCVS_alpha = string(repmat('ACC_alpha',num_lVCVS_trials,1)); 
% % ACC LG 
% l_ACC_SCC_lg = diff_data_left(1:num_lSCC_trials,4,5); 
% l_ACC_VCVS_lg = diff_data_left(num_lSCC_trials+1:total_trials,4,5); 
% ACC_label_SCC_lg = string(repmat('ACC_lg',num_lSCC_trials,1)); 
% ACC_label_VCVS_lg = string(repmat('ACC_lg',num_lVCVS_trials,1)); 
% % ACC HG 
% l_ACC_SCC_hg = diff_data_left(1:num_lSCC_trials,4,6); 
% l_ACC_VCVS_hg = diff_data_left(num_lSCC_trials+1:total_trials,4,6); 
% ACC_label_SCC_hg = string(repmat('ACC_hg',num_lSCC_trials,1)); 
% ACC_label_VCVS_hg = string(repmat('ACC_hg',num_lVCVS_trials,1)); 

%DPF LG 
l_DPF_SCC_lg = diff_data_left(1:num_lSCC_trials,6,5); 
l_DPF_VCVS_lg = diff_data_left(num_lSCC_trials+1:total_trials,6,5); 
DPF_label_SCC_lg = string(repmat('DPF_lg',num_lSCC_trials,1)); 
DPF_label_VCVS_lg = string(repmat('DPF_lg',num_lVCVS_trials,1)); 
%DPF HG 
l_DPF_SCC_hg = diff_data_left(1:num_lSCC_trials,6,6); 
l_DPF_VCVS_hg = diff_data_left(num_lSCC_trials+1:total_trials,6,6); 
DPF_label_SCC_hg = string(repmat('DPF_hg',num_lSCC_trials,1)); 
DPF_label_VCVS_hg = string(repmat('DPF_hg',num_lVCVS_trials,1)); 

%MOF HG 
l_MOF_SCC_hg = diff_data_left(1:num_lSCC_trials,1,6); 
l_MOF_VCVS_hg = diff_data_left(num_lSCC_trials+1:total_trials,1,6); 
MOF_label_SCC_hg = string(repmat('MOF_hg',num_lSCC_trials,1)); 
MOF_label_VCVS_hg = string(repmat('MOF_hg',num_lVCVS_trials,1)); 



%% concatenate all data and labels. 
all_left_SCC = vertcat(l_DPF_SCC_lg,l_DPF_SCC_hg,l_MOF_SCC_hg); 
all_left_VCVS = vertcat(l_DPF_VCVS_lg,l_DPF_VCVS_hg,l_MOF_VCVS_hg); 
all_left_matrix = vertcat(all_left_SCC,all_left_VCVS); 

DBS_target_SCC = string(repmat('SCC',length(all_left_SCC),1)); 
DBS_target_VCVS = string(repmat('VCVS',length(all_left_VCVS),1)); 

foiroi_label_SCC = vertcat(DPF_label_SCC_lg,DPF_label_SCC_hg,MOF_label_SCC_hg); 
foiroi_label_VCVS = vertcat(DPF_label_VCVS_lg,DPF_label_VCVS_hg,MOF_label_VCVS_hg); 

%% 
tbl_left = array2table(all_left_matrix, 'VariableNames',{'PSD'}); 
tbl_left{:,2} = vertcat(DBS_target_SCC,DBS_target_VCVS);
tbl_left{:,3} = vertcat(foiroi_label_SCC,foiroi_label_VCVS); 
tbl_left{:,4} = string(repmat('L',height(tbl_left),1)); 
tbl_left.Properties.VariableNames = {'PSD','DBS_target','foi_roi','Hemi'}; 


%% Save for python 
addpath(genpath('/Users/anushaallawala/Python_projects/'))
filename = sprintf('/Users/anushaallawala/Python_projects/%s_left_boxplot_data_%s.csv',PatientID,date); 

writetable(tbl_left,filename)



