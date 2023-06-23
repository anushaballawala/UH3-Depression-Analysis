
%% Load data 
clear all 
clc 
clear
PatientID = 'DBSTRD001'; 
load(sprintf('%s_PSD_tstats_25-May-2022.mat',PatientID));
load(sprintf('%s_zscore_bl_poststim_left_SCC_VCVS_052522',PatientID)); 
mult_compare = load(sprintf('%s_mult_compare_data.mat',PatientID));
%% 001 - left 


%create a table. 

%left SCC vs VCVS 
%headings: PSD, DBS_target, foi_roi_label

num_lSCC_trials = length(labels_lSCC);
num_lVCVS_trials = length(labels_lVCVS); 
total_trials = num_lSCC_trials + num_lVCVS_trials; 

% ACC low gamma
l_ACC_SCC = diff_data_left(1:num_lSCC_trials,1,5); 
l_ACC_VCVS = diff_data_left(num_lSCC_trials+1:total_trials,1,5); 
ACC_label_SCC = string(repmat('ACC_beta',num_lSCC_trials,1)); 
ACC_label_VCVS = string(repmat('ACC_beta',num_lVCVS_trials,1)); 

%DPF alpha 
l_DPF_SCC_alpha = diff_data_left(1:num_lSCC_trials,3,3); 
l_DPF_VCVS_alpha = diff_data_left(num_lSCC_trials+1:total_trials,3,3);
DPF_label_alpha_SCC = string(repmat('DPF_alpha',num_lSCC_trials,1)); 
DPF_label_alpha_VCVS = string(repmat('DPF_alpha',num_lVCVS_trials,1)); 

%DPF beta 
l_DPF_SCC_beta = diff_data_left(1:num_lSCC_trials,3,4); 
l_DPF_VCVS_beta = diff_data_left(num_lSCC_trials+1:total_trials,3,4);
DPF_label_beta_SCC = string(repmat('DPF_beta',num_lSCC_trials,1)); 
DPF_label_beta_VCVS = string(repmat('DPF_beta',num_lVCVS_trials,1)); 

% MOF delta 
l_MOF_SCC_delta = diff_data_left(1:num_lSCC_trials,5,1); 
l_MOF_VCVS_delta = diff_data_left(num_lSCC_trials+1:total_trials,5,1);
MOF_label_SCC_delta = string(repmat('MOF_delta',num_lSCC_trials,1)); 
MOF_label_VCVS_delta = string(repmat('MOF_delta',num_lVCVS_trials,1)); 

%MOF theta 
l_MOF_SCC_theta = diff_data_left(1:num_lSCC_trials,5,2); 
l_MOF_VCVS_theta = diff_data_left(num_lSCC_trials+1:total_trials,5,2);
MOF_label_SCC_theta = string(repmat('MOF_theta',num_lSCC_trials,1)); 
MOF_label_VCVS_theta = string(repmat('MOF_theta',num_lVCVS_trials,1)); 

% SFG alpha 
l_SFG_SCC_delta = diff_data_left(1:num_lSCC_trials,7,1);
l_SFG_VCVS_delta = diff_data_left(num_lSCC_trials+1:total_trials,7,1);
SFG_label_delta_SCC = string(repmat('SFG_delta',num_lSCC_trials,1)); 
SFG_label_delta_VCVS = string(repmat('SFG_delta',num_lVCVS_trials,1)); 

%SFG beta 
l_SFG_SCC_beta = diff_data_left(1:num_lSCC_trials,7,4);
l_SFG_VCVS_beta = diff_data_left(num_lSCC_trials+1:total_trials,7,4);
SFG_label_beta_SCC = string(repmat('SFG_beta',num_lSCC_trials,1)); 
SFG_label_beta_VCVS = string(repmat('SFG_beta',num_lVCVS_trials,1)); 

% SFG high gamma 
l_SFG_SCC_hg = diff_data_left(1:num_lSCC_trials,7,6);
l_SFG_VCVS_hg = diff_data_left(num_lSCC_trials+1:total_trials,7,6);
SFG_label_hg_SCC = string(repmat('SFG_hg',num_lSCC_trials,1)); 
SFG_label_hg_VCVS = string(repmat('SFG_hg',num_lVCVS_trials,1)); 


%MTG alpha 
l_MTG_SCC_alpha = diff_data_left(1:num_lSCC_trials,6,3);
l_MTG_VCVS_alpha = diff_data_left(num_lSCC_trials+1:total_trials,6,3);
MTG_label_alpha_SCC = string(repmat('MTG_alpha',num_lSCC_trials,1));
MTG_label_alpha_VCVS = string(repmat('MTG_alpha',num_lVCVS_trials,1)); 

% MTG beta 
l_MTG_SCC_beta = diff_data_left(1:num_lSCC_trials,6,4);
l_MTG_VCVS_beta = diff_data_left(num_lSCC_trials+1:total_trials,6,4);
MTG_label_beta_SCC= string(repmat('MTG_beta',num_lSCC_trials,1)); 
MTG_label_beta_VCVS= string(repmat('MTG_beta',num_lVCVS_trials,1)); 

all_left_SCC = vertcat(l_ACC_SCC,l_DPF_SCC_alpha,l_DPF_SCC_beta,l_MOF_SCC_delta,l_MOF_SCC_theta,l_SFG_SCC_delta,l_SFG_SCC_beta,l_SFG_SCC_hg,l_MTG_SCC_alpha,l_MTG_SCC_beta); 
all_left_VCVS = vertcat(l_ACC_VCVS,l_DPF_VCVS_alpha,l_DPF_VCVS_beta,l_MOF_VCVS_delta,l_MOF_VCVS_theta,l_SFG_VCVS_delta,l_SFG_VCVS_beta,l_SFG_VCVS_hg,l_MTG_VCVS_alpha,l_MTG_VCVS_beta); 
all_left_matrix = vertcat(all_left_SCC,all_left_VCVS); 

DBS_target_SCC = string(repmat('SCC',length(all_left_SCC),1)); 
DBS_target_VCVS = string(repmat('VCVS',length(all_left_VCVS),1)); 

foiroi_label_SCC = vertcat(ACC_label_SCC,DPF_label_alpha_SCC,DPF_label_beta_SCC,MOF_label_SCC_delta,MOF_label_SCC_theta,SFG_label_delta_SCC,...
    SFG_label_beta_SCC,SFG_label_hg_SCC,MTG_label_alpha_SCC,MTG_label_beta_SCC); 
foiroi_label_VCVS = vertcat(ACC_label_VCVS,DPF_label_alpha_VCVS,DPF_label_beta_VCVS,MOF_label_VCVS_delta,MOF_label_VCVS_theta,SFG_label_delta_VCVS,...
    SFG_label_beta_VCVS,SFG_label_hg_VCVS, MTG_label_alpha_VCVS,MTG_label_beta_VCVS); 

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



