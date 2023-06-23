
clear all 
clc 
PatientID = 'DBSTRD002'; 
load(sprintf('%s_PSD_tstats_25-May-2022.mat',PatientID));
load(sprintf('%s_zscore_bl_poststim_right_SCC_VCVS_052522.mat',PatientID)); 
mult_compare = load(sprintf('%s_mult_compare_data.mat',PatientID));

%% 002 Right 

num_rSCC_trials = length(labels_rSCC);
num_rVCVS_trials = length(labels_rVCVS); 
total_trials = num_rSCC_trials + num_rVCVS_trials; 

%ACC delta 
r_ACC_SCC = diff_data_right(1:num_rSCC_trials,4,1); 
r_ACC_VCVS = diff_data_right(num_rSCC_trials+1:total_trials,4,1);
ACC_label_SCC = string(repmat('ACC_delta_r',num_rSCC_trials,1)); 
ACC_label_VCVS = string(repmat('ACC_delta_r',num_rVCVS_trials,1)); 

%MOF delta 
r_MOF_SCC_delta = diff_data_right(1:num_rSCC_trials,1,1);
r_MOF_VCVS_delta = diff_data_right(num_rSCC_trials+1:total_trials,1,1);
MOF_label_SCC_delta = string(repmat('MOF_delta_r',num_rSCC_trials,1)); 
MOF_label_VCVS_delta = string(repmat('MOF_delta_r',num_rVCVS_trials,1)); 

%MOF LG 
r_MOF_SCC_lg = diff_data_right(1:num_rSCC_trials,1,5);
r_MOF_VCVS_lg = diff_data_right(num_rSCC_trials+1:total_trials,1,5);
MOF_label_SCC_lg = string(repmat('MOF_lg_r',num_rSCC_trials,1)); 
MOF_label_VCVS_lg = string(repmat('MOF_lg_r',num_rVCVS_trials,1)); 

% SFG delta 
r_SFG_SCC_delta = diff_data_right(1:num_rSCC_trials,8,1);
r_SFG_VCVS_delta = diff_data_right(num_rSCC_trials+1:total_trials,8,1);
SFG_label_SCC_delta = string(repmat('SFG_delta_r',num_rSCC_trials,1)); 
SFG_label_VCVS_delta = string(repmat('SFG_delta_r',num_rVCVS_trials,1)); 

% STG delta 
r_STG_SCC_delta = diff_data_right(1:num_rSCC_trials,7,1);
r_STG_VCVS_delta = diff_data_right(num_rSCC_trials+1:total_trials,7,1);
STG_label_SCC_delta = string(repmat('STG_delta_r',num_rSCC_trials,1)); 
STG_label_VCVS_delta = string(repmat('STG_delta_r',num_rVCVS_trials,1));

%STG theta 
r_STG_SCC_theta = diff_data_right(1:num_rSCC_trials,7,2);
r_STG_VCVS_theta = diff_data_right(num_rSCC_trials+1:total_trials,7,2);
STG_label_SCC_theta = string(repmat('STG_theta_r',num_rSCC_trials,1)); 
STG_label_VCVS_theta = string(repmat('STG_theta_r',num_rVCVS_trials,1));
% STG alpha 
r_STG_SCC_alpha = diff_data_right(1:num_rSCC_trials,7,3);
r_STG_VCVS_alpha = diff_data_right(num_rSCC_trials+1:total_trials,7,3);
STG_label_SCC_alpha = string(repmat('STG_alpha_r',num_rSCC_trials,1)); 
STG_label_VCVS_alpha = string(repmat('STG_alpha_r',num_rVCVS_trials,1));
% STG beta 
r_STG_SCC_beta = diff_data_right(1:num_rSCC_trials,7,4);
r_STG_VCVS_beta = diff_data_right(num_rSCC_trials+1:total_trials,7,4);
STG_label_SCC_beta = string(repmat('STG_beta_r',num_rSCC_trials,1)); 
STG_label_VCVS_beta = string(repmat('STG_beta_r',num_rVCVS_trials,1));
% STG low gamma 
r_STG_SCC_lg = diff_data_right(1:num_rSCC_trials,7,5);
r_STG_VCVS_lg = diff_data_right(num_rSCC_trials+1:total_trials,7,5);
STG_label_SCC_lg = string(repmat('STG_lg_r',num_rSCC_trials,1)); 
STG_label_VCVS_lg = string(repmat('STG_lg_r',num_rVCVS_trials,1));

%% Concatenate. 

all_right_SCC = vertcat(r_ACC_SCC,r_MOF_SCC_delta,r_MOF_SCC_lg,r_SFG_SCC_delta,r_STG_SCC_delta,r_STG_SCC_alpha,r_STG_SCC_beta,r_STG_SCC_lg); 
all_right_VCVS = vertcat(r_ACC_VCVS,r_MOF_VCVS_delta,r_MOF_VCVS_lg,r_SFG_VCVS_delta,r_STG_VCVS_delta,r_STG_VCVS_alpha,r_STG_VCVS_beta,r_STG_VCVS_lg); 
all_right_matrix = vertcat(all_right_SCC,all_right_VCVS); 

DBS_target_SCC = string(repmat('SCC',length(all_right_SCC),1)); 
DBS_target_VCVS = string(repmat('VCVS',length(all_right_VCVS),1)); 

foiroi_label_SCC = vertcat(ACC_label_SCC,MOF_label_SCC_delta,MOF_label_SCC_lg,SFG_label_SCC_delta,STG_label_SCC_delta,STG_label_SCC_alpha,STG_label_SCC_beta,STG_label_SCC_lg); 
foiroi_label_VCVS = vertcat(ACC_label_VCVS,MOF_label_VCVS_delta,MOF_label_VCVS_lg,SFG_label_VCVS_delta,STG_label_VCVS_delta,STG_label_VCVS_alpha,STG_label_VCVS_beta,STG_label_VCVS_lg);  

%% 
tbl_right = array2table(all_right_matrix, 'VariableNames',{'PSD'}); 
tbl_right{:,2} = vertcat(DBS_target_SCC,DBS_target_VCVS);
tbl_right{:,3} = vertcat(foiroi_label_SCC,foiroi_label_VCVS); 
tbl_right{:,4} = string(repmat('R',height(tbl_right),1)); 
tbl_right.Properties.VariableNames = {'PSD','DBS_target','foi_roi','Hemi'}; 

%% Save for python 
addpath(genpath('/Users/anushaallawala/Python_projects/'))
filename = sprintf('/Users/anushaallawala/Python_projects/%s_right_boxplot_data_%s.csv',PatientID,date); 

writetable(tbl_right,filename)

%% 
clear 
PatientID = 'DBSTRD002';
left_tbl = readtable(sprintf('%s_left_boxplot_data_23-Jun-2022.csv',PatientID));
right_tbl = readtable(sprintf('%s_right_boxplot_data_23-Jun-2022.csv',PatientID));
 
all_tbl = vertcat(left_tbl,right_tbl); 
filename_new = sprintf('/Users/anushaallawala/Python_projects/%s_all_boxplot_data_%s.csv',PatientID,date); 

writetable(all_tbl,filename_new)

