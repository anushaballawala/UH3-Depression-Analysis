
clear all 
clc 
PatientID = 'DBSTRD001'; 
load(sprintf('%s_PSD_tstats_25-May-2022.mat',PatientID));
load(sprintf('%s_zscore_bl_poststim_right_SCC_VCVS_052522.mat',PatientID)); 
mult_compare = load(sprintf('%s_mult_compare_data.mat',PatientID));

%% 

num_rSCC_trials = length(labels_rSCC);
num_rVCVS_trials = length(labels_rVCVS); 
total_trials = num_rSCC_trials + num_rVCVS_trials; 

%DPF low gamma  
r_DPF_SCC = diff_data_right(1:num_rSCC_trials,3,5); 
r_DPF_VCVS = diff_data_right(num_rSCC_trials+1:total_trials,3,5);
DPF_label_SCC = string(repmat('DPF_lg_r',num_rSCC_trials,1)); 
DPF_label_VCVS = string(repmat('DPF_lg_r',num_rVCVS_trials,1)); 

%MOF beta 
r_MOF_SCC = diff_data_right(1:num_rSCC_trials,5,4); 
r_MOF_VCVS = diff_data_right(num_rSCC_trials+1:total_trials,5,4);
MOF_label_SCC = string(repmat('MOF_beta_r',num_rSCC_trials,1)); 
MOF_label_VCVS = string(repmat('MOF_beta_r',num_rVCVS_trials,1)); 

all_right_SCC = vertcat(r_DPF_SCC,r_MOF_SCC); 
all_right_VCVS = vertcat(r_DPF_VCVS,r_MOF_VCVS); 
all_right_matrix = vertcat(all_right_SCC,all_right_VCVS); 

DBS_target_SCC = string(repmat('SCC',length(all_right_SCC),1)); 
DBS_target_VCVS = string(repmat('VCVS',length(all_right_VCVS),1)); 

foiroi_label_SCC = vertcat(DPF_label_SCC,MOF_label_SCC); 
foiroi_label_VCVS = vertcat(DPF_label_VCVS,MOF_label_VCVS); 

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
PatientID = 'DBSTRD001';
left_tbl = readtable(sprintf('%s_left_boxplot_data_23-Jun-2022.csv',PatientID));
right_tbl = readtable(sprintf('%s_right_boxplot_data_20-Jun-2022.csv',PatientID));
 
all_tbl = vertcat(left_tbl,right_tbl); 
filename_new = sprintf('/Users/anushaallawala/Python_projects/%s_all_boxplot_data_%s.csv',PatientID,date); 

writetable(all_tbl,filename_new)



