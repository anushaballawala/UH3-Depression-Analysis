clear all 
clc 
close all 

%% Load data 

PatientID = 'DBSTRD001'; 
subjecttype = 'DBSTRD' ; 
PatientID = 'DBSTRD001'; baseline_type = 'BaselineFix'; run_num = 7; blk_num = 1; 
experiment_name = sprintf('%s_run-0%d_blk-0%d',baseline_type,run_num,blk_num); 

load(sprintf('%s_%s_Indivtrial_ROI_PSD.mat',PatientID,experiment_name)); 
disp(PatientID); 
disp(experiment_name); 

num_trials = size(ROI_poststim_alpha,1); 
num_ROI = size(ROI_poststim_alpha,2); 
%% Combine all stim states into one 3-D array 

delta = cat(3,ROI_prestim_delta,ROI_stim_delta,ROI_poststim_delta); 
theta = cat(3,ROI_prestim_theta,ROI_stim_theta,ROI_poststim_theta); 
alpha = cat(3,ROI_prestim_alpha,ROI_stim_alpha,ROI_poststim_alpha); 
beta = cat(3,ROI_prestim_beta,ROI_stim_beta,ROI_poststim_beta); 
gamma = cat(3,ROI_prestim_gamma,ROI_stim_gamma,ROI_poststim_gamma); 

%% Generate labels for stim states 

stimstate_labels = zeros(num_trials,num_ROI,3); 

stimstate_labels(:,:,1) = -1; %prestim
stimstate_labels(:,:,2) = 1; %stim
stimstate_labels(:,:,3) = 0; %poststim

%% Save data 

outputdir = sprintf('E:/DBSTRD/%s/Experiments/BaselineFix/Processed Data/PSD',PatientID); 
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end 
filename = sprintf('%s/%s_%s_%s_data_for_stats.mat',...
    outputdir,PatientID,experiment_name,date); 
save(filename,'alpha','beta','theta','delta','gamma','stimstate_labels')





