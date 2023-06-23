%% 

% This script relies on the run_setup_15s script,
% generate_heatmap_mtspectrumc_power script 
%% setup input parameters & run setup file 
clear all 
close all
PatientID = 'DBSTRD001'; 
FOI = input('enter FOI ','s'); 
subjecttype = 'DBSTRD' ; 
PatientID = 'DBSTRD001'; baseline_type = 'BaselineFix'; run_num = 7; blk_num = 1; 
experiment_name = sprintf('%s_run-0%d_blk-0%d',baseline_type,run_num,blk_num); 

disp(PatientID); 
disp(experiment_name); 
%% Load data 

% Load single trial time series data 
timeseriesfile = sprintf('E:/%s/%s/Experiments/%s/Preprocessed data/%s/referencedsig_%s.mat',...
    subjecttype, PatientID, baseline_type, experiment_name, experiment_name); 
timeseries_data = load(timeseriesfile); 
fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 

% Load epoched data 
epocheddatafile = sprintf('E:/%s/%s/Experiments/%s/Epoched Data/%s/baseline_all_singletrial_%s.mat',...
    subjecttype, PatientID, baseline_type, experiment_name, experiment_name); 
epoched_data = load(epocheddatafile); 

%get metadata & initialize for processing 
metadata = epoched_data.metadata; 
metadata.processing.PSD = []; 
%% Compute z-score mean and std from full recording 


timeseries_data = timeseries_data.all_ref_signals;
z_mean = @(data)(mean(data,2));
z_std  = @(data)(std(data,0,2)); 


z_mean_baseline = z_mean(timeseries_data); 
z_std_baseline = z_std(timeseries_data); 


%% run through analysis for each configuration 

data = epoched_data.output; 
data = permute(data,[2,3,1]); % changed dims to match what chronux needs
tic

ROI_prestim = []; 
ROI_stim = []; 
ROI_poststim = [] ; 

[ROI_labels,ROI_prestim,ROI_stim,...
            ROI_poststim] = generate_mtspectrumc_power_withstats_baseline(PatientID,...
            metadata,data,z_mean,z_std,lower_freq,...
            upper_freq, ROI_labels); 


toc

%% save data 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/BaselineFix/Processed Data/',PatientID); 
filename = sprintf('%s/%s_%s_Indivtrial_PSD.mat',...
    outputdir,PatientID,experiment_name); 
save(filename,'FOI_prestim_avg_gamma','FOI_stim_avg_gamma','FOI_poststim_gamma',...
    'FOI_prestim_avg_beta','FOI_stim_avg_beta','FOI_poststim_avg_beta',...
    'FOI_prestim_avg_alpha','FOI_stim_avg_alpha','FOI_poststim_avg_alpha',...
    'FOI_prestim_avg_theta','FOI_stim_avg_theta','FOI_poststim_avg_theta',...
    'FOI_prestim_avg_delta','FOI_stim_avg_delta','FOI_poststim_avg_delta',...
    'S_prestim','S_poststim','S_stim','f_prestim','f_stim',...
    'f_poststim','metadata','timeseriesfile')

%% save data for ROI data 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/BaselineFix/Processed Data/',PatientID); 
filename = sprintf('%s/%s_%s_Indivtrial_ROI_PSD.mat',...
    outputdir,PatientID,experiment_name); 
save(filename,'ROI_prestim_delta','ROI_stim_delta','ROI_poststim_delta',...
    'ROI_prestim_theta','ROI_stim_theta','ROI_poststim_theta',...
    'ROI_prestim_alpha','ROI_stim_alpha','ROI_poststim_alpha',...
    'ROI_prestim_beta','ROI_stim_beta','ROI_poststim_beta',...
    'ROI_prestim_gamma','ROI_stim_gamma','ROI_poststim_gamma',...
    'metadata','timeseriesfile')






