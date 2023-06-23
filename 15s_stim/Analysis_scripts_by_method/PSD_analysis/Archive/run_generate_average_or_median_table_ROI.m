% run_generate_average_or_median_table_ROI
clear all 
PatientID = 'DBSTRD001'; 
hemi = 'r' ; 
DBS_target = 'SCC'; 
freq = 130; 
stim_freq = freq; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
% Load single trial time series data 
decompdatafile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s/15s_stim_all_currdir_singletrial_%s.mat',...
    PatientID,experiment_name,experiment_name); 
decomp_data = load(decompdatafile); 
fs = decomp_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(decomp_data.metadata.preprocessing.GoodChannelLabels); 

%%  
contact_config = 'elec25'; 

generate_average_or_median_table_ROI(decompdatafile,contact_config,PatientID,experiment_name,...
    DBS_target); 

%% generate bar chart with error bars for DBS targets in separate script 






