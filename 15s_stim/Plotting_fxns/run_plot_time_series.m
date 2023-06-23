%INPUT: Epoched time series data 
%Output: Plots of time series data for each trial 
clear all 
PatientID = 'DBSTRD002'; 
hemi = 'r' ; 
DBS_target = 'SCC'; 
freq = 130; 
stim_freq = freq; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
% Load single trial time series data 
epocheddatafile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s/15s_stim_all_currdir_timeseries_singletrial_%s.mat',...
    PatientID,experiment_name,experiment_name); 
epoched_data = load(epocheddatafile); 
fs = epoched_data.metadata.preprocessing.New_SamplingRate; 

switch PatientID 
    case 'DBSTRD001' 
        ch_labels = deblank(epoched_data.metadata.preprocessing.GoodChannelLabels); 
    case 'DBSTRD002' 
        ch_labels = load(sprintf('%s_goodch.mat',experiment_name); 
        ch_labels = deblank(ch_labels); 
    case 'DBSTRD003'
        ch_labels = deblank(epoched_data.metadata.preprocessing.GoodChannelLabels); 
end 
%% 
switch DBS_target 
    case 'VCVS'
        [elec1, elec25, elec36, elec47, elec8] = assign_VCVS_conditions(epoched_data,hemi); 
    case 'SCC'
        [elec1, elec234, elec25, elec36, elec47,...
            elec567,elec8] = assign_SCC_conditions(epoched_data,hemi); 
end
%% Plot time-series data 
data = elec1; 
contact_config = 'elec1'; 
plot_time_series_channel(data , PatientID,contact_config, fs, ch_labels, experiment_name) 


%% 

        