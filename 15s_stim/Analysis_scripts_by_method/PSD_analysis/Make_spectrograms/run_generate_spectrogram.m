%% load info about dataset 
clear all 
PatientID = 'DBSTRD002'; 
contact_config = 'elec25'; 
DBS_target = 'SCC'; 
hemi = 'l'; 
stim_freq = 130; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
figure_type = 'Spectrogram'; 

%% load epoched data 
epocheddatafile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s/15s_stim_all_currdir_singletrial_%s.mat',...
    PatientID,experiment_name,experiment_name); 
epoched_data = load(epocheddatafile); 
fs = epoched_data.metadata.preprocessing.New_SamplingRate; 

switch PatientID 
    case 'DBSTRD001' 
        ch_labels = deblank(epoched_data.metadata.preprocessing.GoodChannelLabels); 
    case 'DBSTRD002' 
        ch_labels = load(sprintf('%s_goodch.mat',experiment_name)); 
        ch_labels = deblank(ch_labels.good_ch_labels); 
    case 'DBSTRD003'
        ch_labels = deblank(epoched_data.metadata.preprocessing.GoodChannelLabels); 
end 

switch DBS_target 
    case 'VCVS'
        [elec1, elec25, elec36, elec47, elec8] = assign_VCVS_conditions(epoched_data,hemi); 
    case 'SCC'
        [elec1, elec234, elec25, elec36, elec47,...
            elec567,elec8] = assign_SCC_conditions(epoched_data,hemi); 
end
%% 
switch contact_config 
    case 'elec25'
        data = elec25; 
    case 'elec1'
        data = elec1; 
    case 'elec8'
        data = elec8; 
    case 'elec36'
        data = elec36; 
    case 'elec47'
        data = elec47; 
    case '234'
        data = elec234; 
    case '567'
        data = elec567; 
end 
        
%% 

num_samples = size(data,4); 
num_ch = size(data,2); 
num_trials = size(data,1); 
t = (0:(num_samples-1))./fs; 
freqs = epoched_data.metadata.decomp_parameters.freqs; 
data = log10(data); 
%% Generate spectrograms 


lower_lim = -2; 
upper_lim = 10; 
freq_range = 1:16; 
outputdir = make_directory(PatientID, figure_type, ...
    experiment_name, contact_config) ; 

generate_spectrogram(num_ch,outputdir,data,ch_labels, t, freqs, ...
    lower_lim, upper_lim, freq_range); 


