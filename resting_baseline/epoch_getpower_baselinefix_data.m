%% Sricpt for segmenting baseline data based on lags from autocorrelation analysis. 
% Output: power for each segmented window from baseline data 


% load bl data 
clear all 
clc 
addpath(genpath('/Users/anushaallawala/Data/')); 
PatientID = 'DBSTRD003'; 
%bl_data = load('/Users/anushaallawala/Data/Baseline_Data/003/referencedsig_BaselineFix_run-01_blk-01.mat'); 
bl_data = load('/Users/anushaallawala/Data/Baseline_Data/003/referencedsig_BaselineFix_date-04-18-2021.mat');
data = bl_data.all_ref_signals ; 
experiment_name = 'BaselineFix_date-04-18-2021'; %** extract from filename 
%% 
% segment based on 75th percentile # of lags 

switch PatientID 
    case 'DBSTRD001'
        lags = 9800; 
    case 'DBSTRD002' 
        lags = 7900; 
    case 'DBSTRD003'
        lags = 6400; 
end 

%% Calculating # of trials that baseline rec will give after dealing with temporal correlation 
% 
fs = 1000; 
lag_secs = lags/fs; 

trial_time_sec = lag_secs+5; 
trial_time_samples = lags+5000; 
total_samples = length(data); 
total_time_sec = total_samples/fs; 
total_time_min = total_time_sec/60; 

num_trials = round(total_time_sec/trial_time_sec); 

disp(num_trials)

%metadata

metadata.epoching.trial_time = trial_time_samples;
metadata.epoching.total_samples = total_samples;
%% Segment data 

    trial_length_sec = 5;
    lag_sec = lag_secs;
    metadata.epoching.lag_sec = lag_sec;
    full_trial_length_sec = trial_length_sec + lag_sec;
    
    sec2samples = @(x) round(x * fs);
    trial_length_samples = sec2samples(trial_length_sec);
    full_trial_length_samples = sec2samples(full_trial_length_sec);
     metadata.epoching.full_trial_length_samples = full_trial_length_samples; 


    output = []; 
    num_channels = size(data, 1);
    for start_idx = 1:full_trial_length_samples:total_samples
        end_idx = start_idx + trial_length_samples - 1; 
        trial_data = data(:, start_idx:end_idx);
        trial_data = reshape(trial_data, 1, num_channels, trial_length_samples);        
        output = vertcat(output, trial_data);
    end
    
    disp(size(output)) % trials X channels X samples 

    indivtr_data = output; 
%  
%% Compute PSD using multitaper for each contact configuration 
 
ch_labels = deblank(bl_data.metadata.preprocessing.GoodChannelLabels); 
metadata = bl_data.metadata; 
metadata.processed.date = date; 
addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/Toolboxes/chronux_2_12.v03')); 
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/chronux_2_12'));
disp(size(data))
params.Fs = fs;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 1;
params.tapers = [4 7]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
S = []; 
for i=1:num_channels
    for j = 1:num_trials 
        [S(j,i,:),f,Serr(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_data(j,i,:),[3,1,2])),params);
    end    
end 

%% save PSD

outputdir = sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/%s/%s',PatientID,experiment_name, 'PSD'); 
if ~exist(outputdir,'dir'); mkdir(outputdir); end   
    filename = sprintf('%s/%s_mtspectrum_pow_indivtr.mat',outputdir,experiment_name); 

    save(filename, 'f', 'S', 'Serr','metadata') 
    disp('PSD for all channels saved');
    
%-------------- END OF CODE ---------------% 

%% Save epoched time series data 

metadata.preprocessing.num_lags = lags; 

outputdir = '/Volumes/TRD_Project/DBSTRD/DBSTRD003/Experiments/BaselineFix/Epoched Data/'; 
if ~exist(outputdir,'dir'); mkdir(outputdir); end   
    filename = sprintf('%s/baseline_all_singletrial_%s.mat',outputdir,experiment_name); 

    save(filename, 'indivtr_data','metadata') 
    disp('Epoched time series for all channels saved');












