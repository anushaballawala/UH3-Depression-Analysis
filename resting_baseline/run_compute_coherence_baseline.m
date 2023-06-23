% Functions: Used to generate coherence matrices during baseline recording, using
% chronux function for coherence
% Input: PatientID, experiment info
% Output:* ocherence matrix 
% Dependencies: UH3-Depression-Preprocessing repo, chronux repo 
%
% Anusha Allawala, 09/21
% Updates: 

%% Info  

clear all 
clc 
close all
subjecttype = 'DBSTRD'; 
PatientID = 'DBSTRD001'; baseline_type = 'BaselineFix'; 
run_num = '07'; 
blk_num = '01'; 
experiment_name = sprintf('BaselineFix_run-%s_blk-%s',run_num,blk_num);  
stim_state = 'poststim'; % change *** 

disp(PatientID); 
disp(experiment_name); 
%% Load data 

% Load single trial time series data 
timeseriesfile = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Epoched Data/BaselineFix_run-%s_blk-%s/baseline_all_singletrial_%s.mat',...
   PatientID,run_num,blk_num,experiment_name); 

%timeseries_data = load('15s_stim_all_currdir_singletrial_rVCVS_f130.mat'); 
timeseries_data = load(timeseriesfile); 
fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 

%get metadata & initialize for processing 
metadata = timeseries_data.metadata; 
metadata.processing = []; 


%% make sure structure is the same 

data = timeseries_data.output;
disp(size(data))
switch PatientID 
    case 'DBSTRD001'
        data = data; %*** make sure its the same size 
    case 'DBSTRD002'
        data = permute(data,[3,2,1]); 
end 
disp('new data size is:')
disp(size(data))

%% Data struct info 

% data = timeseries_data.all_ref_signals; 
num_ch = size(data,3); 
sprintf('Num of channels %d',num_ch) 
total_time = size(data,1); 
sprintf('Total time %d',total_time)
time_win = 30000; %time window in samples - 30s to resemble longstim 
num_trials = size(data,2); 
sprintf('Num of trials %d', num_trials)

%% 
fs = 1000; 

win_len = 1500; % 1500 samples or 1.5 s hamming window 
noverlap = 500; % 500 samples 
nfft = win_len; % not currently using this 
num_channels = num_ch; 

pre_stim_win = 1:4900; 
stim_win = 5000:2000;
stim_win1 = 5010:10000;
stim_win2 = 10000:15000;
stim_win3 = 15000:19990;
post_stim_win = 21000:26000; 

prestim_data = data(pre_stim_win,:,:); 
stim1_data = data(stim_win1,:,:); 
stim2_data = data(stim_win2,:,:); 
stim3_data = data(stim_win3,:,:); 
poststim_data = data(post_stim_win,:,:); 

%% 

switch stim_state 
    case 'prestim'
        data_for_coh = prestim_data; 
    case 'stim1' 
        data_for_coh = stim1_data; 
    case 'stim2'
        data_for_coh = stim2_data; 
    case 'stim3'
        data_for_coh = stim3_data; 
    case 'poststim'
        data_for_coh = poststim_data; 
end 


%% 
%initialize 
matrix_coh_delta = cell(1,1); matrix_err_delta = cell(1,1); 
matrix_coh_alpha = cell(1,1); matrix_err_alpha = cell(1,1); 
matrix_coh_beta = cell(1,1); matrix_err_beta = cell(1,1); 
matrix_coh_gamma = cell(1,1);  matrix_err_gamma = cell(1,1); 
freqs_delta = {}; freqs_theta = {}; freqs_alpha = {}; freqs_beta = {}; freqs_gamma = {}; 


[matrix_coh_delta,matrix_err_delta,...
    matrix_coh_theta, matrix_err_theta,...
    matrix_coh_alpha,matrix_err_alpha,...
    matrix_coh_beta,matrix_err_beta,...
    matrix_coh_gamma, matrix_err_gamma,...
    freqs_delta,freqs_theta,freqs_alpha,...
    freqs_beta,freqs_gamma] = compute_coherence(num_channels,num_trials,data_for_coh,win_len,noverlap,[],fs);

disp('computed coherence')
%% save data 

metadata.coherence.win_length = win_len; 
metadata.coherence.noverlap = noverlap; 
metadata.coherence.nfft = win_len; 
metadata.coherence.num_trials = num_trials; 
metadata.coherence.num_channels = num_channels; 
metadata.coherence.window_type = 'Hamming'; 
metadata.coherence.stimstate = stim_state; 


outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/Coherence',PatientID); 
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end 

filename = sprintf('%s/%s_%s_matlabcohfxn.mat',outputdir,experiment_name,stim_state); 


save(filename,'matrix_coh_delta','matrix_err_delta',...
    'matrix_coh_theta','matrix_err_theta',...
    'matrix_coh_alpha','matrix_err_alpha',...
    'matrix_coh_beta','matrix_err_beta',...
    'matrix_coh_gamma', 'matrix_err_gamma',...
    'freqs_delta','freqs_theta','freqs_alpha',...
    'freqs_beta','freqs_gamma','metadata'); 

disp('done')

%% Plot coherence matrices 
trial_type = 'Averaged_across30s'; 
caxis_limits = [0 0.5]; 
t_delta = generate_coherence_maps_baseline(matrix_coh_delta,...
    'delta',PatientID,experiment_name,caxis_limits,trial_type); 
t_theta = generate_coherence_maps_baseline(matrix_coh_theta,...
    'theta',PatientID,experiment_name,caxis_limits,trial_type); 
t_alpha = generate_coherence_maps_baseline(matrix_coh_alpha,...
    'alpha',PatientID,experiment_name,caxis_limits,trial_type); 
t_beta  = generate_coherence_maps_baseline(matrix_coh_beta,...
    'beta',PatientID,experiment_name,caxis_limits,trial_type); 
t_gamma = generate_coherence_maps_baseline(matrix_coh_gamma,...
    'gamma',PatientID,experiment_name,caxis_limits,trial_type); 
disp('done plotting') 










