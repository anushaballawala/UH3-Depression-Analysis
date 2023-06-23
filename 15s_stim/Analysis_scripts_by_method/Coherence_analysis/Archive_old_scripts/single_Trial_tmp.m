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
PatientID = 'DBSTRD001'; baseline_type = 'BaselineFix'; run_num = 7; blk_num = 1; 
experiment_name = sprintf('%s_run-0%d_blk-0%d',baseline_type,run_num,blk_num); 

disp(PatientID); 
disp(experiment_name); 
%% Load artificially epoched data 
epocheddatafile = sprintf('E:/%s/%s/Experiments/%s/Epoched Data/%s/baseline_all_singletrial_%s.mat',...
    subjecttype, PatientID, baseline_type, experiment_name, experiment_name); 
epoched_data = load(epocheddatafile); 


%% Set up metadata and fake time windows to match stim 

fake_prestim_win = 1:5000;
fake_stim_win = 5010:20000; 
%fake_stim_win = 5010:10010; 
%fake_stim_win = 10010:15010; 
%fake_stim_win = 15010:20000; 
fake_poststim_win = 21000:26000;

%get metadata & initialize for processing 
metadata = epoched_data.metadata; 
metadata.processing = []; 
fs = 1000; 
num_ch = size(epoched_data.output,3); 
total_time = size(epoched_data.output,1); 
time_win = total_time; %time window in samples - 30s to resemble longstim 
num_trials = 9; 

%setup fake data 
fake_prestim_data = epoched_data.output(fake_prestim_win,:,:); 
fake_stim_data = epoched_data.output(fake_stim_win,:,:); 
fake_poststim_data = epoched_data.output(fake_poststim_win,:,:); 
%% Run coherence and get single trial output 

win_len = 1500; % 1500 samples or 1.5 s hamming window 
noverlap = 500; % 500 samples 
nfft = win_len; % not currently using this 
num_channels = num_ch;  




%% 

matrix_coh_theta = {}; theta_freqs = {}; cell_coh_theta = {}; 
matrix_coh_delta = {}; delta_freqs = {}; cell_coh_delta = {}; 
matrix_coh_beta = {}; beta_freqs = {}; cell_coh_beta = {}; 
matrix_coh_alpha = {}; alpha_freqs = {}; cell_coh_alpha = {}; 
matrix_coh_gamma = {}; gamma_freqs = {}; cell_coh_gamma = {}; 

for i = 1:num_trials 
data = fake_poststim_data(:,i,:); 
[matrix_coh_delta,...
    matrix_coh_theta,...
    matrix_coh_alpha,...
    matrix_coh_beta,...
    matrix_coh_gamma,...
    delta_freqs,theta_freqs,alpha_freqs,...
    beta_freqs,gamma_freqs] = test_fxn_theta(num_ch,data,win_len,noverlap,[],fs); 

cell_coh_theta{i} = matrix_coh_theta; 
cell_coh_delta{i} = matrix_coh_delta; 
cell_coh_alpha{i} = matrix_coh_alpha; 
cell_coh_beta{i} = matrix_coh_beta; 
cell_coh_gamma{i} = matrix_coh_gamma; 

end 
%% Average across all trials 

avg_theta = mean(cat(4,cell_coh_theta{:}),4); 
avg_delta = mean(cat(4,cell_coh_delta{:}),4);
avg_alpha = mean(cat(4,cell_coh_alpha{:}),4);
avg_beta = mean(cat(4,cell_coh_beta{:}),4);
avg_gamma = mean(cat(4,cell_coh_gamma {:}),4);


%% save data 

data_type = 'avgtrials'; 
fake_data_type = 'fake_poststim'; 
metadata.coherence.win_length = win_len; 
metadata.coherence.noverlap = noverlap; 
metadata.coherence.nfft = win_len; 
metadata.coherence.num_trials = num_trials; 
metadata.coherence.num_channels = num_channels; 
metadata.coherence.window_type = 'Hamming'; 

outputdir = sprintf('E:/DBSTRD/%s/Experiments/BaselineFix/Processed Data/Coherence/%s',PatientID,experiment_name); 
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end 

filename = sprintf('%s/%s_%s_%s_matlabcohfxn.mat',outputdir,experiment_name,...
    data_type,fake_data_type); 


save(filename,'avg_delta',...
    'avg_theta',...
    'avg_alpha',...
    'avg_beta',...
    'avg_gamma',...
    'delta_freqs','theta_freqs','alpha_freqs',...
    'beta_freqs','gamma_freqs','metadata'); 

disp('done')
%% Plot avg coherence matrices 
caxis_limits = [0 0.5]; 
trial_type = 'fakepoststim_avg'; 
FOI = 'delta'; 
generate_coherence_maps_baseline(avg_delta,FOI,PatientID,...
    experiment_name,caxis_limits,trial_type)



%% Plot coherence matrices for each trial 

trial_type = fake_data_type; 
caxis_limits = [0 0.5]; 

t_delta = generate_coherence_maps_baseline_singletr(cell_coh_delta,...
    'delta',PatientID,experiment_name,caxis_limits,trial_type,num_trials); 

t_theta = generate_coherence_maps_baseline_singletr(cell_coh_theta,...
    'theta',PatientID,experiment_name,caxis_limits,trial_type,num_trials); 

t_alpha = generate_coherence_maps_baseline_singletr(cell_coh_alpha,...
    'alpha',PatientID,experiment_name,caxis_limits,trial_type,num_trials); 

t_beta  = generate_coherence_maps_baseline_singletr(cell_coh_beta,...
    'beta',PatientID,experiment_name,caxis_limits,trial_type,num_trials); 

t_gamma = generate_coherence_maps_baseline_singletr(cell_coh_gamma,...
    'gamma',PatientID,experiment_name,caxis_limits,trial_type,num_trials); 
disp('done plotting') 










