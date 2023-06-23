%%% 1. test to see how long coherence takes for one pair 
    %%% takes 4 days to do one state (poststim, e.g.), 129 x 129 channels
    %%% takes 1.25 minutes to do one state, 3 x 3 channels  
%%% 2. calculate how long it would take for multiple pairs 
%%% 3a. test to see if alternate script from cohen textbook 
%%% 3b. or matlab fxn would take less time 
    %%% mscohere takes 0.16 seconds for 3 x3 channels for one state 

%% setup input parameters & run setup file 
%**** generate single trial plot of baseline data ***** 
clear all 
clc 
close all
PatientID = 'DBSTRD001'; DBS_target = 'VCVS'; hemi = 'l'; stim_freq = 130; 
stim_state = 'prestim';  
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
disp(PatientID); 
disp(experiment_name); 
%check 
% check_question = input('is this the right experiment?y/n ','s'); 
% if strcmp(check_question,'n') == 1 
%     error('not the right exp')
% end 

[contact_configs,num_contact_configs] = output_contact_config_conditions(PatientID, DBS_target); 
%stim_state  = input("enter stimstate 'poststim', 'stim' or' prestim'",'s'); 

%% get data 
[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,ROI_labels,...
    elec1,elec25,elec36,elec47,...
    elec8] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name); 


%% Define window of interest. 
prestim_param = 'full window'; 
stim_param = 'full window'; 
poststim_param = 'full window'; 

[pre_stim_win, stim_win, post_stim_win] = define_window_of_interest(PatientID,...
    prestim_param, stim_param, poststim_param); 

% add to metadata 
metadata.processing.pre_stim_window = pre_stim_win; 
metadata.processing.stim_window = stim_win; 
metadata.processing.post_stim_window = post_stim_win;
metadata.processing.stim_window_notes = '5 second duration 1-second post stim picked as post-stim window'; 

%% 3. Lets try a different way to calculate coherence since chronux takes 4ever
%define coherence params 

win_len = 1500; % 1500 samples or 1.5 s hamming window 
noverlap = 500; % 500 samples 
nfft = win_len; % not currently using this 
%num_trials = size(elec1,1); 
num_trials = 1; 
num_channels = size(elec1,2); 
%initialize 
matrix_coh_delta = cell(1,num_contact_configs); matrix_err_delta = cell(1,num_contact_configs); 
matrix_coh_alpha = cell(1,num_contact_configs); matrix_err_alpha = cell(1,num_contact_configs); 
matrix_coh_beta = cell(1,num_contact_configs); matrix_err_beta = cell(1,num_contact_configs); 
matrix_coh_gamma = cell(1,num_contact_configs);  matrix_err_gamma = cell(1,num_contact_configs); 
freqs_delta = {}; freqs_theta = {}; freqs_alpha = {}; freqs_beta = {}; freqs_gamma = {}; 

switch stim_state 
    case 'prestim'
        time_win = pre_stim_win; 
        disp(pre_stim_win)
    case 'stim'
        time_win = stim_win; 
        disp(stim_win)
    case 'poststim'
        time_win = post_stim_win; 
        disp(post_stim_win)
end 


tic 

for i = 1:num_contact_configs
    
   
    contact_config = all_config_names{i}; 
    switch contact_config
        case 'elec25'
            data = permute(elec25,[3,1,2]); 
            data = data(time_win,:,:); 
        case 'elec1'
            data = permute(elec1,[3,1,2]); 
            data = data(time_win,:,:); 
        case 'elec8'
           data = permute(elec8,[3,1,2]); 
            data = data(time_win,:,:); 
        case 'elec36'
            data = permute(elec36,[3,1,2]); 
            data = data(time_win,:,:); 
        case 'elec47'
             data = permute(elec47,[3,1,2]); 
            data = data(time_win,:,:); 
        case 'elec234'
            data = permute(elec234,[3,1,2]); 
            data = data(time_win,:,:); 
        case 'elec567'
            data = permute(elec567,[3,1,2]); 
            data = data(time_win,:,:); 
    end 
    disp(all_config_names{i}); 
    names{i} = all_config_names{i}; 

[matrix_coh_delta{i},matrix_err_delta{i},...
    matrix_coh_theta{i},matrix_err_theta{i},...
    matrix_coh_alpha{i},matrix_err_alpha{i},...
    matrix_coh_beta{i},matrix_err_beta{i},...
    matrix_coh_gamma{i}, matrix_err_gamma{i},...
    freqs_delta{i},freqs_theta{i},freqs_alpha{i},...
    freqs_beta{i},freqs_gamma{i}] = compute_coherence(num_channels,num_trials,data,win_len,noverlap,[],fs); 

end 

toc 

%% save data 

metadata.coherence.win_length = win_len; 
metadata.coherence.noverlap = noverlap; 
metadata.coherence.nfft = win_len; 
metadata.coherence.num_trials = num_trials; 
metadata.coherence.num_channels = num_channels; 
metadata.coherence.window_type = 'Hamming'; 
metadata.coherence.contact_config_order = names;  

outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/%s_firsttrial',PatientID,experiment_name); 
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end 

filename = sprintf('%s/%s_%s_matlabcohfxn.mat',outputdir,stim_state,experiment_name); 


save(filename,'matrix_coh_delta','matrix_err_delta',...
    'matrix_coh_theta','matrix_err_theta',...
    'matrix_coh_alpha','matrix_err_alpha',...
    'matrix_coh_beta','matrix_err_beta',...
    'matrix_coh_gamma', 'matrix_err_gamma',...
    'freqs_delta','freqs_theta','freqs_alpha',...
    'freqs_beta','freqs_gamma','metadata'); 

disp('done')



%% Figures 

caxis_limits = [0 0.5]; 
t_delta = generate_coherence_maps(matrix_coh_delta,names,...
    'delta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
t_theta = generate_coherence_maps(matrix_coh_theta,names,...
    'theta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
t_alpha = generate_coherence_maps(matrix_coh_alpha,names,...
    'alpha',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
t_beta  = generate_coherence_maps(matrix_coh_beta,names,...
    'beta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
t_gamma = generate_coherence_maps(matrix_coh_gamma,names,...
    'gamma',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
disp('done plotting') 
%% 




