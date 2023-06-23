%%%% run_compute_coherence_oscar_poststim_DBStarget.m %%%%

%%% GOAL: compute coherence for all trials of VCVS stim and SCC stim and
%%% generate matrix for each FOI across all channels and combine contact
%%% configurations (VCVS vs SCC) 

%%% OUTPUT: Ch X Ch for each FOI 


%% 
function [] = run_compute_coherence_oscar_poststim_DBStarget(filename)


%filename = '/gpfs/data/dborton/TRD_Project/DBSTRD/DBSTRD002/Experiments/15s_stim/Epoched Data/lVCVS_f130/15s_stim_all_currdir_timeseries_singletrial_lVCVS_f130.mat';
timeseriesfile = filename;


%% from filename, get PatientID, DBS_target, hemisphere 
stim_state = 'poststim'; 
subjectID_pattern= '/DBSTRD00*\w*'; 
[start_idx_subid end_idx_subid] = regexp(timeseriesfile,subjectID_pattern);
PatientID = extractBetween(timeseriesfile,start_idx_subid+1,end_idx_subid);
PatientID = PatientID{1}; 

braintarget = '/(l|r)(VCVS|SCC)_';
[start_idx_br, end_idx_br] = regexp(timeseriesfile,braintarget);
DBS_target = extractBetween(timeseriesfile,start_idx_br+1,end_idx_br-1);
DBS_target = DBS_target{1}; 

hemi = DBS_target(1); 

DBS_target = DBS_target(2:end); 
stim_freq = 130; 

experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq);  

disp(PatientID)
disp(experiment_name)
%% 

[contact_configs,num_contact_configs] = output_contact_config_conditions(PatientID, DBS_target); 

%% get data 

%Specify any conditions where there are <7 current steering configs. 
is_VCVS_DBSTRD001 = strcmp(DBS_target,'VCVS') == 1  && strcmp(PatientID, 'DBSTRD001') == 1; 
is_DBSTRD003 = strcmp(PatientID,'DBSTRD003') == 1; 
    
% Ring contacts not tested with VCVS 001 but with subsequent patients. 
if is_VCVS_DBSTRD001
[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,...
    elec1,elec25,elec36,elec47,...
    elec8] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name,filename); 

elseif is_DBSTRD003
    disp('003')
    %add code here since subset of conditions tested *****
else
    [timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,...
    elec1,elec25,elec36,elec47,...
    elec8,elec234,elec567] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name,filename); 
    
end 
%% Define window of interest. 
prestim_param = 'full window'; 
stim_param = 'split into 3'; 
poststim_param = 'full window'; 

[~, ~, post_stim_win, ~, ~,...
    ~] = define_window_of_interest(PatientID,...
    prestim_param, stim_param, poststim_param); 

% add to metadata 

metadata.processing.post_stim_window = post_stim_win;
metadata.processing.stim_window_notes = '5 second duration 1-second post stim picked as post-stim window,... three 1s stim windows chosen'; 


%% Concatenate all contact configs into one matrix 


if is_VCVS_DBSTRD001 
    %concatenate 5 contacts into matrix and label them 
     all_data = vertcat(elec1,elec25,elec36,elec47,elec8); 
     config_order = {'elec1','elec25','elec36','elec47','elec8'}; 
elseif is_DBSTRD003
    % concatenate 3 contacts into matrix and label them %*** dont know
    % which configs for this yet 
     all_data = vertcat(elec1,elec234,elec25); 

else
    % concatenate 7 contacts into matrix and label them 
    all_data = vertcat(elec1,elec234,elec25,elec36,elec47,elec567,elec8); 
    config_order = {'elec1','elec234','elec25','elec36','elec47','elec567','elec8'}; 
end 

size(all_data)
%% 3. Lets try a different way to calculate coherence since chronux takes 4ever
%define coherence params 

win_len = 1000; % 1000 samples or 1 s hamming window 
noverlap = [];
nfft = win_len; % not currently using this 
num_trials = size(all_data,1); 
num_channels = size(all_data,2);


%initialize 
matrix_coh_delta = []; matrix_err_delta = []; 
matrix_coh_alpha = []; matrix_err_alpha = []; 
matrix_coh_beta = []; matrix_err_beta = []; 
matrix_coh_lowgamma = [];  matrix_err_lowgamma = []; 
matrix_coh_highgamma = [];  matrix_err_highgamma = []; 

time_win = post_stim_win; 
disp(size(post_stim_win))

tic 

data_for_coh = permute(all_data,[3,1,2]); 
data_for_coh = data_for_coh(time_win,:,:); 

[matrix_coh_delta,matrix_err_delta,...
    matrix_coh_theta,matrix_err_theta,...
    matrix_coh_alpha,matrix_err_alpha,...
    matrix_coh_beta,matrix_err_beta,...
    matrix_coh_lowgamma, matrix_err_lowgamma,matrix_coh_highgamma,matrix_err_highgamma,...
    delta_freqs,theta_freqs,alpha_freqs,...
    beta_freqs,lowgamma_freqs,highgamma_freqs,...
    delta_indiv_tr,theta_indiv_tr,alpha_indiv_tr,...
    beta_indiv_tr,lowgamma_indiv_tr,highgamma_indiv_tr] = compute_coherence(num_channels,num_trials,data_for_coh,win_len,noverlap,[],fs); 

toc 

%% save data 

metadata.coherence.win_length = win_len; 
metadata.coherence.noverlap = noverlap; 
metadata.coherence.nfft = win_len; 
metadata.coherence.num_trials = num_trials; 
metadata.coherence.num_channels = num_channels; 
metadata.coherence.window_type = 'Hamming'; 
metadata.coherence.contact_config_order = config_order;  
metadata.coherence.date = date; 

outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence',PatientID); 
%outputdir = sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence',PatientID); 

if ~exist(outputdir,'dir')
    mkdir(outputdir)
end 

filename = sprintf('%s/%s_%s_matlabcohfxn_alltrials.mat',outputdir,stim_state,experiment_name); 

disp('saving....')
save(filename,'data_for_coh','matrix_coh_delta','matrix_err_delta',...
    'matrix_coh_theta','matrix_err_theta',...
    'matrix_coh_alpha','matrix_err_alpha',...
    'matrix_coh_beta','matrix_err_beta',...
    'matrix_coh_lowgamma', 'matrix_err_lowgamma',...
    'matrix_coh_highgamma', 'matrix_err_highgamma',...
    'delta_indiv_tr','theta_indiv_tr','alpha_indiv_tr',...
    'beta_indiv_tr','lowgamma_indiv_tr','highgamma_indiv_tr',...
    'delta_freqs','theta_freqs','alpha_freqs',...
    'beta_freqs','lowgamma_freqs','highgamma_freqs','metadata','all_config_names','config_order'); 

disp('saved')



% %% Figures 
% 
% caxis_limits = [0 0.5]; 
% t_delta = generate_coherence_maps(matrix_coh_delta,names,...
%     'delta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
% t_theta = generate_coherence_maps(matrix_coh_theta,names,...
%     'theta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
% t_alpha = generate_coherence_maps(matrix_coh_alpha,names,...
%     'alpha',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
% t_beta  = generate_coherence_maps(matrix_coh_beta,names,...
%     'beta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
% t_gamma = generate_coherence_maps(matrix_coh_gamma,names,...
%     'gamma',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits); 
% disp('done plotting') 
%% 


