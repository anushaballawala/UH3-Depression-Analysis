%%% 1. test to see how long coherence takes for one pair 
    %%% takes 4 days to do one state (poststim, e.g.), 129 x 129 channels
    %%% takes 1.25 minutes to do one state, 3 x 3 channels  
%%% 2. calculate how long it would take for multiple pairs 
%%% 3a. test to see if alternate script from cohen textbook 
%%% 3b. or matlab fxn would take less time 
    %%% mscohere takes 0.16 seconds for 3 x3 channels for one state 

function [] = run_compute_coherence_oscar_stim3(filename)
timeseriesfile = filename;
stim_state = 'stim3'; 
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


%%%%% old code 
% %get subjectID
% subjectID_pattern= 'DBSTRD00\d*'; 
% [start_idx_subid end_idx_subid] = regexp(filename,subjectID_pattern);
% subjectID = extractBetween(filename,start_idx_subid,end_idx_subid);
% PatientID = subjectID{1}; 
% 
% % get DBS_target 
% braintarget = '(l|r)(VCVS|SCC)';
% [start_idx_br, end_idx_br] = regexp(filename,braintarget);
% start_idx_br = start_idx_br(1); end_idx_br = end_idx_br(1); 
% DBS_target = extractBetween(filename,start_idx_br+1,end_idx_br);
% DBS_target = DBS_target{1}; 
% 
% %get hemisphere. 
% hemi = '(l|r)(VCVS|SCC)';
% [start_idx_hemi, end_idx_hemi] = regexp(filename,hemi);
% start_idx_hemi = start_idx_hemi(1); end_idx_hemi = end_idx_hemi(1); 
% 
% switch DBS_target
%     case 'SCC'
%         hemi = extractBetween(filename,start_idx_hemi,end_idx_hemi-3);
%         hemi = hemi{1};
%         
%     case 'VCVS'       
%         hemi = extractBetween(filename,start_idx_hemi,end_idx_hemi-4);
%         hemi = hemi{1};
% end
%% 

disp(PatientID); 
disp(experiment_name); 

[contact_configs,num_contact_configs] = output_contact_config_conditions(PatientID, DBS_target); 

%% get data 

if strcmp(DBS_target,'VCVS') == 1  && strcmp(PatientID, 'DBSTRD001') == 1
[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,ROI_labels,...
    elec1,elec25,elec36,elec47,...
    elec8] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name,filename); 
else
    [timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,ROI_labels,...
    elec1,elec25,elec36,elec47,...
    elec8,elec234,elec567] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name,filename); 
    
end 
%% Define window of interest. 
prestim_param = 'full window'; 
stim_param = 'split into 3'; 
poststim_param = 'full window'; 

[pre_stim_win, stim_win, post_stim_win, stim_win1, stim_win2,...
    stim_win3] = define_window_of_interest(PatientID,...
    prestim_param, stim_param, poststim_param); 

% add to metadata 
metadata.processing.pre_stim_window = pre_stim_win; 
metadata.processing.stim_window1 = stim_win1; 
metadata.processing.stim_window2 = stim_win2; 
metadata.processing.stim_window3 = stim_win3; 
metadata.processing.post_stim_window = post_stim_win;
metadata.processing.stim_window_notes = '5 second duration 1-second post stim picked as post-stim window,... three 1s stim windows chosen'; 

%% 3. Lets try a different way to calculate coherence since chronux takes 4ever
%define coherence params 

win_len = 1500; % 1500 samples or 1.5 s hamming window 
noverlap = 500; % 500 samples 
nfft = win_len; % not currently using this 
num_trials = size(elec1,1); 
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
    case 'stim1' 
        time_win = stim_win1;
    case 'stim2'
        time_win = stim_win2;
    case 'stim3'
        time_win = stim_win3;         
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
metadata.coherence.date = date; 

outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence',PatientID); 

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
    'freqs_beta','freqs_gamma','metadata','all_config_names'); 

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
    end 



