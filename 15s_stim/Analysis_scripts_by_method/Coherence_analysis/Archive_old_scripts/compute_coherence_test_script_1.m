%%% 1. test to see how long coherence takes for one pair 
    %%% takes 4 days to do one state (poststim, e.g.), 129 x 129 channels
    %%% takes 1.25 minutes to do one state, 3 x 3 channels  
%%% 2. calculate how long it would take for multiple pairs 
%%% 3a. test to see if alternate script from cohen textbook 
%%% 3b. or matlab fxn would take less time 
    %%% mscohere takes 0.16 seconds for 3 x3 channels for one state 


%% setup input parameters & run setup file 
clear all 
close all
PatientID = 'DBSTRD001'; DBS_target = 'SCC'; hemi = 'r'; stim_freq = 130; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
disp(experiment_name); 
%check 
check_question = input('is this the right experiment?y/n','s'); 
if strcmp(check_question,'n') == 1 
    error('not the right exp')
end 

contact_configs = {'elec1','elec25','elec234','elec36','elec47','elec234','elec8'}; 
num_contact_configs = numel(contact_configs); 

%% get data 
[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,ROI_labels,...
    elec1,elec25,elec36,elec47,...
    elec8,elec234,elec567] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name); 


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


% %% compute coherence for a handful of channels 
% 
% % Start with coherence between one channel and the rest 
% window = 2.5; %2.5 second window 
% winstep = 0.001; % 10 ms steps because sampling rate is 1000 Hz 
% movingwin = [window winstep]; 
% params.tapers = [4 6]; 
% params.pad = 0; %padding 
% params.Fs = 1000; %1000 Hz 
% params.fpass = [4 7]; %theta band 
% params.err = [1 0.05]; % 95 % confidence interva; 
% params.trialave = 1; 
% 
% prestim = elec1(:,:,pre_stim_win); 
% prestim_data = permute(prestim,[3,1,2]); 
% 
% poststim_data = elec1(:,:,pre_stim_win);
% poststim_data = permute(poststim_data,[3,1,2]);
% 
% %num_channels = size(prestim_data,3);
% num_channels = 3; 
% matrix_prestim = zeros(num_channels, num_channels);
% matrix_poststim = zeros(num_channels, num_channels);
% %% 2. Lets try coherence with chronux 
% tic 
% num_channels = 5; 
% for i = 1:num_channels
%     for j = i:num_channels
%         %%% for prestim 
%         first_ch = prestim_data(:,:,i);
%         second_ch = prestim_data(:,:,j);
%         C=cohgramc(first_ch,...
%            second_ch,movingwin,params);
%        % take average across coh values. 
%         val = mean(mean(C)); 
%         % generate matrix of coherence values. 
%         matrix_prestim(i,j) = val;
%         matrix_prestim(j,i) = val;
%         
%         %%%% for poststim 
%         first_ch_post = poststim_data(:,:,i);
%         second_ch_post = poststim_data(:,:,j);
%         [C_post,phi_post,S12_post,S1_post,S2_post,t_post,...
%             f_post,confC_post,phistd_post]=cohgramc(first_ch_post,...
%            second_ch_post,movingwin,params);
%        % take average across coh values. 
%         val_post = mean(mean(C_post)); 
%         % generate matrix of coherence values. 
%         matrix_poststim(i,j) = val_post;
%         matrix_poststim(j,i) = val_post;
%     end
% end
% toc 
%% 3. Lets try a different way to calculate coherence since chronux takes 4ever
%define coherence params 

win_length = 1500; % 1500 samples or 1.5 s hamming window 
noverlap = 500; % 500 samples 
nfft = win_length; % not currently using this 
num_trials = size(elec1,1); 
num_channels = size(elec1,2); 
%initialize 
matrix_coh_delta = cell(1,num_contact_configs); matrix_err_delta = cell(1,num_contact_configs); 
matrix_coh_alpha = cell(1,num_contact_configs); matrix_err_alpha = cell(1,num_contact_configs); 
matrix_coh_beta = cell(1,num_contact_configs); matrix_err_beta = cell(1,num_contact_configs); 
matrix_coh_gamma = cell(1,num_contact_configs);  matrix_err_gamma = cell(1,num_contact_configs); 
freqs_delta = {}; freqs_theta = {}; freqs_alpha = {}; freqs_beta = {}; freqs_gamma = {}; 


stim_state  = input("enter stimstate 'poststim', 'stim' or' prestim'",'s'); 

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
% for prestim 
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
    freqs_beta{i},freqs_gamma{i}] = compute_coherence(num_channels,num_trials,data,window,noverlap,[],fs); 

end 

toc 

%% save data 

metadata.coherence.win_length = win_length; 
metadata.coherence.noverlap = noverlap; 
metadata.coherence.nfft = win_length; 
metadata.coherence.num_trials = num_trials; 
metadata.coherence.num_channels = num_channels; 
metadata.coherence.window_type = 'Hamming'; 
metadata.coherence.contact_config_order = names;  

outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/%s',PatientID,experiment_name); 
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

t_delta = generate_coh_maps(matrix_coh_delta,names,'delta',PatientID,experiment_name,ch_labels,stim_state); 
t_theta = generate_coh_maps(matrix_coh_theta,names,'theta',PatientID,experiment_name,ch_labels,stim_state); 
t_alpha = generate_coh_maps(matrix_coh_alpha,names,'alpha',PatientID,experiment_name,ch_labels,stim_state); 
t_beta  = generate_coh_maps(matrix_coh_beta,names,'beta',PatientID,experiment_name,ch_labels,stim_state); 
t_gamma = generate_coh_maps(matrix_coh_gamma,names,'gamma',PatientID,experiment_name,ch_labels,stim_state); 

t_delta_cb = generate_coh_maps_cb(matrix_coh_delta,names,'delta',PatientID,experiment_name,ch_labels,stim_state); 
t_theta_cb = generate_coh_maps_cb(matrix_coh_theta,names,'theta',PatientID,experiment_name,ch_labels,stim_state); 
t_alpha_cb = generate_coh_maps_cb(matrix_coh_alpha,names,'alpha',PatientID,experiment_name,ch_labels,stim_state); 
t_beta_cb  = generate_coh_maps_cb(matrix_coh_beta,names,'beta',PatientID,experiment_name,ch_labels,stim_state); 
t_gamma_cb = generate_coh_maps_cb(matrix_coh_gamma,names,'gamma',PatientID,experiment_name,ch_labels,stim_state); 



% function 
function t = generate_coh_maps(data_for_fig,names,FOI,PatientID,experiment_name,ch_labels,stim_state)
f = figure; 
f.Position = [100 100 1400 800]; 
t = tiledlayout(2,4); 
for i = 1:7 
nexttile 
imagesc(data_for_fig{i}) 
title(sprintf('%s %s',names{i},FOI))
colormap redbluecmap ; 
colorbar 
t.TileSpacing = 'compact'; 
t.Padding = 'compact'; 

filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081621_Coherence/%s_%s_%s',...
    PatientID,stim_state,experiment_name,FOI); 
filename_png = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081621_Coherence/%s_%s_%s.png',...
    PatientID,stim_state,experiment_name,FOI); 
saveas(gcf,filename); 
saveas(gcf,filename_png); 
end 
end 


% function 
function t = generate_coh_maps_cb(data_for_fig,names,FOI,PatientID,experiment_name,ch_labels,stim_state)
f = figure; 
f.Position = [100 100 1400 800]; 
t = tiledlayout(2,4); 
for i = 1:7 
nexttile 
imagesc(data_for_fig{i}) 
title(sprintf('%s %s',names{i},FOI))
colormap redbluecmap ; 
colorbar 
caxis([0 0.3]) 
t.TileSpacing = 'compact'; 
t.Padding = 'compact'; 

filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081621_Coherence/%s_%s_%s_colorbar.fig',...
    PatientID,stim_state,experiment_name,FOI); 
saveas(gcf,filename); 

end 
end 


%% 





