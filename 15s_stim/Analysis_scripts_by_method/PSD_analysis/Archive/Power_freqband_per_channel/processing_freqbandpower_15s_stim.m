%%% The goal of this script is to average and normalize preprocessed data.
%%% The following actions are performed within this script 
%%% 1. normalization using classic method(normalize using pre-stim window,
%%% then convert to log transform and dB units)
%%% 2. average values across time, and across FOI and export values into
%%% table to generate 3-D plots of power across electrodes 
%%% 3. average values across FOI, but not time (to make time-domain plot)
%%% 4. compute standard confidence intervals across trials to plot power
%%% over time 
% Input: preprocessed data, trials x ch x frequencies x time 
% Output: 

clear all 
clc
%function[tbl] = processing_freqbandpower_15s_stim(file,filename,freqs)
file = load('15s_stim_elec16_lVCVS_f130.mat');
% target = 
% Load the data for lVCVS, contact 8.
% Warning: SLOW.  ~1 minute or more to load.
%%
% Sampling rate.
fs = file.metadata.preprocessing.New_SamplingRate  ;
% Data matrix.
data = file.indiv2; %selected indiv7 for SCC 
% Trial averaged data matrix. 
data_avg = file.condition2; 

%Define frequency band. 
freqs = 8:13; 
% Get other channel info from metadata; 
metadata = file.metadata; 
channel_label = metadata.preprocessing.GoodChannelLabels; 
ch_idx = metadata.preprocessing.GoodChannelsIdx; 
%%
% Describe the data as a sanity check.
disp("The data consists of the power in *physical units* for each bin");
disp("Size of data (trials x channels x freq x time):");
disp(size(data));
num_trials = size(data, 1);
num_channels = size(data, 2);
num_freq = size(data, 3);
num_samples = size(data, 4);
%% 
% Define the different time windows of interest, in samples.
pre_stim_win = 1:5000; 
stim_win = 5001:20000;
% TODO: use post_stim windows in analysis.
post_stim_win1 = 20001:25000; 
post_stim_win2 = 25001:30000; 
post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);

%% 
% Normalize data by the average power in the pre-stim window.
% (1) Compute average baseline power across trials and time window
% This results in 1 average spectral power datapt for each freq 
bl_pw = mean(mean(data(:,:,:,pre_stim_win),4),1);

% (2) Divide power by baseline to get normalized power (in dB)
bl_normalize_db = @(x) 10 * log10(bsxfun(@rdivide, x, bl_pw));

% (3) Actually perform the normalization
for trial = 1:num_trials
    norm_pw_dB(trial, :, :, :) = bl_normalize_db(data(trial, :, :, :));
end 
%%
% index freqband power for each trial
% Define the range of frequency indexes for the freqband band.
% By design, these indexes match the actual frequencies (in Hz
average_over_freqband = @(x) mean(x(:, :, freqs, :), 3);

average_freqband_norm_pw_dB = average_over_freqband(norm_pw_dB);

% TLDR; Input x to average_over_time is assumed to have exactly 4 different 
% dimensions
% Matlab will remove all trailing dimensions of size 1.  
% But it does not remove internal dimensions.  This can lead to bugs.  
% Be careful.

average_over_time = @(x, win) mean(x(:, :, :, win), 4);

% pre stim 
average_pre_stim_freqband_norm_pw_dB = ...
    average_over_time(average_freqband_norm_pw_dB, pre_stim_win);

% during stim 
average_stim_freqband_norm_pw_dB = ...
    average_over_time(average_freqband_norm_pw_dB, stim_win);

% first post-stim window  
average_post_stim1_freqband_norm_pw_dB = ...
    average_over_time(average_freqband_norm_pw_dB, post_stim_win1);

% second post-stim window 
average_post_stim2_freqband_norm_pw_dB = ...
    average_over_time(average_freqband_norm_pw_dB, post_stim_win2);


%% Compute means and confidence intervals etc. across trials.

Initialize.
mu = zeros(num_channels, num_samples);
sd = mu;
sdmn = mu;
CI_upper = mu;
CI_lower = mu;

parfor ch = 1:num_channels
    mu(ch, :) = squeeze(mean(average_freqband_norm_pw_dB(:, ch, :, :), 1)); 
    sd(ch, :) = squeeze(std(average_freqband_norm_pw_dB(:, ch, :, :), 0, 1)); 
    sdmn(ch, :) = sd(ch, :) ./ sqrt(num_trials); 
    CI_upper(ch, :) = mu(ch, :) + 2 * sdmn(ch, :); 
    CI_lower(ch, :) = mu(ch, :) - 2 * sdmn(ch, :); 
end 
%%
% compute confidence intervals using fxn 

%% Compute average power across all trials for each time window.   

average_over_trials = @(x) mean(x,1); 

% average across all trials for pre-stim window. 
average_pre_stim_freqband_alltrials = ...
    average_over_trials(average_pre_stim_freqband_norm_pw_dB); 

% average across all trials during stimulation window. 
average_stim_freqband_alltrials = ...
    average_over_trials(average_stim_freqband_norm_pw_dB); 

% average across all trials for first post-stimulation window 
average_post_stim1_freqband_alltrials = ...
    average_over_trials(average_post_stim1_freqband_norm_pw_dB); 

%average across all trials for second post-stimulation window 
average_post_stim2_freqband_alltrials = ...
    average_over_trials(average_post_stim2_freqband_norm_pw_dB); 

%% Compile all values into table. 

[x,y,z] = get_elec_coordinates(ch_idx); 
table = []; 
table = [x,y,z,average_pre_stim_freqband_alltrials',...
    average_pre_stim_freqband_alltrials',...
    average_post_stim1_freqband_alltrials',...
    average_post_stim2_freqband_alltrials']; 


% convert array to table 
tbl = array2table(table,'VariableNames',{'x','y','z','pre-stim','stim',...
    'post_stim1','post_stim2'}); 
%end 
filename_csv_trial = 'E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/050621/lVCVS_f130_elec8_alpha_singletrialblnorm.csv';

writetable(tbl,filename_csv_trial);
%% Calculate trial_averaged baseline data. 

% normalize trialaveraged data to averaged baseline 
bl_pw_trialavg = mean(data_avg(:,:,pre_stim_win),3);

% (2) Divide power by baseline to get normalized power (in dB)
bl_normalize_db_trialavg = @(x) 10 * log10(bsxfun(@rdivide, x, bl_pw_trialavg));

% (3) Actually perform the normalization

norm_pw_dB_trialavg(:, :, :) = bl_normalize_db_trialavg(data_avg(:, :, :));

%% %% Compute average power across time for each time window.   

% average across frequency band of interest 
norm_pw_dB_trialavg = norm_pw_dB_trialavg(:,freqs,:); 

norm_pw_dB_trialavg_freq = squeeze(mean(norm_pw_dB_trialavg,2)); 

% average across all trials during pre-stimulation window. 
avg_pre_stim_freqband_trialavg = mean(norm_pw_dB_trialavg_freq(:,pre_stim_win),2); 

% average across all trials during stimulation window. 
avg_stim_freqband_trialavg = mean(norm_pw_dB_trialavg_freq(:,stim_win),2); 

% average across all trials for first post-stimulation window 
avg_post_stim1_freqband_trialavg = mean(norm_pw_dB_trialavg_freq(:,post_stim_win1),2); 

%average across all trials for second post-stimulation window 
avg_post_stim2_freqband_trialavg = mean(norm_pw_dB_trialavg_freq(:,post_stim_win2),2); 


%% Generate table with trial_averaged baseline data. 

table_trialavg = [x,y,z,avg_pre_stim_freqband_trialavg,...
    avg_stim_freqband_trialavg,...
    avg_post_stim1_freqband_trialavg,...
    avg_post_stim2_freqband_trialavg]; 

% convert array to table 
tbl_trialavg = array2table(table_trialavg,'VariableNames',{'x','y','z','pre-stim','stim',...
    'post_stim1','post_stim2'}); 
%end 

filename_csv = 'E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/050621/lVCVS_f130_elec8_alpha.csv';

writetable(tbl_trialavg,filename_csv);

