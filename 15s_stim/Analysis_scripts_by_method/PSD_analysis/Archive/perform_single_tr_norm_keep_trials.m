%%  function that performs normalization and keeps single trial data 

function [percent_stim_avg_freq,...
    percent_poststim1_avg_freq, percent_poststim2_avg_freq] = perform_single_tr_norm_keep_trials(data,freqs)

%% Describe the data as a sanity check.
disp("The data consists of the power in *physical units* for each bin");
disp("Size of data (trials x channels x freq x time):");
disp(size(data));
num_trials = size(data, 1);
num_channels = size(data, 2);
num_freq = size(data, 3);
num_samples = size(data, 4);

%% 
% Define the different time windows of interest, in samples.
pre_stim_win = 1:4800; %1:5000
stim_win = 5500:20000; %5001:20000;
% TODO: use post_stim windows in analysis.
post_stim_win1 = 20200:25000; 
post_stim_win2 = 25001:30000; 
post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
%% Average across freq band power  

% extract power in freq band. 

data_theta = data(:,:,freqs,:); 
disp(size(data_theta));
data_theta_avg = mean(data_theta,3); 
disp(size(data_theta_avg)); 
data_theta_avg = squeeze(data_theta_avg); 

%* adding line here for heatmap across time 
data_theta_avg_trials = squeeze(mean(data_theta_avg,1)); 

%% calculate perent baseline using single trial normalization 

% 1. Calculate power in baseline window in freq band of interest  
baseline = data_theta(:,:,:,pre_stim_win); 
avg_baseline = mean(baseline,4); % average across TIME (not trials)  
avg_baseline_theta = mean(avg_baseline,3); %average across freq band 

% 2. Convert signal to % change from baseline.  
percent_change = 100*(data_theta - avg_baseline)./avg_baseline; 

% 3. Extract signal for different stimulation windows of interest.  
percent_prestim = percent_change(:,:,:,pre_stim_win); 
percent_stim = percent_change(:,:,:,stim_win); 
percent_post_stim1 = percent_change(:,:,:,post_stim_win1); 
percent_post_stim2 = percent_change(:,:,:,post_stim_win2); 

% normalize using percent of baseline.  
% percent_prestim = 100*(baseline - avg_baseline)./avg_baseline; % pre_stim 
% percent_stim = 100*(stim - avg_baseline)./avg_baseline; % during stim 
% percent_poststim1 = 100*(post_stim1 - avg_baseline)./avg_baseline; % post stim 1
% percent_poststim2 = 100*(post_stim2 - avg_baseline)./avg_baseline; % post stim 2

% 4. Average signal across frequency band of interest. 
percent_prestim_avg_freq = squeeze(mean(percent_prestim,3)); 
percent_stim_avg_freq = squeeze(mean(percent_stim,3)); 
percent_poststim1_avg_freq = squeeze(mean(percent_post_stim1,3)); 
percent_poststim2_avg_freq = squeeze(mean(percent_post_stim2,3)); 