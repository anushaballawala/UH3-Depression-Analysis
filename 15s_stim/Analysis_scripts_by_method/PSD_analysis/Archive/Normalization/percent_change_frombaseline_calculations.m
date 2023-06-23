%% Percent change from baseline calculations 

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
%Output:

%% load data 
clear all 

target = 'rSCC'; 

%file = load(sprintf('15s_stim_cond7_%s_f130.mat',target));

freqs = 4:7; 
freqband = 'theta'; 
%% 


% Add if else statement ** 
if strcmp(target,'lVCVS') == 1 
    cond = 16; 
    file = load(sprintf('15s_stim_elec%d_%s_f130.mat',cond,target));
    data = file.indiv2; % data matrix. 
    data_avg = file.condition2; % trial averaged data matrix. 
    
elseif strcmp(target,'rVCVS') == 1
    cond = 16; 
    file = load(sprintf('15s_stim_elec%d_%s_f130.mat',cond,target));   
    data = file.indiv2; 
    data_avg = file.condition2;   
    
elseif strcmp(target,'lSCC') == 1
    cond = 7; 
    file = load(sprintf('15s_stim_cond%d_%s_f130.mat',cond,target));    
    data = file.indiv7; 
    data_avg = file.condition7; 
    
elseif strcmp(target,'rSCC') == 1
    cond = 7; 
    file = load(sprintf('15s_stim_cond%d_%s_f130.mat',cond,target));    
    data = file.indiv7; 
    data_avg = file.condition7; 
    
end 

fs = file.metadata.preprocessing.New_SamplingRate  ;

% Data matrix.
% data = file.indiv7; %selected indiv7 for SCC %indiv2 for VCVS
% % Trial averaged data matrix. 
% data_avg = file.condition7; 


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
pre_stim_win = 1:4800; %1:5000
stim_win = 5500:20000; %5001:20000;
% TODO: use post_stim windows in analysis.
post_stim_win1 = 20200:25000; 
post_stim_win2 = 25001:30000; 
post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);

%% Average across theta power. 

% extract power in theta band. 

data_theta = data(:,:,freqs,:); 
disp(size(data_theta));
data_theta_avg = mean(data_theta,3); 
disp(size(data_theta_avg)); 
data_theta_avg = squeeze(data_theta_avg); 

%* adding line here for heatmap across time 
data_theta_avg_trials = squeeze(mean(data_theta_avg,1)); 


%% METHOD 1 - SINGLE TRIAL BASED NORMALIZATION 
%% calculate perent baseline using single trial normalization 

% caclulate baseline. 
baseline = data_theta(:,:,:,pre_stim_win); 
avg_baseline = mean(baseline,4); % average across time. 
avg_baseline_theta = mean(avg_baseline,3); %average across freq band. 

%normalize the signal 

percent_change = 100*(data_theta - avg_baseline)./avg_baseline; 

% get stim and post stim windows. 
percent_prestim = percent_change(:,:,:,pre_stim_win); 
percent_stim = percent_change(:,:,:,stim_win); 
percent_post_stim1 = percent_change(:,:,:,post_stim_win1); 
percent_post_stim2 = percent_change(:,:,:,post_stim_win2); 

% normalize using percent of baseline.  
% percent_prestim = 100*(baseline - avg_baseline)./avg_baseline; % pre_stim 
% percent_stim = 100*(stim - avg_baseline)./avg_baseline; % during stim 
% percent_poststim1 = 100*(post_stim1 - avg_baseline)./avg_baseline; % post stim 1
% percent_poststim2 = 100*(post_stim2 - avg_baseline)./avg_baseline; % post stim 2

% take average of everything across freqband.  
percent_prestim_avg_freq = squeeze(mean(percent_prestim,3)); 
percent_stim_avg_freq = squeeze(mean(percent_stim,3)); 
percent_poststim1_avg_freq = squeeze(mean(percent_post_stim1,3)); 
percent_poststim2_avg_freq = squeeze(mean(percent_post_stim2,3)); 

% take average across trials 
pc_prestim_avg_trials = squeeze(mean(percent_prestim_avg_freq, 1));
pc_stim_avg_trials = squeeze(mean(percent_stim_avg_freq, 1));
pc_poststim1_avg_trials = squeeze(mean(percent_poststim1_avg_freq, 1));
pc_poststim2_avg_trials = squeeze(mean(percent_poststim2_avg_freq, 1));

%take average across time 
pc_prestim_avg_time = mean(pc_prestim_avg_trials, 2); 
pc_stim_avg_time = mean(pc_stim_avg_trials, 2); 
pc_poststim1_avg_time = mean(pc_poststim1_avg_trials, 2); 
pc_poststim2_avg_time = mean(pc_poststim1_avg_trials, 2); 

%% plot heatmap of all channels 

% take average of trials and theta across all time 
pc_change_trials = squeeze(mean(percent_change,1)); 
pc_change_freqband = squeeze(mean(pc_change_trials,2)); 

%% 
figure('Position',[500 100 900 700])
imagesc(pc_change_freqband)
colormap(redblue)
caxis([-100 200])
colorbar
    xline(5*1000,'linewidth', 2.5); 
    xline(20*1000, 'linewidth', 2.5); 
    
xlabel('Time in samples')
ylabel('Channel number')
set(gca,'YTick',1:4:136,'YTickLabel',channel_label(1:4:136)); 
filename = sprintf('%s_contact8_f130_singletrial_heatmap_theta.png',target); 
saveas(gcf, filename)
%% 
figure('Position',[500 100 900 700])
subplot(311)
scatter(1:136,pc_stim_avg_time,40,'filled')
ylim([-50 200])
title('single trial norm - during stim - percent change from BL, theta')
subplot(312)
scatter(1:136,pc_poststim1_avg_time,40,'filled')
ylim([-50 200])
title('single trial norm - post stim1 - pc change from bl, theta')
subplot(313)
scatter(1:136,pc_poststim2_avg_time,40,'filled')
ylim([-50 200])
title('single trial norm - post stim2 - pc change from bl, theta')
filename = sprintf('%s_contact8_f130_singletrial_pc_change_dist_theta.png',target); 
saveas(gcf, filename)



%% UNCOMMENT EVEYERHTING BELOW 
% % Plot the data across time for each channel.
% time_axis = 0 : (1/fs) : (num_samples/fs); 
% time_axis = time_axis(1:end-1);
% 
% % parfor ch = 1:num_channels
% %     figure('Position',[200, 200, 2200, 600]);
% %     hold on;
% % 
% %     plot(time_axis, squeeze(data_theta_avg_dB(1,ch,:)),'linewidth',2);
% %     xline(5,'linewidth', 2.5); 
% %     xline(20, 'linewidth', 2.5); 
% %     ylabel('Log-transformed power');
% %     xlabel('Time[s]');
% %     title(sprintf("Channel %s", channel_label(ch)));
% %     
% %     plot(time_axis, squeeze(data_theta_avg_dB(2,ch,:)),'linewidth',2); 
% %     plot(time_axis, squeeze(data_theta_avg_dB(3,ch,:)),'linewidth',2); 
% %     plot(time_axis, squeeze(data_theta_avg_dB(4,ch,:)),'linewidth',2); 
% %     plot(time_axis, squeeze(data_theta_avg_dB(5,ch,:)),'linewidth',2); 
% % 
% %     name = sprintf('%s_logtransformthetapowerovertime.png',deblank(channel_label(ch)));
% %     saveas(gcf,name); 
% % end
% %% average data across time 
% 
% 
%% Plot 3D data 

good_ch_idx = metadata.preprocessing.GoodChannelsIdx  ; 
[x,y,z] = get_elec_coordinates(good_ch_idx);

%% 

upper_lim = 50; 
lower_lim = -50; 
figure()

h = scatter3(x,y,z,70,'filled','CData',pc_prestim_avg_time,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('pre-stim')

figure()
h = scatter3(x,y,z,70,'filled','CData',pc_stim_avg_time,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('stim')

figure()
h = scatter3(x,y,z,70,'filled','CData',pc_poststim1_avg_time,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim1')

figure()
h = scatter3(x,y,z,70,'filled','CData',pc_poststim2_avg_time,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim2')
