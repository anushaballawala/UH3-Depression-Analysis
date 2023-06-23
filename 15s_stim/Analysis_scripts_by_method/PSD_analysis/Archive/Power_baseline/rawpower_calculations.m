%% RawPowercalculations 

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
file = load('15s_stim_elec16_lVCVS_f130.mat');

%% 

fs = file.metadata.preprocessing.New_SamplingRate  ;
% Data matrix.
data = file.indiv2; %selected indiv7 for SCC 
% Trial averaged data matrix. 
data_avg = file.condition2; 

%Define frequency band. 
freqs = 4:7; 
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

%% Average across theta power. 

% extract power in theta band. 

data_theta = data(:,:,freqs,:); 
disp(size(data_theta));
data_theta_avg = mean(data_theta,3); 
disp(size(data_theta_avg)); 
data_theta_avg = squeeze(data_theta_avg); 

%% convert everything to dB using log transform 
data_theta_avg_dB = 10 * log10(data_theta_avg); 

%% plot data across time-domain 
% Plot the data across time for each channel.
time_axis = 0 : (1/fs) : (num_samples/fs); 
time_axis = time_axis(1:end-1);

% parfor ch = 1:num_channels
%     figure('Position',[200, 200, 2200, 600]);
%     hold on;
% 
%     plot(time_axis, squeeze(data_theta_avg_dB(1,ch,:)),'linewidth',2);
%     xline(5,'linewidth', 2.5); 
%     xline(20, 'linewidth', 2.5); 
%     ylabel('Log-transformed power');
%     xlabel('Time[s]');
%     title(sprintf("Channel %s", channel_label(ch)));
%     
%     plot(time_axis, squeeze(data_theta_avg_dB(2,ch,:)),'linewidth',2); 
%     plot(time_axis, squeeze(data_theta_avg_dB(3,ch,:)),'linewidth',2); 
%     plot(time_axis, squeeze(data_theta_avg_dB(4,ch,:)),'linewidth',2); 
%     plot(time_axis, squeeze(data_theta_avg_dB(5,ch,:)),'linewidth',2); 
% 
%     name = sprintf('%s_logtransformthetapowerovertime.png',deblank(channel_label(ch)));
%     saveas(gcf,name); 
% end
%% average data across time 

avg_time_prestim = mean(data_theta_avg_dB(3:5,:,pre_stim_win),3); 
avg_time_stim = mean(data_theta_avg_dB(3:5,:,stim_win),3);
avg_time_poststim1 = mean(data_theta_avg_dB(3:5,:,post_stim_win1),3);
avg_time_poststim2 = mean(data_theta_avg_dB(3:5,:,post_stim_win2),3);
%% 

parfor i = 1:num_channels
    figure() 
    scatter(1:5,avg_time_prestim(:,i),50,'filled')
    hold on 
    scatter(1:5,avg_time_stim(:,i),50,'filled')
    scatter(1:5,avg_time_poststim1(:,i),50,'filled')
    scatter(1:5,avg_time_poststim2(:,i),50,'filled')
    xticks([1 2 3 4 5])
    xticklabels({'Trial 1','Trial2','Trial3','Trial4','Trial5'})
    legend('Pre-stim','stim','Poststim1','Poststim2','Location','best')
    name = sprintf('%s_logtransformthetapower.png',deblank(channel_label(i)));
    saveas(gcf,name); 
end 
%% Calculate averaged data for pre-, during and post-stim 

avg_prestim = mean(avg_time_prestim,1); 
avg_stim = mean(avg_time_stim,1); 
avg_poststim1 = mean(avg_time_poststim1,1);
avg_poststim2 = mean(avg_time_poststim2,1); 

%% Plot 3D data 

good_ch_idx = metadata.preprocessing.GoodChannelsIdx  ; 
[x,y,z] = get_elec_coordinates(good_ch_idx);

%% 

upper_lim = 70; 
lower_lim = 40; 
figure()

h = scatter3(x,y,z,70,'filled','CData',avg_prestim,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('pre-stim')

figure()
h = scatter3(x,y,z,70,'filled','CData',avg_stim,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
caxis('auto')
title('stim')

figure()
h = scatter3(x,y,z,70,'filled','CData',avg_poststim1,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
caxis('auto')
title('post-stim1')

figure()
h = scatter3(x,y,z,70,'filled','CData',avg_poststim2,'MarkerEdgeColor','k'); 
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
caxis('auto')
title('post-stim2')
