%% METHOD 2 - TRIAL AVERAGED BASELINE NORMALIZATION 
%% calculate percent baseline using averaged data 

%100* (signal-baseline/baseline) *** continue code here 
mean_trials = squeeze(mean(data_theta_avg,1)); 
baseline = mean_trials(:,pre_stim_win); 
avg_baseline = mean(baseline,2); 
stim = mean_trials(:,stim_win); 
post_stim1 = mean_trials(:,post_stim_win1); 
post_stim2 = mean_trials(:,post_stim_win2); 

% percent change for pre stim window. 
percent_prestim = 100*((baseline - avg_baseline)./avg_baseline); 
avg_percent_prestim = mean(percent_prestim,2); 

% percent change for stim.
percent_stim = 100*((stim - avg_baseline)./avg_baseline); 
avg_percent_stim = mean(percent_stim,2); 

% percent change for post stim window 1
percent_poststim1 = 100*((post_stim1 - avg_baseline)./avg_baseline);
avg_percent_poststim1 = mean(percent_poststim1,2); 
% percent change for post stim window 2 
percent_poststim2 = 100*((post_stim2 - avg_baseline)./avg_baseline);
avg_percent_poststim2 = mean(percent_poststim2,2); 


%% calculate percent baseline across time 

normalized_data_time = mean_trials./avg_baseline; 

%% 
figure('Position',[500 100 900 700])
imagesc(normalized_data_time)
colormap(redblue)
caxis([-10 20])
colorbar
    xline(5*1000,'linewidth', 2.5); 
    xline(20*1000, 'linewidth', 2.5); 
    
xlabel('Time in samples')
ylabel('Channel number')
set(gca,'YTick',1:4:136,'YTickLabel',channel_label(1:4:136)); 
filename = sprintf('%s_contact8_f130_avgtrial_heatmap_theta.png',target); 
saveas(gcf, filename)
%% 
figure('Position',[500 100 900 700])
subplot(311)
scatter(1:136,avg_percent_stim,40,'filled')
ylim([-50 200])
title('trial averaged norm - during stim - percent change from BL, theta')
subplot(312)
scatter(1:136,avg_percent_poststim1,40,'filled')
ylim([-50 200])
title('trial averaged norm - post stim1 - pc change from bl, theta')
subplot(313)
scatter(1:136,avg_percent_poststim2,40,'filled')
ylim([-50 200])
title('trial averaged norm - post stim2 - pc change from bl, theta') 
filename = sprintf('%s_contact8_f130_avgtrial_pc_change_dist_theta.png',target); 
saveas(gcf, filename)

