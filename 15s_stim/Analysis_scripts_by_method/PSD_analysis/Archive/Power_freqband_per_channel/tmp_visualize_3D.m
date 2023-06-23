%load trial_averaged files for one Run 
%baseline correct data for each trial type 
%look at average across blocks 
%look at average in individual block 
% might be configured to run with multiple blocks, not just one right now
% ** code edits in progress

clear 
% identify files/blocks for Run# being analyzed.  
files = dir('E:\DBSTRD\DBSTRD001\Experiments\15s_stim\Epoched Data\lSCC_f130\15s_stim_cond3_lSCC_f130.mat'); 

[neuraldatafile, neuraldata_signal] = load_files_from_dir(files); 

num_files = numel(files); 


%% Baseline correction   
data = neuraldata_signal{1, 1}.condition3 ; 
ch_remove = [1,15,31,47,59,73,89,105,121,135];
ch_keep = [2:14,16:30,32:46,48:58,60:72,74:88,90:104,106:120,122:134]; 
data = data(ch_keep,:,:); 
total_time = 25000; 
% stim_win = 0:15000; 
% post_stim = 15100:25000; 
% post_stim_bl = 20000:25000;
% bl_win = 20000:25000; %last five seconds used as baseline 
% stim_off_win = 16000:25000;

pre_stim = 1:5000; 
stim_win = 5001:20000; 
post_stim_win1 = 20001:25000; 
post_stim_win2 = 25001:30000; 
post_stim_total_win = 20001:30000; 

%get average power across baseline window - using pre-stim as baseline  
bl_win_pw = mean(data(:,:,pre_stim),3);

norm_pw = 10*log10( bsxfun(@rdivide,data(:,:,:),bl_win_pw)); %divide power by baseline 
%avg_norm_pw = mean(norm_pw,2);     



%% look at averaged theta power 

%Theta 
theta = 4:7; 

%pre_stim theta 
theta_norm_pw = squeeze(mean(norm_pw(:,theta,post_stim_total_win),2));
theta_norm_pw = squeeze(mean(theta_norm_pw,2));

theta_pre = squeeze(mean(norm_pw(:,theta,pre_stim),2));
theta_pre = squeeze(mean(theta_pre,2)); 

%post_stim theta 
theta_post = squeeze(mean(norm_pw(:,theta,post_stim_total_win),2));
theta_post = squeeze(mean(theta_post,2)); 

theta_post2 = squeeze(mean(norm_pw(:,theta,post_stim_win2),2));
theta_post2 = squeeze(mean(theta_post2,2)); 

%% 
alpha = 8:13; 
alpha_norm_pw = squeeze(mean(norm_pw(:,alpha,post_stim_total_win),2)); 
alpha_norm_pw = mean(alpha_norm_pw,2);

alpha_pre = squeeze(mean(norm_pw(:,alpha,pre_stim),2));
alpha_pre = squeeze(mean(alpha_pre,2)); 

%post_stim theta 
alpha_post = squeeze(mean(norm_pw(:,alpha,post_stim_total_win),2));
alpha_post = squeeze(mean(theta_post,2)); 

alpha_post2 = squeeze(mean(norm_pw(:,alpha,post_stim_win2),2));
alpha_post2 = squeeze(mean(alpha_post2,2));

%% 
highgamma = 24:28;
highgamma_norm_pw = squeeze(mean(norm_pw(:,highgamma,post_stim_total_win),2)); 
highgamma_norm_pw = mean(highgamma_norm_pw,2);

highgamma_pre = squeeze(mean(norm_pw(:,highgamma,pre_stim),2));
highgamma_pre = squeeze(mean(highgamma_pre,2)); 

%post_stim theta 
highgamma_post = squeeze(mean(norm_pw(:,highgamma,post_stim_total_win),2));
highgamma_post = squeeze(mean(theta_post,2)); 

highgamma_post2 = squeeze(mean(norm_pw(:,highgamma,post_stim_win2),2));
highgamma_post2 = squeeze(mean(highgamma_post2,2));


%% create table to export 
% 
% table_as_matrix = [alpha_post,theta_post,highgamma_post]; 
%    tbl = array2table(table_as_matrix, 'VariableNames', ...
%     {'Alpha','Theta','Highgamma'});
% 
% writetable(tbl,'lSCC_f130_elec25.csv'); 





%% create 3-D plot

good_ch_idx = neuraldata_signal{1, 1}.metadata.preprocessing.GoodChannelsIdx  ; 
good_ch_labels = neuraldata_signal{1, 1}.metadata.preprocessing.GoodChannelLabels  ; 
good_ch_labels = good_ch_labels(ch_keep); 
%[x,y,z] = get_elec_coordinates(good_ch_idx); 

[~,~, ChannelInfo]  = xlsread('electrodes.xlsx');  

ChannelInfo = ChannelInfo(2:end,:); 
remove = [59:64,79:128,191:192,207:256]; 
keep = [1:58,65:78,129:190,193:206]; 
ChannelInfo = ChannelInfo(keep,:);


x = cell2mat(ChannelInfo(1:end,2)); 
y = cell2mat(ChannelInfo(1:end,3)); 
z = cell2mat(ChannelInfo(1:end,4));

x2 = x(ch_keep); 
x = x2; 
y = y(ch_keep); 
z = z(ch_keep); 
%%
%plot 

upper_lim = 1;
lower_lim = -1;

figure('Position',[200 200 2200 600])
subplot(1,3,1)
scatter3(x,y,z,'filled','CData',tbl_trialavg.stim)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('pre-stim')

subplot(1,3,2)
scatter3(x,y,z,'filled','CData',theta_post)

colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim')
subplot(1,3,3)
scatter3(x,y,z,'filled','CData',theta_post2)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim-end')
sgtitle('theta')
%% 
upper_lim = 1;
lower_lim = -1;

figure('Position',[200 200 2200 600])
subplot(1,3,1)
scatter3(x,y,z,'filled','CData',highgamma_pre)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('pre-stim')

subplot(1,3,2)
scatter3(x,y,z,'filled','CData',highgamma_post)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim')
subplot(1,3,3)
scatter3(x,y,z,'filled','CData',highgamma_post2)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim-end')
sgtitle('highgamma')

%% 

figure('Position',[200 200 2200 600])
subplot(1,3,1)
scatter3(x,y,z,'filled','CData',alpha_pre)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('pre-stim')

subplot(1,3,2)
scatter3(x,y,z,'filled','CData',alpha_post)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim')
subplot(1,3,3)
scatter3(x,y,z,'filled','CData',alpha_post2)
colormap redblue 
colorbar
caxis([lower_lim upper_lim])
%caxis('auto')
title('post-stim-end')
sgtitle('alpha')

%% 







%******************************************
    srate_new = 1000; fs = 1000;
    dataDuration = 30000;
timeAx = [ 0: (1/fs):dataDuration/srate_new ]; 
timeAx = timeAx(1:end-1);
    freqs = [1 2 3 4 5 6 7 8 9 10 11 12 15 18 21 24 27 30 35 40 45 50 60 70 90 110 130 150]; %center frequencies 



%% 
% good_ch_idx = neuraldata_signal{1, 1}.good_ch  ; 
% good_ch_labels = neuraldata_signal{1, 1}.good_ch_labels;  
lower_lim = -15; 
upper_lim =15; 
% 
%%% TO DO :  save 4. fix loop 
%%% 4. 3-D plots with just baseline averaged 
good_ch_labels = deblank(good_ch_labels); 
 

for ch_idx = 1:125
  figure('Position',[200 200 2200 600])
    contourf(timeAx,freqs(1:15),squeeze(norm_pw(ch_idx,:,:)),'linecolor','none');
        
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        set(0,'defaultTextInterpreter','none'); 
        set(gca, 'YScale', 'log')
        caxis([lower_lim upper_lim])
        colormap redblue
        colorbar('westoutside');
        title(sprintf('%s',good_ch_labels{ch_idx})); 
        saveas(gcf, sprintf('spectrogram_lowfreqs_DBSTRD001_ch%s.png',good_ch_labels{ch_idx})); 
    
end 
    

%% 

h = scatter3(x,y,z,'filled','CData',theta_post);
caxis([-2 2])
Cdata = h.CData; 
cmap = colormap(redblue); 

% cmin = min(Cdata(:)); 
% cmax = max(Cdata(:)); 
cmin = -2;
cmax = 2; 
m = length(cmap); 
idx = fix((Cdata-cmin)/(cmax-cmin)*m)+1; 


RGB = squeeze(ind2rgb(idx,cmap)); 






