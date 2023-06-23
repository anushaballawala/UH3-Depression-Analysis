% load data 
%% load data 
clear all 

target = 'lVCVS'; 

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
%% get theta power 

data_freq = data(:,:,freqs,:); 
data_raw = squeeze(mean(data_freq,3)); 
data_log = squeeze(10*log10(mean(data_freq,3))); 

%avg values for each time window 
pre_stim_win = 1:4800; %1:5000
stim_win = 5500:20000; %5001:20000;
post_stim_win1 = 20200:25000; 
post_stim_win2 = 25001:30000; 

prestim_avg = mean(data_raw(:,:,pre_stim_win),3); 
stim_avg = mean(data_raw(:,:,stim_win),3); 
poststim1_avg = mean(data_raw(:,:,post_stim_win1),3); 
poststim2_avg = mean(data_raw(:,:,post_stim_win2),3); 

% average across all channels 
prestim_avg_allchan = mean(prestim_avg(:,131),2); 
stim_avg_allchan = mean(stim_avg(:,131),2); 
poststim1_avg_allchan = mean(poststim1_avg(:,131),2); 
poststim2_avg_allchan = mean(poststim2_avg(:,131),2); 

%%
% get single trial data
figure()
scatter(1:5, prestim_avg_allchan,'filled'); 
hold on 
scatter(1:5, stim_avg_allchan,'filled'); 
scatter(1:5, poststim1_avg_allchan,'filled');
scatter(1:5,poststim2_avg_allchan,'filled'); 
ylabel('Raw theta power ')
xticks([1 2 3 4 5])
xticklabels({'Trial 1', 'Trial 2','Trial 3','Trial 4','Trial 5'}); 
legend('prestim', 'stim', 'poststim1', 'poststim2'); 
%% log transofmred theta power acroaa ll trials 


prestim_avg_log = mean(data_log(:,:,pre_stim_win),3); 
stim_avg_log = mean(data_log(:,:,stim_win),3); 
poststim1_avg_log = mean(data_log(:,:,post_stim_win1),3); 
poststim2_avg_log = mean(data_log(:,:,post_stim_win2),3); 

% average across all channels 
prestim_avg_allchan_log = mean(prestim_avg_log(:,131),2); 
stim_avg_allchan_log = mean(stim_avg_log(:,131),2); 
poststim1_avg_allchan_log = mean(poststim1_avg_log(:,131),2); 
poststim2_avg_allchan_log = mean(poststim2_avg_log(:,131),2); 

%% plot log transformed data 
% get single trial data
figure()
scatter(1:5, prestim_avg_allchan_log,'filled'); 
hold on 
scatter(1:5, stim_avg_allchan_log,'filled'); 
scatter(1:5, poststim1_avg_allchan_log,'filled');
scatter(1:5,poststim2_avg_allchan_log,'filled'); 
ylabel('Log-transformed theta power (dB)')
xticks([1 2 3 4 5])
xticklabels({'Trial 1', 'Trial 2','Trial 3','Trial 4','Trial 5'}); 
legend('prestim', 'stim', 'poststim1', 'poststim2'); 

% log-transform data 


% average across leads (RO)


% generate scatter plot of ROIS  


