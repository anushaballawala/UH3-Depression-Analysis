% ***TODO: make bar plots for each channel 
%** repeat ROI barplots for each configuration 
% **figure out why ROI 14 looks weird 


%load mtspectrumc data for all contact configurations 

PatientID = 'DBSTRD002'; 
hemi = 'l'; 
DBS_target = 'SCC'; 
stim_freq = 130; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
maindirectory = dir('E:\DBSTRD\DBSTRD002\Experiments\15s_stim\Processed Data\multitaper_spectrum_all_ch\lSCC_f130\*\*time_series_ROI*.mat'); 
num_files = numel(maindirectory); 
FOI = input('enter FOI ','s'); 

switch FOI 
    case 'delta'
        lower_freq = 1; 
        upper_freq = 4; 
    case 'theta'
        lower_freq = 4; 
        upper_freq = 7; 
    case 'alpha'
        lower_freq = 8; 
        upper_freq = 13; 
    case 'beta'
        lower_freq = 13; 
        upper_freq = 30; 
end 
%% load data 

%%%%%% load precalculated mtspectrum power 
for i = 1: num_files 
    all_data(i) = load(maindirectory(i).name);
    fprintf('Loading.....%s',maindirectory(i).name); 
    %append electrode name at the end of all_data 
end 
num_ch = size(all_data(1).S_all_ch,2); 

%%%%%% Load single trial time series data 
timeseriesfile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
    PatientID,hemi,DBS_target,stim_freq,hemi,DBS_target,stim_freq); 
timeseries_data = load(timeseriesfile); 
fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 


%%  Average across frequency bands. 

pre_stim_freqs = all_data(1).f_all_ch{1,1}; 
stim_freqs = all_data(1).f_stim_all_ch{1, 1} ; 
post_stim_freqs = all_data(1).f_poststim_all_ch{1, 1} ; 

%% find indices for prestim 
%prestim 
freq_prestim_idx = find(pre_stim_freqs >lower_freq & pre_stim_freqs < upper_freq); 
freq_prestim_f = pre_stim_freqs(freq_prestim_idx); 
%stim 
freq_stim_idx = find(stim_freqs >lower_freq & stim_freqs < upper_freq);
freq_stim_f = stim_freqs(freq_stim_idx); 
%poststim
freq_poststim_idx = find(post_stim_freqs >lower_freq & post_stim_freqs < upper_freq);
freq_poststim_f = post_stim_freqs(freq_poststim_idx); 

%% Average across frequency band  

for i = 1:num_ch 
    %prestim 
    prestim_S(i) = mean(all_data(1).S_all_ch{1,i}(freq_prestim_idx),2); 
    %stim
    stim_S(i) = mean(all_data(1).S_stim_all_ch{1,i}(freq_stim_idx),2); 
    %poststim 
    poststim_S(i) = mean(all_data(1).S_poststim_all_ch{1,i}(freq_poststim_idx),2); 
end

%% Obtain SEM by performing multitaper spectrum fxn for each trial 
%Define the different time windows of interest, in samples.

switch PatientID 
    case 'DBSTRD001' 
        pre_stim_win = 1:5000; 
        stim_win = 5001:20000; 
        post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
        post_stim_win2 = 25001:30000; 
        post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
    case 'DBSTRD002'     
        pre_stim_win = 1:4500; 
        stim_win = 5001:20000; 
        post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
        post_stim_win2 = 25001:30000; 
        post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
        %post_stim_total_win = 21000:26000; 
end

% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data = data(:,:,pre_stim_win); 
 indivtr_stim_data = data(:,:,stim_win); 
 indivtr_poststim_total_data = data(:,:,post_stim_total_win); 
 
 %% Perform multitaper spectrum analysis to get indivtr power 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0; %dont average across all trials 
params.tapers = [4 7]; 

for i = 1:num_ch
    for j = 1:num_trials 
%prestim
[S_prestim,f_prestim,~] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(j,i,:),[3,1,2])),params);
f_prestim = f_prestim(2:end); S_prestim = S_prestim(2:end);
mean_pow_prestim = 10*log(S)'; %convert to dB
%stim 
[S_stim,f_stim,~] = mtspectrumc(squeeze(permute(indivtr_stim_data(j,i,:),[3,1,2])),params);
f_stim = f_stim(2:end); S_stim = S_stim(2:end);
mean_pow_stim = 10*log(S_stim)'; %convert to dB
%poststim 
[S_poststim,f_poststim,~] = mtspectrumc(squeeze(permute(indivtr_post_stim_data(j,i,:),[3,1,2])),params);
f_poststim = f_poststim(2:end); S_poststim = S_poststim(2:end);
mean_pow_poststim = 10*log(S_poststim)'; %convert to dB
    end
end

%% group into ROI.

ROI_prestim = generate_ROI_mtspectrum_pow(theta_prestim_S,PatientID); 

%stim 
ROI_stim = generate_ROI_mtspectrum_pow(theta_stim_S,PatientID); 

%poststim 
ROI_poststim = generate_ROI_mtspectrum_pow(theta_poststim_S,PatientID); 

%% subtract to get pre-stim and stim-post 

% prestim - stim 
ROI_preminusstim = ROI_prestim - ROI_stim; 
% stim - poststim 
ROI_stimminuspost = ROI_stim - ROI_poststim; 

%% restructure array to do grouped barplots 

ROI_all = [ROI_preminusstim ROI_stimminuspost];

%% Make bar plots 



figure() 
b = bar(1:17,ROI_all); 
hold on 
[ngroups,nbars] = size(ROI_all); 

%get the xcoordinate of the bars 
x = nan(nbars, ngroups); 
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end 

%plot errorbars 
errorbar(x',ROI_all,ROI_error,'k','linestyle','none')
hold off 
















