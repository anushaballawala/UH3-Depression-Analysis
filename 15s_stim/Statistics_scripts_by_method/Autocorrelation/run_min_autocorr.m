%% FUNCTION: calculating temporal autocorrelation across baseline spectral data
%  for use against 15-s stim data in TRD subjects. The minimum # of lags 
% it takes for autocorrelation to drop to zero or close to zero is what is 
% used to generate artificial windows for analysis in baseline data. The # of
% lags is skipped when generating these windows (.ppt for reference is on
% Desktop). 

% INPUTS: baseline data struct, with 1xch struct, containing freqs X time
% power 
% OUTPUTS: figures using autocorr fxn in MATLAB, minimum # of lags and
% distribution of lag# for each ROI and FOI(this data is not saved, only figs are) 
% NOTE: Code can be buggy and not tested in a while, this was created just to be able to get # of lags
% info 

% Written by AA 11/2021 (with input from MH). 

%% Load baseline data 
clc 
clear 
close all 

addpath(genpath('/Users/anushaallawala/Data/Baseline_Data'))

PatientID = 'DBSTRD003';

% load data after performing morelet wavelet analysis
%bl_data = load('/Users/anushaallawala/Data/Baseline_Data/002/decomp_signal_BaselineFix_run-01_blk-01.mat'); 
bl_data = load('decomp_signal_BaselineFix_date-04-18-2021.mat'); 
data = bl_data.decomp_signal; 
metadata = bl_data.metadata  ; 
%% Get FOI from this signal 

freqband = {'delta', 'theta','alpha','beta','lowgamma','highgamma'}; 

%initialize. 
FOI_final_data = {}; 
matrix_ROI = {}; 
bothhemi_matrix_ROI = {}; 

for j = 1:numel(freqband)
    freqname = freqband{j}; 
    switch freqname
        case 'delta'
            foi_band = 1:4;
        case 'theta'
            foi_band = 4:7;
        case 'alpha'
            foi_band = 8:12;
        case 'beta'
            foi_band = 12:18;
        case 'lowgamma'
            foi_band = 19:24;
        case 'highgamma'
            foi_band = 24:28;
    end
    
    % for a single freqband of interest, compute avg in that band. 
    for i = 1:numel(data)
        FOI_data(i).data = data(i).amplitude(foi_band,:);
        FOI_avg(i).data = squeeze(mean(FOI_data(i).data,1));
    end
    tmp = squeeze(cell2mat(struct2cell(FOI_avg)));
    FOI_final_data{j} = permute(tmp,[2,1]);
    
    % convert channels to ROI
    ch_dim = 1;
    [matrix_ROI{j},ROI_labels] = generate_ROI_from_ch(FOI_final_data{j},PatientID,ch_dim);
end


num_ROI = numel(ROI_labels);
%% 

first_lag = {} ; 
first_lag_allROI = {}; 
for j = 1:numel(freqband)
    [first_lag{j},first_lag_allROI{j}] = find_min_autocorr(matrix_ROI{j}, ROI_labels,PatientID,freqband{j}); 
end


%% Make histogram with distribution of values 

all_freqband_lags = cell2mat(first_lag); 
histogram(all_freqband_lags,12)

%% histograms for each FOI 

for i = 1:numel(freqband) 
    figure
    histogram(first_lag{i},8)
    title(freqband{i})
    xlabel('Number of samples')
    ylabel('Count')
    saveas(gcf, sprintf('FOI_%s_autocorr_%s.fig',freqband{i},PatientID))

end 

%% histogram for each ROI 

ROI = {}; 
for i = 1:numel(freqband)
    for j = 1:numel(ROI_labels)
    ROI{j}(i) = first_lag{1, i}(:,j);     
    end 
end 
%% ** uncomment 
for i = 1:numel(ROI_labels)
    all_ROI_lags = ROI{i}; 
    figure() 
    histogram(all_ROI_lags,12)
    title(ROI_labels{i});
    xlabel('Number of samples')
    ylabel('Count')
    saveas(gcf, sprintf('ROI_%s_autocorr_%s.fig',ROI_labels{i},PatientID))

end 
%% get median 

median_ROI = median(sort(all_freqband_lags,'ascend')); 
mean_ROI = mean(all_freqband_lags); 
%%% for 002 median is 2 seconds 
%%% for 001 mean is 6.6 seconds 

%%% for 001 median is 3.56 seconds 
%%% for 001 mean is 9.6 seconds 

%%% for 003 median is 4.04 seconds 
%%% for 003 mean is 4.43 seconds 

%% Get # of lags based on percentile chosen 

percentile = 70; 
time_for_lag = prctile(all_freqband_lags,percentile); 
disp(time_for_lag)

% Calculating # of trials that baseline rec will give after dealing with temporal correlation 
%% 
fs = 1000; 
lags = time_for_lag; 
lag_secs = lags/fs; 

trial_time = lag_secs+5; 

total_samples = length(FOI_avg(1).data); 
total_time_sec = total_samples/fs; 
total_time_min = total_time_sec/60; 

num_trials = total_time_sec/trial_time; 

disp(num_trials)



