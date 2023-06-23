%% preprocessing_15s_stim_calculations_PARRM.m
% FUNCTION: main script for performing filtering, re-referencing,
% spectral decomposition for 15-second stimulation data  
% INPUTS:converted .mat files containing NS3 and NS5 structs 
% OUTPUTS: saves chXtime at multiple points in script (
% note: %make sure FT toolbox is not in your path, nanmean fxn screws up
% demeaning  
% note: order of preprocessing is changed since PARRM is implemented 
% Created for UH3 Depression Project by AA 06/2021

%function preprocessing_15s_calculations(montageInfo,ch,metadata,outputdir,good_ch_idx,good_ch_labels)

%% %%%%%%% temp data %%%%% - delete after troubleshooting script 
ch = data; 
srate_raw = 2000; 

%% extract info from metadata

behav_type = metadata.general.Experiment;  
ChannelLabel = montageInfo.ElectrodeLabels;  
contacts_perprobe  = montageInfo.ContactsPerProbe;  
num_chan = montageInfo.RecordedChannelCount ; %* check that this 256 for idx

%% initialize preprocessing metadata steps

metadata.preprocessing.Demeaned = false;
metadata.preprocessing.DownSampled = false;
metadata.preprocessing.Referenced = false;
metadata.preprocessing.BandPassed = false;
metadata.preprocessing.NotchFiltered = false;
metadata.preprocessing.ArtifactRemoval = false; 

%% get githash info - until we switch to Oscar this is probs useless 

[~,git_hash_string] = system('git rev-parse HEAD');
metadata.preprocessing.GitHash = git_hash_string;

%% 1. Set up timestamps from epoch table 

tbl = epochdata.tbl; 
trialStartTimes = tbl.Time;

%% 2. Set up variables w/relevant info for PARRM 

trialDuration = 15; % 15 seconds of stim ON 
fs = srate_raw; % sampling rate (no downsampling) 
stimFreq = 130; % 130 Hz stim frequency 
guessPer = fs/stimFreq; %Guess Period 
buffer = round(guessPer/2); 

% Filtered=highpass(data(1,:),1,2000);

Filtered = data(1,:); 
    for i=1:length(trialStartTimes)
        temp=Filtered(1,round(trialStartTimes(i)*fs)-fs:round(trialStartTimes(i)*fs)+trialDuration*fs-1+fs);
        [pks,locs]=findpeaks(-temp,'MinPeakDistance',0.9*guessPer);
        idx=kmeans(pks',2);
        m1=abs(mean(pks(idx==1)));
        m2=abs(mean(pks(idx==2)));
        if m1>m2
            locs=locs(idx==1);
        else
            locs=locs(idx==2);
        end
        temp=temp(1,locs(1)-buffer:locs(end)+buffer);
        if mod(i,5)==1

            Period=FindPeriodLFP(temp,[1,length(temp)-1],2000/130);
            PARRM=PeriodicFilter(Period,2000,0.04,0); %changed from 4000 and 0.01
            %PeriodicFilter(Period,TemporalWidth,PeriodWidth,OmitWidth,Direction)
        end
        baseline=movmean(temp,4000);
        temp=temp-baseline;
        start=round(trialStartTimes(i)*fs)-fs+locs(1)-buffer-1;
        Filtered(start:start+length(temp)-1)=((filter2(PARRM.',temp','same')-temp')./(1-filter2(PARRM.',ones(size(temp')),'same'))+temp')'+baseline;
    end

%Filtered=lowpass(Filtered,100,2000);

%% Bandpass data 

order = 4;
passBand = [0.5 250]; %0.5-250 Hz
BPsignal = zeros(size(Filtered,1),size(Filtered,2)); %initialize
for ichan = 1:size(Filtered,1)   %run Bandpass filter for all channels
    BPsignal(ichan,:) = butterworthBPFilter(Filtered(ichan,:),passBand,order,srate_raw);
end


%% 3. Demean data 
ch = BPsignal;  
tic
ch_baseline = nanmean(ch,2);
ch_baseline_demean = ch-repmat(ch_baseline, 1, length(ch));
num_ch = size(ch,1); 

disp('done demeaning') 
%add to metadata
metadata.preprocessing.Demeaned = true;
toc
%% 2. Bandpass data  
tic
order = 4;
passBand = [0.5 250]; %0.5-250 Hz
BPsignal = zeros(size(ch_baseline_demean,1),size(ch_baseline_demean,2)); %initialize
for ichan = 1:size(ch_baseline_demean,1)   %run Bandpass filter for all channels
    BPsignal(ichan,:) = butterworthBPFilter(ch_baseline_demean(ichan,:),passBand,order,srate_raw);
end

disp('done bandpass')
% metadata
metadata.preprocessing.PassBand = passBand;
metadata.preprocessing.BandPassOrder = order;
metadata.preprocessing.BandPassed = true;
toc

%% 3. Notch Filter 
tic

%notch = butterworthNotchFilter(BPsignal,60,4,srate_raw);
notch = butterworthNotchFilter(ch_baseline_demean,60,4,srate_raw);

notch = butterworthNotchFilter(notch,120,4,srate_raw);
notch = butterworthNotchFilter(notch,180,4,srate_raw); 
notch = butterworthNotchFilter(notch,240,4,srate_raw); 

notch_240 = notch;
clear notch; 

disp('finished notch filter')
% add to metadata
metadata.preprocessing.NotchFiltered = true;
metadata.preprocessing.NotchFilterMethod = 'butterworth';
metadata.preprocessing.NotchFilterFreqs = [60 120 180 240];
toc
%% save notch filter data for Matt 
tic
thisfile = sprintf('preproc_PARRM%s.mat', behav_type); 
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
save (fulldestination,'notch_240', 'metadata','-v7.3');

disp('--- saving notch filter data ---')
toc
%% 3a. Find dominant freq (001 and 002 not at exact stim frequency) 
plot_periodogram(data(1,:), srate_raw,'log') 
hold on 
plot_periodogram(Filtered(1,:), srate_raw,'log')
hold on 
plot_periodogram(notch_240(1,:),srate_raw,'log')
legend('raw','PARRM','additional proc') 

%% 4. PARRM (run only for one channel for now) 
% % 1. Define PARRM parameters. 
% stim_freq = 120; % supposed to be 130 but 120 because of cerestim problem
% guessPeriod = srate_raw/stim_freq ; %always samples by stim freq?
% span = [1,length(notch_240)-1]; 
% windowSize = 4500; % Width of the window in sample space for PARRM filter
% skipSize = 20; % Number of samples to ignore in each window in sample space 
% windowDirection = 'both';
% 
% % 2. Find Period of Signal. 
% Period = FindPeriodLFP(notch_240(101,:),span,guessPeriod); 
% % Plot to see over which timescale features are constant 
% plot(mod(1:length(notch_240(1,1:80000))-1,Period),diff(notch_240(1,1:80000)),'o'); 
% 
% % 3. Define window in period space for which samples will be averaged 
% PeriodDist = Period/120; % what is 120 * 
% 
% % 4. Create the linear filter
% PARRM=PeriodicFilter(Period,windowSize,PeriodDist,skipSize,windowDirection);
% 
% % 5. Filter using the linear filter and remove edge effects
% %stim_removed_ch = zeros(size(notch_240,1),size(notch_240,2)); %initialize
% %stim_removed_ch = []; 
%     stim_removed_ch =  (filter2(PARRM.',notch_240(1,:)','same')-notch_240(1,:)'./(1-filter2(PARRM.',ones(size(notch_240(1,:)')),'same'))+ notch_240(1,:)')'; 
% 
% % for i = 1:1 %num_ch
% %     stim_removed_ch(i,:) =  (filter2(PARRM.',notch_240(i,:)','same')-notch_240(i,:)'./(1-filter2(PARRM.',ones(size(notch_240(i,:)')),'same'))+ notch_240(i,:)')'; 
% % end 
% 
% %% 
% % metadata
% metadata.preprocessing.PARRM_Parameters.guessPeriod = guessPeriod;
% metadata.preprocessing.PARRM_Parameters.span = span;
% metadata.preprocessing.PARRM_Parameters.windowSize = windowSize;
% metadata.preprocessing.PARRM_Parameters.skipSize = skipSize;
% metadata.preprocessing.PARRM_Parameters.windowDirection = windowDirection;
% metadata.preprocessing.PARRM_Parameters.PeriodDist = PeriodDist;

%% 5. Downsample 
tic
srate_new = 1000; %new sampling rate - going to 1000 Hz from 2000 
r = 2; %factor for downsampling 
ds_Ch=[];

%downsample data to 1 kHz 
for k = 1:num_ch 
    ds_Ch(k,:) = decimate(ch_baseline_demean(k,:),r); %downsampling 
end 

disp('done downsampling')
% add to metadata
metadata.preprocessing.DownSampled = true;
metadata.preprocessing.New_SamplingRate = srate_new;
toc
%% 6. Visualize data before re-referencing 

%% 7. Re-reference 
tic
ref_method = 'Bipolar';
all_ref_signals = rereference(BPsignal,ref_method,contacts_perprobe,montageInfo.OrderedContactIndicesByProbe);  

disp('finished re-referencing')
%metadata
metadata.preprocessing.Referenced = true;
metadata.preprocessing.ReferenceMethod = ref_method;
toc
%% 8. Visualize data after re-referencing 

tic
end_ch = size(all_ref_signals,1); 
dir_for_figs = [outputdir '/Figures'];
if ~exist(dir_for_figs, 'dir')
    mkdir(dir_for_figs)
end

path_for_figs = [dir_for_figs '/Inspect_Ch%d_%s.png'];
close all 

disp('saved figures after referencing')
toc

%% Plot re-referenced signals to identify bad channels 



%% index of good and bad channels 

all_ch = 1:num_chan;
bad_ch = setdiff(all_ch,good_ch_idx);
bad_ch_labels = ChannelLabel(bad_ch); 

% metadata
metadata.preprocessing.GoodChannelsIdx = good_ch_idx;
metadata.preprocessing.GoodChannelLabels = good_ch_labels;


% PART 1 DONE ************************************************ %% 
%% choose method for spectral decomposition for neural data 
tic
data_cwt = all_ref_signals; 
end_ch = size(all_ref_signals,1); 
method = 'Wavelet'; 

if strcmp('Wavelet',method) == 1 
    wave_num = [3 3 4 4 4 5 5 5 6 6 6 6 7 7 8 8 9 9 10 10 11 11 12 13 14 16 17 18]; 
    freqs = [1 2 3 4 5 6 7 8 9 10 11 12 15 18 21 24 27 30 35 40 45 50 60 70 90 110 130 150]; %center frequencies 
    %metadata
    metadata.decomp_parameters.freqs = freqs; 
    metadata.decomp_parameters.wave_num = wave_num; 
    metadata.decomp_parameters.method = method; 
    
    [decomp_signal] = choose_decomp(data_cwt,method,freqs,wave_num,end_ch,srate_new);

elseif strcmp('Hilbert',method) == 1
    freqs = [1:12 15 18 21 24 27 30 35 40 45 50 60 70 90 110 130]; %center freqs
    % set up frequency range
    freqs_range = zeros(length(freqs), 2);
    freqs_range(1,:) = [0.5 1.5]; %set first freq band
    for i = 2:(length(freqs))
        % 50% overlap with previous frequency band
        half_range_length = freqs(i) - freqs(i - 1);
        freqs_range(i,1) = freqs(i) - half_range_length;
        freqs_range(i,2) = freqs(i) + half_range_length;
    end
    
    %metadata
    metadata.decomp_parameters.freqs = freqs;  
    metadata.decomp_parameters.freqs_range = freqs_range;
    metadata.decomp_parameters.method = method;
    
    [decomp_signal] = choose_decomp(data_cwt,method,freqs_range,[],end_ch,srate_new);
    
end 
  
disp('finished spectral decomposition')
toc
%% save
tic
if strcmp(method,'Wavelet') == 1
    thisfile = sprintf('decomp_signal_%s.mat', behav_type);
elseif strcmp(method,'Hilbert') == 1
    thisfile = sprintf('decomp_signal_Hilbert_%s.mat', behav_type);
end

fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory

metadata.files.PreprocessedData.StrongholdDirectory = outputdir;
metadata.files.PreprocessedData.FileName = thisfile;

save(fulldestination,'decomp_signal','metadata','-v7.3');

fprintf('finished saving %s\n', behav_type)
toc
%% 

clear
%end 





