% FUNCTION: main script for performing filtering, re-referencing,
% spectral decomposition for 15-second stimulation data  
% INPUTS:converted .mat files containing NS3 and NS5 structs 
% OUTPUTS: saves chXtime at multiple points in script (notch filtering,
% re-referencing, bandpass), and 1x ch struct with freqsxtime after waveleting or Hilbert transform   
% Dependencies: none 
% note: %make sure FT toolbox is not in your path, nanmean fxn screws up
% demeaning  
% Created for UH3 Depression Project by AA 11/2020
% Modified by AA 11/2020 : changed source for ChannelLabels 
% Modified by AA 03/22: changes inputs for summary_ch_tbl to include ROI
% info 
%%
function preprocessing_15s_calculations(montageInfo,ch,metadata,outputdir,good_ch_idx,good_ch_labels,summary_ch_info)

%% extract info from metadata

behav_type = metadata.general.Experiment;  
ChannelLabel = montageInfo.ElectrodeLabels;  
contacts_perprobe  = montageInfo.ContactsPerProbe;  
num_chan = montageInfo.RecordedChannelCount ;

try 
    ROI_info = summary_ch_info.csv_montage.ROI; 
    matter_info = summary_ch_info.csv_montage.Matter; 
catch 
    ROI_info = summary_ch_info.csv_montage.ROI_vis; 
    matter_info = summary_ch_info.csv_montage.Matter_vis; 
end 

%% initialize preprocessing metadata steps

metadata.preprocessing.Demeaned = false;
metadata.preprocessing.DownSampled = false;
metadata.preprocessing.Referenced = false;
metadata.preprocessing.BandPassed = false;
metadata.preprocessing.NotchFiltered = false;
metadata.preprocessing.ArtifactRemoval = false; 

%% Demeaning data to see true mean signal %%
tic
ch_baseline = nanmean(ch,2);
ch_baseline_demean = ch-repmat(ch_baseline, 1, length(ch));
num_ch = size(ch,1); 

disp('done demeaning') 
%add to metadata
metadata.preprocessing.Demeaned = true;
toc

%% downsample/decimate data 
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
%% Notch filter to remove line noise at 60 and 120,240(harmonics) 
tic
notch = butterworthNotchFilter(ds_Ch,60,4,srate_new);
notch = butterworthNotchFilter(notch,120,4,srate_new);
notch = butterworthNotchFilter(notch,180,4,srate_new); 
notch = butterworthNotchFilter(notch,240,4,srate_new); 

notch_240 = notch;
clear notch; 

disp('finished notch filter')
% add to metadata
metadata.preprocessing.NotchFiltered = true;
metadata.preprocessing.NotchFilterMethod = 'butterworth';
metadata.preprocessing.NotchFilterFreqs = [60 120 180 240];
toc

%% save data 
% tic
% thisfile = sprintf('notchsig_%s.mat', behav_type); 
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% 
% save (fulldestination,'notch_240', 'metadata','-v7.3');
% 
% disp('--- saving notch filter data ---')
% toc

%% Bandpass data  
tic
signal =  notch_240; 
clear notch_240 ; 
order = 4;
passBand = [1 250]; %0.5-250 Hz
BPsignal = zeros(size(signal,1),size(signal,2)); %initialize
for ichan = 1:size(signal,1)   %run Bandpass filter for all channels
    BPsignal(ichan,:) = butterworthBPFilter(signal(ichan,:),passBand,order,srate_new);
end

disp('done bandpass')
% metadata
metadata.preprocessing.PassBand = passBand;
metadata.preprocessing.BandPassOrder = order;
metadata.preprocessing.BandPassed = true;
toc
%% save bandpass signal data 
% tic
% thisfile = sprintf('bandpasssig_%s.mat', behav_type); 
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save (fulldestination,'BPsignal','metadata','-v7.3');
% 
% disp('--- saving bandpass data ---')
% toc

%% reference (default method Bipolar, but other methods can be used) 
tic

BPsignal_noempty = BPsignal;  
empty_idx = [157,208,209]; 
BPsignal_noempty(empty_idx,:) = []; 
summary_ch_info.csv_montage(empty_idx,:) = []; 
summary_ch_info.NSP_labels(empty_idx) = []; 

ref_method = 'Bipolar';
all_ref_signals = rereference(BPsignal_noempty,ref_method,contacts_perprobe,montageInfo.OrderedContactIndicesByProbe,...  
montageInfo.ElectrodeLabels,ROI_info,matter_info); 
disp('finished re-referencing')
%metadata
metadata.preprocessing.Referenced = true;
metadata.preprocessing.ReferenceMethod = ref_method;
toc
%% Plot re-referenced signals to identify bad channels 
tic
end_ch = size(all_ref_signals,1); 
dir_for_figs = [outputdir '/Figures'];
if ~exist(dir_for_figs, 'dir')
    mkdir(dir_for_figs)
end
filename = metadata.general.Experiment ;   
path_for_figs = [dir_for_figs '/Inspect_Ch%d_%s.fig'];

ElectrodeID = montageInfo.MacroContactIndices; %Channel number 
% 
% plot_PSD_withref(ds_Ch,all_ref_signals,srate_new,1,end_ch,filename,montageInfo.ElectrodeLabels,...
%     ElectrodeID,path_for_figs)
% 
% disp('channels after re-referencing, no empties')
% 
% 
% close all 
% 
% disp('saved figures after referencing')
% toc


%% index of good and bad channels 

% metadata
metadata.preprocessing.GoodChannelsIdx = good_ch_idx;
metadata.preprocessing.GoodChannelLabels = good_ch_labels;
metadata.preprocessing.ChannelTbl = summary_ch_info; 
%% save referenced signals 

%index referenced signals 
all_ref_signals = all_ref_signals(good_ch_idx',:); 

tic
thisfile = sprintf('referencedsig_%s.mat', behav_type);
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
save (fulldestination,'all_ref_signals', 'metadata','-v7.3')

disp('saved re-referenced signals')
toc
disp('----- preprocessing part 1 complete -----')

% %% choose method for spectral decomposition for neural data 
% tic
% data_cwt = all_ref_signals; 
% end_ch = size(all_ref_signals,1); 
% method = 'Wavelet'; 
% 
% if strcmp('Wavelet',method) == 1 
%     wave_num = [3 3 4 4 4 5 5 5 6 6 6 6 7 7 8 8 9 9 10 10 11 11 12 13 14 16 17 18]; 
%     freqs = [1 2 3 4 5 6 7 8 9 10 11 12 15 18 21 24 27 30 35 40 45 50 60 70 90 110 130 150]; %center frequencies 
%     %metadata
%     metadata.decomp_parameters.freqs = freqs; 
%     metadata.decomp_parameters.wave_num = wave_num; 
%     metadata.decomp_parameters.method = method; 
%     
%     [decomp_signal] = choose_decomp(data_cwt,method,freqs,wave_num,end_ch,srate_new);
% 
% elseif strcmp('Hilbert',method) == 1
%     freqs = [1:12 15 18 21 24 27 30 35 40 45 50 60 70 90 110 130]; %center freqs
%     % set up frequency range
%     freqs_range = zeros(length(freqs), 2);
%     freqs_range(1,:) = [0.5 1.5]; %set first freq band
%     for i = 2:(length(freqs))
%         % 50% overlap with previous frequency band
%         half_range_length = freqs(i) - freqs(i - 1);
%         freqs_range(i,1) = freqs(i) - half_range_length;
%         freqs_range(i,2) = freqs(i) + half_range_length;
%     end
%     
%     %metadata
%     metadata.decomp_parameters.freqs = freqs;  
%     metadata.decomp_parameters.freqs_range = freqs_range;
%     metadata.decomp_parameters.method = method;
%     
%     [decomp_signal] = choose_decomp(data_cwt,method,freqs_range,[],end_ch,srate_new);
%     
% end 
%   
% disp('finished spectral decomposition')
% toc
% %% save
% tic
% if strcmp(method,'Wavelet') == 1
%     thisfile = sprintf('decomp_signal_%s.mat', behav_type);
% elseif strcmp(method,'Hilbert') == 1
%     thisfile = sprintf('decomp_signal_Hilbert_%s.mat', behav_type);
% end
% 
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% 
% metadata.files.PreprocessedData.StrongholdDirectory = outputdir;
% metadata.files.PreprocessedData.FileName = thisfile;
% 
% save(fulldestination,'decomp_signal','metadata','-v7.3');
% 
% fprintf('finished saving %s\n', behav_type)
% toc
% %% 
% 
% clear
end 
