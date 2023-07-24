% set up directory 
clear 
clc 

PatientID = 'DBSTRD010'; 
filename = 'bandpasssig_rSCC_f130_blk-03.mat';
subfolder = 'rSCC_f130_blk-03';
outputdir = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/%s/EXP/15s_stim/Preprocessed Data/%s',...
    PatientID,subfolder); 
fullfile = sprintf('%s/%s', outputdir, filename); 

%load data 
neuraldata = load(fullfile); 

%get metadata. 
metadata = neuraldata.metadata; 
behav_type = metadata.general.Experiment; 
montageInfo = metadata.general.montageInfo; 
ChannelLabel = montageInfo.ElectrodeLabels;  
contacts_perprobe  = montageInfo.ContactsPerProbe;  
num_chan = montageInfo.RecordedChannelCount ;

%% reference (default method Bipolar, but other methods can be used) 
tic

BPsignal_noempty = neuraldata.BPsignal;  
empty_idx = montageInfo.EmptyChannelIndices; 
BPsignal_noempty(empty_idx,:) = []; 

% get good channel and other montage info about channkels. 
ch_info = load(sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/%s/EXP/15s_stim/%s_15stim_goodch.mat',...
    PatientID, PatientID)); 
good_ch_idx = ch_info.good_ch_idx; 
good_ch_labels = ch_info.good_ch_labels;

%load csv montage with grey matter info. 
load(sprintf('%s_electrode_montage_15stim.mat',PatientID)); 

summary_ch_info = channel_info; 
%remove empty channels. 
summary_ch_info.csv_montage(empty_idx,:) = []; 
summary_ch_info.NSP_labels(empty_idx) = []; 
try 
    ROI_info = summary_ch_info.csv_montage.ROI; 
    matter_info = summary_ch_info.csv_montage.Matter; 
catch 
    ROI_info = summary_ch_info.csv_montage.ROI_vis; 
    matter_info = summary_ch_info.csv_montage.Matter_vis; 
end 

%% re-reference 
ref_method = 'Bipolar';
all_ref_signals = rereference(BPsignal_noempty,ref_method,contacts_perprobe,montageInfo.OrderedContactIndicesByProbe,...  
montageInfo.ElectrodeLabels,ROI_info,matter_info); 
disp('finished re-referencing')
%metadata
metadata.preprocessing.Referenced = true;
metadata.preprocessing.ReferenceMethod = ref_method;
toc

%% index of good and bad channels 

% metadata
metadata.preprocessing.GoodChannelsIdx = good_ch_idx;
metadata.preprocessing.GoodChannelLabels = good_ch_labels;
metadata.preprocessing.ChannelTbl = summary_ch_info; 


filename = metadata.general.Experiment ;   
ElectrodeID = montageInfo.MacroContactIndices; %Channel number 


%% save referenced signals 

%index referenced signals 
all_ref_signals = all_ref_signals(good_ch_idx',:); 

clear fullfile
tic 
thisfile = sprintf('referencedsig_%s.mat', behav_type);
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
save (fulldestination,'all_ref_signals', 'metadata','-v7.3')

disp('saved re-referenced signals')
toc 











