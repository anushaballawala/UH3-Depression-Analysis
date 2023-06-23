
%% script to complete initial preprocessing steps

%function [outputdir,montageInfo,metadata,data,good_ch_idx,good_ch_labels] = Preprocessing_setup_15s_stim_DBSTRD006(subjecttype,PatientID,experiment_type,hemi,DBStarget,run,freq)
%% load data
%Stim15s_lSCC_f130
clear 
PatientType = 'DBSTRD'; 
PatientID = 'DBSTRD006'; 
hemi = 'r'; freq = '130'; 

DBStarget = 'SCC'; 
block_full = '01'; block_num = '1-2'; 

neuraldatafile = sprintf('sub-%s_task-Stim15s_run-%s%s_blk-f%s_subblk-%s.rawData.mat',...
PatientID, hemi, DBStarget, freq, block_num); 
neuraldatafolder = sprintf('/users/aallawa1/data/TRD_Project/%s/%s/Raw Data/%s%sstim15s_%s_f%s/',...
    PatientType, PatientID, hemi, DBStarget, block_full, freq); 

neuraldata = load(neuraldatafile);

%% get time-series data

data = double(neuraldata.NS3.Data); % get data matrix from NS3 structure

%% get channel info

montageInfo = getMontageInfo(neuraldata); 
%%  patient ID

metadata.general.PatientID = PatientID; 

fprintf('Patient ID: %s\n', PatientID);

%% get behavioral task info from filename 
 
stim_info = sprintf('%s%s_f%s_blk-%s',hemi,DBStarget,freq,block_full); 

metadata.general.Experiment = stim_info; 

fprintf('Experiment: %s\n', stim_info);
%% Obtain good channel indices and convert signal to voltage 

ch_info = load('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/DBSTRD006_15stim_goodch.mat'); 
good_ch_idx = ch_info.good_ch_idx; 
good_ch_labels = ch_info.good_ch_label;
%data = data(montageInfo.MacroContactIndices,:)*0.25; % convert from bit depth to voltage 
load(sprintf('%s_electrode_montage_15stim.mat',PatientID)); 

data = data*0.25; % convert from bit depth to voltage 

%% set up metadata structure

metadata.general.DateTime = neuraldata.NS3.MetaTags.DateTime; 
metadata.files.RawDataFile = neuraldatafile;
metadata.general.montageInfo = montageInfo; 
metadata.general.SamplingRate = neuraldata.NS3.MetaTags.SamplingFreq;
metadata.general.DateTime = neuraldata.NS3.MetaTags.DateTime; 

if isfield(neuraldata.NS3.MetaTags, 'DataDurationSec')
    metadata.general.Duration = neuraldata.NS3.MetaTags.DataDurationSec;
    metadata.general.TotalSamples = neuraldata.NS3.MetaTags.DataPoints; 
elseif isfield(neuraldata.NS3.MetaTags, 'DataPointsSec')
    metadata.general.Duration = neuraldata.NS3.MetaTags.DataPointsSec;
    metadata.general.TotalSamples = neuraldata.NS3.MetaTags.DataPoints; 
end 

disp('set up metadata structure') 

%% output directory info

computer = 'Oscar'; experiment = '15s_stim'; 
outputdir = make_dir_preprocessing_task(computer,stim_info,PatientType,PatientID,experiment); 
% 
% outputdir = [outputdir '/Bipolar'];
if ~exist(outputdir,'dir'), mkdir(outputdir), end
%end 

%% run preprocessing 

preprocessing_15s_calculations(montageInfo,data,metadata,outputdir,good_ch_idx,good_ch_labels,channel_info)





