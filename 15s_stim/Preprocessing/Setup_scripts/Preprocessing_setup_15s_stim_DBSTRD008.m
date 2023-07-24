function [] = Preprocessing_setup_15s_stim_DBSTRD008(fullfile)
%% script to complete initial preprocessing steps

%function [outputdir,montageInfo,metadata,data,good_ch_idx,good_ch_labels] = Preprocessing_setup_15s_stim_DBSTRD006(subjecttype,PatientID,experiment_type,hemi,DBStarget,run,freq)
%% load data
%Stim15s_lSCC_f130

[PatientID,PatientType,block_file,...
    block_num,freq,hemi, DBStarget] = extract_stim_info(fullfile); 

if contains(fullfile, 'e1') 
    neuraldatafile = fullfile; 
else 
neuraldatafile = sprintf('sub-%s_task-Stim15s_run-%s%s_blk-f%s_subblk-%s_rawData.mat',...
PatientID, hemi, DBStarget, freq, block_file); 
end 
neuraldatafolder = sprintf('/users/aallawa1/data/TRD_Project/%s/%s/Raw Data/%s%sstim15s_%s_f%s/',...
    PatientType, PatientID, hemi, DBStarget, block_num, freq); 

neuraldata = load(neuraldatafile);

%% get time-series data

data = double(neuraldata.NS5.Data); % get data matrix from NS5 structure

%% get channel info

montageInfo = getMontageInfo(neuraldata); 
%%  patient ID

metadata.general.PatientID = PatientID; 

fprintf('Patient ID: %s\n', PatientID);

%% get behavioral task info from filename 
 
stim_info = sprintf('%s%s_f%s_blk-%s',hemi,DBStarget,freq,block_num); 

metadata.general.Experiment = stim_info; 

fprintf('Experiment: %s\n', stim_info);
%% Obtain good channel indices and convert signal to voltage 

ch_info = load(sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/%s/EXP/15s_stim/%s_15stim_goodch.mat',...
    PatientID, PatientID)); 
good_ch_idx = ch_info.good_ch_idx; 
good_ch_labels = ch_info.good_ch_labels;
%data = data(montageInfo.MacroContactIndices,:)*0.25; % convert from bit depth to voltage 
load(sprintf('%s_electrode_montage_15stim.mat',PatientID)); 

data = data*0.25; % convert from bit depth to voltage 

%% set up metadata structure

metadata.general.DateTime = neuraldata.NS5.MetaTags.DateTime; 
metadata.files.RawDataFile = neuraldatafile;
metadata.general.montageInfo = montageInfo; 
metadata.general.SamplingRate = neuraldata.NS5.MetaTags.SamplingFreq;
metadata.general.DateTime = neuraldata.NS5.MetaTags.DateTime; 

if isfield(neuraldata.NS5.MetaTags, 'DataDurationSec')
    metadata.general.Duration = neuraldata.NS5.MetaTags.DataDurationSec;
    metadata.general.TotalSamples = neuraldata.NS5.MetaTags.DataPoints; 
elseif isfield(neuraldata.NS5.MetaTags, 'DataPointsSec')
    metadata.general.Duration = neuraldata.NS5.MetaTags.DataPointsSec;
    metadata.general.TotalSamples = neuraldata.NS5.MetaTags.DataPoints; 
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

end 



