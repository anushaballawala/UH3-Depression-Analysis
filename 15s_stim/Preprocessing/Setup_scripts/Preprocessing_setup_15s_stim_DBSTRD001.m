%% script to complete initial preprocessing steps
% blank version to input parameters for each patient

%* run number will have to be changed for 002 and onwards 
function [outputdir,montageInfo,metadata,data,good_ch_idx,good_ch_labels] = Preprocessing_setup_15s_stim_DBSTRD001(subjecttype,PatientID,experiment_type,hemi,DBStarget,run,freq)
%% load data

neuraldatafolder = sprintf('E:/%s/%s/Raw Data/%s%sstim15s_%s_f%d',subjecttype, PatientID, hemi,DBStarget,run,freq);
neuraldatafile = sprintf('%s/sub-%s_task-%s%sstim15s_run-%s_blk-f%d_rawData',neuraldatafolder,PatientID, hemi,DBStarget,run,freq); 
neuraldata = load(neuraldatafile);

%% get time-series data

data = neuraldata.NS3.Data; % get data matrix from NS3 structure

%% get channel info

montageInfo = getMontageInfo(neuraldata); 
%%  patient ID

metadata.general.PatientID = PatientID; 

fprintf('Patient ID: %s\n', PatientID);

%% get behavioral task info from filename 
 
stimpattern = 'blk-\w*';
[start_idx, end_idx] = regexp(neuraldatafile,  stimpattern);
stim_type = extractBetween(neuraldatafile, start_idx+4, end_idx-8);
stim_type = stim_type{1}; 

% run = 'run-\d*';
% [start_idx, end_idx] = regexp(neuraldatafile,  run);
% run_num = extractBetween(neuraldatafile, start_idx+4, end_idx);
% run_num = run_num{1}; 

dbstarget = '-(l|r)(VCVS|SCC)';
[start_idx, end_idx] = regexp(neuraldatafile, dbstarget);
dbstarget = extractBetween(neuraldatafile, start_idx+1,end_idx); 
dbstarget = dbstarget{1}; 

stim_info = append(dbstarget,'_',stim_type); 

disp('assigned stimtype label')

metadata.general.Experiment = stim_info; 

fprintf('Experiment: %s\n', stim_info);
%% index channels to get good ones and convert to voltage 

%load good channels 
chinfo = load(sprintf('E:/%s/%s/Experiments/%s/Preprocessed data/%s/%s_goodch.mat',subjecttype,PatientID,experiment_type,stim_info,stim_info)); 
good_ch_idx = chinfo.good_ch_idx; 
good_ch_labels = chinfo.good_ch_labels; %** change back to labels 
%data = data(montageInfo.MacroContactIndices,:)*0.25; % convert from bit depth to voltage 
data = data*0.25; % convert from bit depth to voltage 

%% set up metadata structure

% note, might need to edit some stuff here! * 

metadata.general.montageInfo = montageInfo; 
metadata.general.SamplingRate = neuraldata.NS3.MetaTags.SamplingFreq;
metadata.general.DateTime = neuraldata.NS3.MetaTags.DateTime; 
metadata.general.Duration = neuraldata.NS3.MetaTags.DataDurationSec;
metadata.general.TotalSamples = neuraldata.NS3.MetaTags.DataPoints; %this could be wrong
metadata.files.RawDataFile = neuraldatafile;



%% output directory info

SubjectType = 'DBSTRD';
computer = 'Stronghold';
outputdir = make_dir_preprocessing_task(computer,stim_info,SubjectType,PatientID,experiment_type);

% %% deal with bipolar vs laplacian folder for outputs
% 
% outputdir = [outputdir '/Bipolar'];
if ~exist(outputdir,'dir'), mkdir(outputdir), end
end 

