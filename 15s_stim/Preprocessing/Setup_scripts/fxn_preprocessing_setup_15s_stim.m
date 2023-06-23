function [neuraldata] = fxn_preprocessing_setup_15s_stim(neuraldatafile)

% code fixes:
% 1. check that make_dir fxn works 
% 2. Check that datapoints sec bit works 
% 3. add good channel code after getting table 

%% load data

%%%%%% Use when troubleshooting %%%%%%
%neuraldatafile = dir('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD003/Raw Data/15s stim/Stim15s_VCVSf130_567/*.mat');
%neuraldata = load(neuraldatafile.name); 

neuraldata = load(neuraldatafile);
disp('loaded neuraldata file') 

%% get time-series data

data = neuraldata.NS3.Data; % get data matrix from NS3 structure

%% get channel info

montageInfo = getMontageInfo(neuraldata); 

%% get patient type and patient ID
disp('got montage info')
%% get patient ID 

%get patientID info 
start_idx = strfind(neuraldatafile, 'DBSTRD'); 
end_idx = start_idx + 8; 
PatientID = extractBetween(neuraldatafile, start_idx, end_idx);
PatientID = PatientID{1}; 
metadata.general.PatientID = PatientID; 

fprintf('Patient ID: %s\n', PatientID);

%% get behavioral task info from filename 

%set block and run pattern strings for 003 separately.
if strcmp(PatientID,'DBSTRD003') == 1
    block = 'blk-\w*';
    run = 'run-\w*\d';
else
    block = 'blk-\d*';
    run = 'run-\d*';
end 

[start_idx, end_idx] = regexp(neuraldatafile,  block);
block_num = extractBetween(neuraldatafile, start_idx, end_idx-8);
block_num = block_num{1}; 

%extract run num info 
[start_idx, end_idx] = regexp(neuraldatafile,  run);
run_num = extractBetween(neuraldatafile, start_idx, end_idx);
run_num = run_num{1}; 
%% get DBS target info 
if strcmp(PatientID,'DBSTRD003') == 1
    dbstarget_pattern = 'VCVS|SCC';  
    [start_idx, end_idx] = regexp(neuraldatafile, dbstarget_pattern);
    dbstarget = extractBetween(neuraldatafile, start_idx,end_idx); 
    dbstarget = dbstarget{1}; 
    run_info = 'f130'; 
    blk_info = block_num(5:end); 
    stim_info = strcat(dbstarget,'_',run_info,'_',blk_info); 
else
    dbstarget_pattern = '-(l|r)(VCVS|SCC)';
    [start_idx, end_idx] = regexp(neuraldatafile, dbstarget_pattern);
    dbstarget = extractBetween(neuraldatafile, start_idx+1,end_idx); 
    dbstarget = dbstarget{1}; 
    stim_info = append(dbstarget,'_',block_num); 

end

fprintf('Experiment: %s\n', stim_info);
metadata.general.Experiment = stim_info; 

fprintf('Experiment: %s\n', stim_info);

%% convert to voltage 

data = data(montageInfo.MacroContactIndices,:)*0.25; % convert from bit depth to voltage 

%% set up metadata structure

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

metadata.files.RawDataFile = neuraldatafile;

disp('set up metadata structure') 
%% append good channel info to metadata with ROI and gray matter/white matter info 
    
ch_info = load(sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Preprocessed Data/%s/%s_goodch.mat',PatientID,stim_info,stim_info)); 

% path = sprintf('/gpfs/data/dborton/TRD_Project/%s/%s/Documentation',PatientType,PatientID);
% addpath(path)
[summary_ch_info,summary_ch_good_idx_info] = create_ch_summary_tbl_DBSTRD(PatientID, ch_info.good_ch_idx); 

%**** get good_ch_labels here too from table???
good_ch_labels = summary_ch_good_idx_info.Label;
%% output directory info ***** 

computer = 'Oscar';

outputdir = make_dir_preprocessing_task(computer,stim_info,'DBSTRD',PatientID,'15s_stim');
%outputdir = make_dir_preprocessing_task(computer,stim_info,'DBSTRD',PatientID,'BaselineFix');

%% deal with bipolar vs laplacian folder for outputs

outputdir = [outputdir '/Bipolar'];
if ~exist(outputdir,'dir'), mkdir(outputdir), end

%% Finally, now run preprocessing 

ch = data;
preprocessing_15s_calculations(montageInfo,ch,metadata,outputdir,...
    ch_info.good_ch_idx,ch_info.good_ch_labels,summary_ch_good_idx_info)


end 


