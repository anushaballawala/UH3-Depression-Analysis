addpath(genpath('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/'))
clear all 
clc 
PatientID = 'DBSTRD006'; 
experiment_type = '15s_stim';
SubjectType = 'DBSTRD'; 
stim_info = 'rSCC_f130_blk-03'; 
stim_info_epoch = 'rSCC_03_f130'; 
disp(stim_info)
datatype = 'timeseries'; %'timeseries or wavelet'

% directory 
outputdir = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/%s/EXP/%s/Preprocessed Data/%s/',...
    PatientID,experiment_type,stim_info); 

% load epochtable with timestamps 
epochdatafile = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/%s/EXP/%s/Epochs/epoch_%s_%s.mat',...
    PatientID,experiment_type,PatientID,stim_info_epoch);
epochdata = load(epochdatafile);

% load preprocessed data
if strcmp(datatype,'timeseries') == 1
    preprocesseddatafile = sprintf('%s/referencedsig_%s.mat',outputdir,stim_info);
    preprocesseddata = load(preprocesseddatafile);
    neuraldata =  preprocesseddata.all_ref_signals;
    metadata = preprocesseddata.metadata;
elseif strcmp(datatype,'wavelet') == 1
    preprocesseddatafile = sprintf('%s/decomp_signal_%s.mat',outputdir,stim_info);
    %preprocesseddatafile = sprintf('%s/hilbert_forcfc%s.mat',outputdir,stim_info);
    preprocesseddata = load(preprocesseddatafile);
    neuraldata = preprocesseddata.decomp_signal;
    metadata = preprocesseddata.metadata;
    
end

%% Conditionbin 

%if strcmp(datatype,'timeseries') == 1
    output_and_metadata = conditionbin_time_series_15s_per_file(PatientID,metadata,neuraldata,epochdata,stim_info);
% elseif strcmp(datatype,'wavelet') == 1
%     conditionbin_15s_stim_perfile(PatientID,metadata,neuraldata,epochdata,stim_info)
%end

disp('done saving') 




