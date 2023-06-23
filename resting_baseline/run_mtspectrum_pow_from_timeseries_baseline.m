%% run_mtspectrum_pow_from_timeseries_baseline.m - Run this to get spectral
% power for baseline recordings of TRD participants 

% Additional info: N/A
% Inputs: N/A
% Outputs: Spectral power from Chronux multitaper fxn with frequencies
% Dependencies: Chronux toolbox, UH3 github repo 
% Sub-functions: mtspectrum_pow_from_timeseries_baseline, outputdir fxn

% Anusha Allawala, 10/2021
% Edits: Changed directories for Oscar 11/21 - AA


%------------ START OF CODE ---------------% 
%% Define input parameters and load data 

clear all 
PatientID = 'DBSTRD001'; 
run_num = '07'; 
blk_num = '01'; 
experiment_name = sprintf('BaselineFix_run-%s_blk-%s',run_num,blk_num);  

% Load single trial time series data 

timeseriesfile = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Epoched Data/BaselineFix_run-%s_blk-%s/baseline_all_singletrial_%s.mat',...
    PatientID,run_num,blk_num,experiment_name); 
 
timeseries_data = load(timeseriesfile); 

fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 

ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 
metadata = timeseries_data.metadata;
metadata.processed.date = date; 
addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/Toolboxes/chronux_2_12.v03')); 

%% Get data matrices for each contact configuration 

disp(experiment_name) 
% switch PatientID
%     case 'DBSTRD001'
%         load(sprintf('ROI_labels_%s.mat', PatientID)); 
%         
%     case 'DBSTRD002'
%         load(sprintf('ROI_labels_%s.mat', PatientID)); 
%  
%     case 'DBSTRD003'
%         disp('load ROI labels'); 
% 
% end

%% Compute PSD using multitaper for each contact configuration 
 
data = timeseries_data.output;
disp(size(data))
switch PatientID 
    case 'DBSTRD001'
        data = permute(data,[2,3,1]); %*** make sure its the same size 
    case 'DBSTRD002'
        data = permute(data,[2,1,3]); 
end 

[pow_prestim,pow_stim1,pow_stim2,pow_stim3,...
    pow_poststim,f_prestim, S_prestim, f_stim1, S_stim1,...
    S_stim2,S_stim3,f_poststim,S_poststim,...
    Serr_prestim, Serr_stim1, Serr_stim2,...
    Serr_stim3, Serr_poststim,metadata] = mtspectrum_pow_from_timeseries_baseline(PatientID,...
    experiment_name,timeseriesfile, timeseries_data, fs, ch_labels,...
    data,metadata);
       
%% save
    outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/%s/%s',PatientID,experiment_name, 'PSD'); 
    if ~exist(outputdir,'dir'); mkdir(outputdir); end   
    filename = sprintf('%s/%s_mtspectrum_pow_indivtr.mat',outputdir,experiment_name); 
       
    save(filename, 'pow_prestim','pow_stim1','pow_stim2','pow_stim3',...
    'pow_poststim','f_prestim', 'S_prestim', 'f_stim1', 'S_stim1',...
    'S_stim2','S_stim3','f_poststim','S_poststim','Serr_prestim', 'Serr_stim1', 'Serr_stim2',...
    'Serr_stim3', 'Serr_poststim','metadata') 
    disp('PSD for all channels saved');
    
%-------------- END OF CODE ---------------% 
