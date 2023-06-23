%% run_mtspectrum_pow_from_timeseries_baseline_autocorr.m - Run this to get spectral
% power for baseline recordings of TRD participants 

% Additional info: N/A
% Inputs: N/A
% Outputs: Spectral power from Chronux multitaper fxn with frequencies
% Dependencies: Chronux toolbox, UH3 github repo 
% Sub-functions: mtspectrum_pow_from_timeseries_baseline, outputdir fxn

% Anusha Allawala, 10/2021
% Edits: Changed the way data is segmented and samples are skipped to 
% account for autocorrelation 01/31/22

%------------ START OF CODE ---------------% 
%% Define input parameters and load data 

clear all
PatientID = 'DBSTRD001'; 
run_num = '07'; 
blk_num = '01'; 
experiment_name = sprintf('BaselineFix_run-%s_blk-%s',run_num,blk_num);  

% Load single trial time series data 

%timeseriesfile = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Epoched Data/BaselineFix_run-%s_blk-%s/baseline_all_singletrial_%s.mat',...
    %PatientID,run_num,blk_num,experiment_name); %**** UNCOMMENT 
 
timeseries_data = load(timeseriesfile); 

fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 

ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 
metadata = timeseries_data.metadata;
metadata.processed.date = date; 
addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/Toolboxes/chronux_2_12.v03')); 


%% Compute PSD using multitaper for each contact configuration 
 
data = timeseries_data.output;
disp(size(data))
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
S = []; 
for i=1:num_ch
    for j = 1:num_trials 
        [S(j,i,:),f,Serr(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_data(j,i,:),[3,1,2])),params);
    end    
end 
%% convert to dB 

pow = 10*log(S_poststim);

%% save
    outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/%s/%s',PatientID,experiment_name, 'PSD'); 
    if ~exist(outputdir,'dir'); mkdir(outputdir); end   
    filename = sprintf('%s/%s_mtspectrum_pow_indivtr.mat',outputdir,experiment_name); 
       
    save(filename, 'pow', 'f', 'S', 'Serr_poststim','metadata') 
    disp('PSD for all channels saved');
    
%-------------- END OF CODE ---------------% 
