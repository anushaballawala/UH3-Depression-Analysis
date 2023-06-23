%function[] = run_timeseries_mtspectrumc_power_003(timeseriesfile)

clear 
clc 
experiment_name = 'rVCVS_f130_blk-03'; 
timeseriesfile = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Epoched Data/%s/15s_stim_all_currdir_timeseries_singletrial_%s.mat',...
    experiment_name, experiment_name); 
timeseriesdata = load(timeseriesfile); 
%% 
PatientID = getPatientID(timeseriesfile); 

%% Load data 

fs = 1000;
ch_labels = deblank(timeseriesdata.metadata.preprocessing.GoodChannelLabels);

good_ch_idx = timeseriesdata.metadata.preprocessing.GoodChannelsIdx; 
summary_ch_info = timeseriesdata.metadata.preprocessing.ChannelTbl.csv_montage(good_ch_idx,:); 

% Add chronux to filepath.
addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/Toolboxes/chronux_2_12.v03')); 


%% Metadata 

metadata = timeseriesdata.metadata; 
metadata.processing.date = date; 

%% Run multitaper spectrum fxn 
   
%num_trials = size(timeseriesdata.output_and_metadata{1}{1}, 1); 
num_ch = size(timeseriesdata.output_and_metadata{1}{1},2); 
num_samples = size(timeseriesdata.output_and_metadata{1}{1}, 3); 
num_contact_config = length(timeseriesdata.output_and_metadata); 

%% concatenate contact config conditions 

if num_contact_config==2 
all_trials = vertcat(timeseriesdata.output_and_metadata{1}{1},...
    timeseriesdata.output_and_metadata{2}{1}); 
elseif num_contact_config ==3 
    all_trials = vertcat(timeseriesdata.output_and_metadata{1}{1},...
    timeseriesdata.output_and_metadata{2}{1},...
    timeseriesdata.output_and_metadata{3}{1}); 
end 

num_trials = size(all_trials,1); 
%% Define the different time windows of interest, in samples.
 
% pre_stim_win = 1:5000; 
% stim_win = 5001:20000; 
post_stim_win = 18100:23000; % artifact ends a bit after the 20 second mark 

disp('defined stim windows')

% metadata.stim_info.pre_stim_win = pre_stim_win; 
% metadata.stim_info.stim_win = stim_win1; 
metadata.stim_info.post_stim_win = post_stim_win; 
 %% Get different stim windows with indiv trial data 
 
%  indivtr_pre_stim_data_left = data_left(:,:,pre_stim_win); 
%  indivtr_stim1_data_left = data_left(:,:,stim_win1); 
%  indivtr_stim2_data_left = data_left(:,:,stim_win2); 
%  indivtr_stim3_data_left = data_left(:,:,stim_win3); 
 indivtr_poststim_data= all_trials(:,:,post_stim_win); 
 
 %% compute PSD for each trial 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
tic 
S_poststim = []; 

for i=1:num_ch
    for j = 1:num_trials 
%       [S_prestim(j,i,:),f_prestim,Serr_prestim(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data_left(j,i,:),[3,1,2])),params);
%          [S_stim1(j,i,:),f_stim1,Serr_stim1(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim1_data_left(j,i,:),[3,1,2])),params);
%         [S_stim2(j,i,:),f_stim2,Serr_stim2(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim2_data_left(j,i,:),[3,1,2])),params);
%         [S_stim3(j,i,:),f_stim3,Serr_stim3(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim3_data_left(j,i,:),[3,1,2])),params);  
        [S_poststim(j,i,:),f_poststim,Serr_poststim(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_poststim_data(j,i,:),[3,1,2])),params);
    end    
end 
toc 
 
%% 
tic 
outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/EXP/15s_stim/Processed Data/PSD',...
    PatientID);
 if ~exist(outputdir,'dir')
     mkdir(outputdir); end 
 
    filename = sprintf('%s_%s_mtspectrum_pow_indivtr.mat',PatientID,experiment_name); 
    filedir = sprintf('%s/%s',outputdir,filename); 
    save(filedir,'f_poststim','S_poststim', 'Serr_poststim','metadata','summary_ch_info',...
        'ch_labels','params') 

toc 
            
