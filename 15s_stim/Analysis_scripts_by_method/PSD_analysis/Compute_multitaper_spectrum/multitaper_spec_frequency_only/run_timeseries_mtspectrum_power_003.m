%function[] = run_timeseries_mtspectrumc_power_003(timeseriesfile)

clear 
clc 

timeseriesfile = dir('/gpfs/data/dborton/TRD_Project/DBSTRD/DBSTRD003/Experiments/15s_stim/Epoched Data/*/15s_stim_all_currdir_timeseries_singletrial_SCC_f130_567.mat'); 
load(sprintf('%s/%s',timeseriesfile.folder,timeseriesfile.name)); 
timeseriesfile = sprintf('%s/%s',timeseriesfile.folder,timeseriesfile.name)
%% 
% get info from filename 
subjectID_pattern= '/DBSTRD00*\w*'; 
[start_idx_subid end_idx_subid] = regexp(timeseriesfile,subjectID_pattern);
PatientID = extractBetween(timeseriesfile,start_idx_subid+1,end_idx_subid);
PatientID = PatientID{1}; 

braintarget = '/(VCVS|SCC)_';
[start_idx_br, end_idx_br] = regexp(timeseriesfile,braintarget);
DBS_target = extractBetween(timeseriesfile,start_idx_br+1,end_idx_br-1);
DBS_target = DBS_target{1}; 

block_num = timeseriesfile(end-11:end-4); 

experiment_name = sprintf('%s_%s',DBS_target,block_num);  

% % Load single trial time series data 
% timeseriesfile = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
%     PatientID,hemi,DBS_target,stim_freq,hemi,DBS_target,stim_freq); 


%% Load data 
timeseries_data = load(timeseriesfile); 

if strcmp(experiment_name,'rVCVS_f130')==1 && strcmp(PatientID,'DBSTRD002') ==1 
    fs = timeseries_data.metadata.metadata_file1.preprocessing.New_SamplingRate;
    ch_labels = deblank(timeseries_data.metadata.metadata_file1.preprocessing.GoodChannelLabels);
else
    fs = timeseries_data.metadata.preprocessing.New_SamplingRate;
    ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels);
end

% Add chronux to filepath.
addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/Toolboxes/chronux_2_12.v03')); 

%% Run function for all contact configurations 
disp(experiment_name) 

switch DBS_target
    case 'VCVS'       
        lVCVS = timeseries_data.output_and_metadata{1, 1}{1, 1};
        rVCVS = timeseries_data.output_and_metadata{1, 2}{1, 1};
        data_left = lVCVS;
        data_right = rVCVS;
    case 'SCC'       
        lSCC = timeseries_data.output_and_metadata{1, 1}{1, 1};
        rSCC = timeseries_data.output_and_metadata{1, 2}{1, 1};      
        data_left = lSCC;
        data_right = rSCC; 
end


%% Metadata 

metadata = timeseries_data.metadata; 
metadata.processing.date = date; 

%% Run multitaper spectrum fxn 
   
num_trials = size(data_left,1)
num_ch = size(data_left, 2)
num_samples = size(data_left, 3)
%% Define the different time windows of interest, in samples.
 
pre_stim_win = 1:5000; 
stim_win = 5001:20000; 
stim_win1 = 5000:10000; 
stim_win2 = 10000:15000; 
stim_win3 = 15000:20000; 
post_stim_win = 20600:25600; % artifact ends a bit after the 20 second mark 

disp('defined stim windows')

metadata.stim_info.pre_stim_win = pre_stim_win; 
metadata.stim_info.stim_win1 = stim_win1; 
metadata.stim_info.stim_win2 = stim_win2; 
metadata.stim_info.stim_win3 = stim_win3; 
metadata.stim_info.post_stim_win = post_stim_win; 
 %% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data_left = data_left(:,:,pre_stim_win); 
 indivtr_stim1_data_left = data_left(:,:,stim_win1); 
 indivtr_stim2_data_left = data_left(:,:,stim_win2); 
 indivtr_stim3_data_left = data_left(:,:,stim_win3); 
 indivtr_poststim_data_left= data_left(:,:,post_stim_win); 
 
 %% compute PSD for each trial 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 
prestim_color = [0.6706    0.1882    0.1882]; 
stim_color = [0.6549    0.3765    0.7098]; 
poststim_color = [0.5020    0.5020    0.5020]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
S_prestim = []; 
S_stim1 = []; S_poststim = []; 
for i=1:num_ch
    for j = 1:num_trials 
    %prestim 
      [S_prestim(j,i,:),f_prestim,Serr_prestim(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data_left(j,i,:),[3,1,2])),params);
         [S_stim1(j,i,:),f_stim1,Serr_stim1(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim1_data_left(j,i,:),[3,1,2])),params);
        [S_stim2(j,i,:),f_stim2,Serr_stim2(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim2_data_left(j,i,:),[3,1,2])),params);
        [S_stim3(j,i,:),f_stim3,Serr_stim3(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim3_data_left(j,i,:),[3,1,2])),params);  
        [S_poststim(j,i,:),f_poststim,Serr_poststim(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_poststim_data_left(j,i,:),[3,1,2])),params);
    end    
end 

outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD',...
    PatientID);
 if ~exist(outputdir,'dir')
     mkdir(outputdir); end 
 
    filename = sprintf('%s_l%s_mtspectrum_pow_indivtr.mat',PatientID,experiment_name); 
    filedir = sprintf('%s/%s',outputdir,filename); 
    save(filedir,'f_prestim', 'S_prestim', 'f_stim1', 'S_stim1',...
        'S_stim2','S_stim3','f_poststim','S_poststim','Serr_prestim', 'Serr_stim1', 'Serr_stim2',...
        'Serr_stim3', 'Serr_poststim','metadata') 


            
%% Now recompute for the right %%%%%%%%%%%%%

%% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data_right = data_right(:,:,pre_stim_win); 
 indivtr_stim1_data_right = data_right(:,:,stim_win1); 
 indivtr_stim2_data_right = data_right(:,:,stim_win2); 
 indivtr_stim3_data_right = data_right(:,:,stim_win3); 
 indivtr_poststim_data_right= data_right(:,:,post_stim_win); 
 
 %% compute PSD for each trial 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 
prestim_color = [0.6706    0.1882    0.1882]; 
stim_color = [0.6549    0.3765    0.7098]; 
poststim_color = [0.5020    0.5020    0.5020]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
S_prestim = []; 
S_stim1 = []; S_poststim = []; 
for i=1:num_ch
    for j = 1:num_trials 
    %prestim 
      [S_prestim(j,i,:),f_prestim,Serr_prestim(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data_right(j,i,:),[3,1,2])),params);
         [S_stim1(j,i,:),f_stim1,Serr_stim1(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim1_data_right(j,i,:),[3,1,2])),params);
        [S_stim2(j,i,:),f_stim2,Serr_stim2(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim2_data_right(j,i,:),[3,1,2])),params);
        [S_stim3(j,i,:),f_stim3,Serr_stim3(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_stim3_data_right(j,i,:),[3,1,2])),params);  
        [S_poststim(j,i,:),f_poststim,Serr_poststim(:,:,i)] = mtspectrumc(squeeze(permute(indivtr_poststim_data_right(j,i,:),[3,1,2])),params);
    end    
end 

outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD',...
    PatientID);
 if ~exist(outputdir,'dir')
     mkdir(outputdir); end 
 
    filename = sprintf('%s_r%s_mtspectrum_pow_indivtr.mat',PatientID,experiment_name); 
    filedir = sprintf('%s/%s',outputdir,filename); 
    save(filedir,'f_prestim', 'S_prestim', 'f_stim1', 'S_stim1',...
        'S_stim2','S_stim3','f_poststim','S_poststim','Serr_prestim', 'Serr_stim1', 'Serr_stim2',...
        'Serr_stim3', 'Serr_poststim','metadata') 


             

%end 