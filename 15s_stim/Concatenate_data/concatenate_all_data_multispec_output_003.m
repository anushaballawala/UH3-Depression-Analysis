%% Script: Concatenate_all_data_multispec_output_003.m

%Function: this script takes all the spectral analysis outputs from the 
% Chronux mtspectrum.c function (saved on the Oscar server) and
% concatenates all of the data from SCC, VCVS stim, plus baseline
% recordings (5 min pre experiment). It also concatenates all of the labels
% for stim state, DBS target, current configuration. This has only been
% written for 130 Hz stim experiments. 

%  Outputs: 
% 1) single 3D matrix that is trials x channels x frequencies. 
%  This output is used for stats and viz analysis. 
% 2) electrode labels, stim_state (prestim, poststim, etc.), DBS target
% label for each trial. 

% Note: I RUN THIS SCRIPT LOCALLY, NOT ON OSCAR. 

% Warning: script is long. Plz be v careful about making any changes. 
% AA 02/07/22

clear 
clc  
% add paths 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/PSD')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD001/Experiments/BaselineFix/Processed Data/')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD002/Experiments/15s_stim/Processed Data/PSD')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD002/Experiments/BaselineFix/Processed Data/')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD003/Experiments/'));


% ------------------------ START CODE -------------------------------% 
%% Load stim data. 

PatientID = 'DBSTRD003'; 

lSCC_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD/*lSCC*.mat',PatientID));
rSCC_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD/*rSCC*.mat',PatientID));
lVCVS_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD/*lVCVS*.mat',PatientID));
rVCVS_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD/*rVCVS*.mat',PatientID));

%% 

for i = 1:numel(lSCC_files)
    lSCC_data{i} = load(sprintf('%s/%s',lSCC_files(i).folder, lSCC_files(i).name)); 
    rSCC_data{i} = load(sprintf('%s/%s',rSCC_files(i).folder, rSCC_files(i).name)); 
    lSCC_name{i} = lSCC_files(i).name;
    rSCC_name{i} = rSCC_files(i).name;
    
end 

for i = 1:numel(lVCVS_files)
    lVCVS_data{i} = load(sprintf('%s/%s',lVCVS_files(i).folder, lVCVS_files(i).name)); 
    lVCVS_name{i} = lVCVS_files(i).name;
    rVCVS_data{i} = load(sprintf('%s/%s',rVCVS_files(i).folder, rVCVS_files(i).name));  
    rVCVS_name{i} = rVCVS_files(i).name;
end

%% Load baseline data. 

switch PatientID 
    case 'DBSTRD001' 
        bl_dir = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/BaselineFix_run-07_blk-01/PSD/*.mat',PatientID));       
        bl_data = load(sprintf('%s/%s',bl_dir(2).folder,bl_dir(2).name)); 
    
    case 'DBSTRD002' 
        bl_dir = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/BaselineFix_run-01_blk-01/PSD/*.mat',PatientID)); 
        bl_data = load(sprintf('%s/%s',bl_dir.folder,bl_dir.name)); 
        
    case 'DBSTRD003'
        bl_dir = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/BaselineFix_date-04-18-2021/PSD/*.mat',PatientID)); 
        bl_data = load(sprintf('%s/%s',bl_dir.folder,bl_dir.name)); 
end 

%% take out channels that are in baselinefix and not in 15s stim 

%load ROI labels 

ROI_labels = load('/Users/anushaallawala/Data/DBSTRD/15s_stim/Electrodes Info/ROI_labels_mat/DBSTRD003_ROI_labels.mat');

bl_remove_idx = [9, 12, 17:18,29, 52,53, 60:61,71] ;
%remove channels from bl 
baseline_good_ch = bl_data.metadata.preprocessing.GoodChannelLabels; 
baseline_good_ch(bl_remove_idx) = []; 

stim_remove_idx = [6, 29, 49, 50] ;
%remove channels from stim good ch labels 
good_ch_labels = ROI_labels.good_ch_labels;
good_ch_labels(stim_remove_idx) = [] ; 

% take_out_stim = [6;29;49;51];
% 
% take_out_bl = [9;12;17;18;29;52;60;61;71]; 
% 
keep_stim_idx = [1:5,7:28,30:48,51:65];
keep_bl_idx = [1:8,10:11,13:16,19:28,30:51,54:59,62:70];

summary_ch_info = ROI_labels.summary_ch_info  ; 
%remove channels from summary table 
summary_ch_info(stim_remove_idx,:) = []; 
%% get spectral data 

for i = 1:numel(lSCC_files) 
    
    lSCC_prestim{i} = (lSCC_data{1,i}.S_prestim(:,keep_stim_idx,:));  
    rSCC_prestim{i} = (rSCC_data{1,i}.S_prestim(:,keep_stim_idx,:)); 
    
    lSCC_stim1{i} = (lSCC_data{1,i}.S_stim1(:,keep_stim_idx,:));  
    rSCC_stim1{i} = (rSCC_data{1,i}.S_stim1(:,keep_stim_idx,:));  
    
    lSCC_stim2{i} = (lSCC_data{1,i}.S_stim2(:,keep_stim_idx,:));  
    rSCC_stim2{i} = (rSCC_data{1,i}.S_stim2(:,keep_stim_idx,:));   
    
    lSCC_stim3{i} = (lSCC_data{1,i}.S_stim3(:,keep_stim_idx,:));  
    rSCC_stim3{i} = (rSCC_data{1,i}.S_stim3(:,keep_stim_idx,:));     
    
    lSCC_poststim{i} = (lSCC_data{1,i}.S_poststim(:,keep_stim_idx,:));  
    rSCC_poststim{i} = (rSCC_data{1,i}.S_poststim(:,keep_stim_idx,:));       
end 

disp('got spectral data from files')

metadata.preprocessing = lSCC_data{1, 1}.metadata.preprocessing  ; 
metadata.epoched = lSCC_data{1, 1}.metadata.epoched;
metadata.PSDparams = lSCC_data{1, 1}.metadata.PSDparams;

metadata.channelinfo.good_ch_labels = good_ch_labels;
metadata.channelinfo.summary_ch_info = summary_ch_info; 
%% Index baseline data 

bl_data.S  = bl_data.S(:,keep_bl_idx,:); 
%% concatenate the right and left SCC stim data.

%prestim 
all_SCC_prestim = cat(1,lSCC_prestim{:},rSCC_prestim{:}); 
% stim 1
all_SCC_stim1  = cat(1,lSCC_stim1{:},rSCC_stim1{:}); 
% stim 2
all_SCC_stim2 = cat(1,lSCC_stim2{:},rSCC_stim2{:}); 
% stim 3 
all_SCC_stim3 = cat(1,lSCC_stim3{:},rSCC_stim3{:}); 
% poststim 
all_SCC_poststim = cat(1,lSCC_poststim{:},rSCC_poststim{:}); 

assert(size(all_SCC_prestim,1)==size(all_SCC_stim1,1)); 
assert(size(all_SCC_prestim,1)==size(all_SCC_poststim,1) & size(all_SCC_stim2,1)==size(all_SCC_stim3,1)); 

disp('num of trials equal across the board for SCC') 
%% Concatenate all labels describing trials and SCC data. 

% concatenate dBS labels 
lSCC_prestim_all = cat(1,lSCC_prestim{:}); 
DBS_lSCC_prestim = repmat({'lSCC'},size(lSCC_prestim_all,1),1); %repeat dbs target label based on num trials 
rSCC_prestim_all = cat(1,rSCC_prestim{:}); 
DBS_rSCC_prestim = repmat({'rSCC'},size(rSCC_prestim_all,1),1); 

%concatenate 5 times because DBS labels are same for all stim states. 
DBS_SCC_labels = vertcat(DBS_lSCC_prestim,DBS_rSCC_prestim,DBS_lSCC_prestim,DBS_rSCC_prestim,...
    DBS_lSCC_prestim,DBS_rSCC_prestim,DBS_lSCC_prestim,DBS_rSCC_prestim,...
    DBS_lSCC_prestim,DBS_rSCC_prestim); 

% concatenate stim labels 
stim_state_SCC_prestim = repmat({'prestim'},size(all_SCC_prestim,1),1); 
stim_state_SCC_stim1 = repmat({'stim1'},size(all_SCC_stim1,1),1); 
stim_state_SCC_stim2 = repmat({'stim2'},size(all_SCC_stim2,1),1); 
stim_state_SCC_stim3 = repmat({'stim3'},size(all_SCC_stim3,1),1); 
stim_state_SCC_poststim = repmat({'poststim'},size(all_SCC_poststim,1),1); 
disp('stim state labels combined')

% combine all data 
all_SCC_data = cat(1,all_SCC_prestim,all_SCC_stim1,all_SCC_stim2,all_SCC_stim3,...
    all_SCC_poststim); 
disp('combined all SCC data') 

% combine all stim labels 
stim_state_SCC_labels = vertcat(stim_state_SCC_prestim,stim_state_SCC_stim1,...
    stim_state_SCC_stim2,stim_state_SCC_stim3,stim_state_SCC_poststim); 
disp('combine stim state labels across all SCC data') 
%% Get electrode labels for SCC (a bit more complicated) 

elec234_lSCC = repmat({'elec234'},size(lSCC_prestim{1, 1},1),1); 
elec567_lSCC = repmat({'elec567'},size(lSCC_prestim{1, 2},1),1);
elec234_rSCC = repmat({'elec234'},size(rSCC_prestim{1, 1},1),1); 
elec567_rSCC = repmat({'elec567'},size(rSCC_prestim{1, 2},1),1);

% combine all electrode labels into one. 

SCC_eleclabels = vertcat(elec234_lSCC, elec567_lSCC,...
    elec234_rSCC, elec567_rSCC); 

all_SCC_eleclabels = vertcat(SCC_eleclabels,SCC_eleclabels,SCC_eleclabels,...
    SCC_eleclabels,SCC_eleclabels); 

disp('extracted electrode labels and combined them all') 

 %% repeat for VCVS  

for i = 1:numel(lVCVS_files) 
    
    lVCVS_prestim{i} = (lVCVS_data{1,i}.S_prestim(:,keep_stim_idx,:));  
    rVCVS_prestim{i} = (rVCVS_data{1,i}.S_prestim(:,keep_stim_idx,:));     
    
    lVCVS_stim1{i} = (lVCVS_data{1,i}.S_stim1(:,keep_stim_idx,:));  
    rVCVS_stim1{i} = (rVCVS_data{1,i}.S_stim1(:,keep_stim_idx,:));  

    lVCVS_stim2{i} = (lVCVS_data{1,i}.S_stim2(:,keep_stim_idx,:));  
    rVCVS_stim2{i} = (rVCVS_data{1,i}.S_stim2(:,keep_stim_idx,:));  

    lVCVS_stim3{i} = (lVCVS_data{1,i}.S_stim3(:,keep_stim_idx,:));  
    rVCVS_stim3{i} = (rVCVS_data{1,i}.S_stim3(:,keep_stim_idx,:));  

    lVCVS_poststim{i} = (lVCVS_data{1,i}.S_poststim(:,keep_stim_idx,:));  
    rVCVS_poststim{i} = (rVCVS_data{1,i}.S_poststim(:,keep_stim_idx,:)); 
end 

disp('combined all VCVS stim states per hemisphere')

all_VCVS_prestim = cat(1,lVCVS_prestim{:},rVCVS_prestim{:}); 
all_VCVS_stim1= cat(1,lVCVS_stim1{:},rVCVS_stim1{:}); 
all_VCVS_stim2= cat(1,lVCVS_stim2{:},rVCVS_stim2{:}); 
all_VCVS_stim3= cat(1,lVCVS_stim3{:},rVCVS_stim3{:}); 
all_VCVS_poststim = cat(1,lVCVS_poststim{:},rVCVS_poststim{:}); 
disp('cat both VCVS hemis data')

assert(size(all_VCVS_prestim,1)==size(all_VCVS_stim1,1)); 
assert(size(all_VCVS_prestim,1)==size(all_VCVS_poststim,1) & size(all_VCVS_stim2,1)==size(all_VCVS_stim3,1)); 

% concatenate dBS labels 
lVCVS_prestim_all = cat(1,lVCVS_prestim{:}); 
DBS_lVCVS_prestim = repmat({'lVCVS'},size(lVCVS_prestim_all,1),1); 
rVCVS_prestim_all = cat(1,rVCVS_prestim{:}); 
DBS_rVCVS_prestim = repmat({'rVCVS'},size(rVCVS_prestim_all,1),1); 

%concatenate 5 times because DBS labels are same for all stim states. 
DBS_VCVS_labels = vertcat(DBS_lVCVS_prestim,DBS_rVCVS_prestim,DBS_lVCVS_prestim,DBS_rVCVS_prestim,...
    DBS_lVCVS_prestim,DBS_rVCVS_prestim,DBS_lVCVS_prestim,DBS_rVCVS_prestim,...
    DBS_lVCVS_prestim,DBS_rVCVS_prestim); 

stim_state_VCVS_prestim = repmat({'prestim'},size(all_VCVS_prestim,1),1); 
stim_state_VCVS_stim1 = repmat({'stim1'},size(all_VCVS_stim1,1),1);
stim_state_VCVS_stim2 = repmat({'stim2'},size(all_VCVS_stim2,1),1);
stim_state_VCVS_stim3 = repmat({'stim3'},size(all_VCVS_stim3,1),1);
stim_state_VCVS_poststim = repmat({'poststim'},size(all_VCVS_poststim,1),1);

stim_state_VCVS_labels = vertcat(stim_state_VCVS_prestim,stim_state_VCVS_stim1,...
    stim_state_VCVS_stim2,stim_state_VCVS_stim3,stim_state_VCVS_poststim); 

% combine all data 
all_VCVS_data = cat(1,all_VCVS_prestim,all_VCVS_stim1,all_VCVS_stim2,all_VCVS_stim3,...
    all_VCVS_poststim); 
disp('cat all VCVS data') 

%% Make VCVS electrode labels 


elec234_lVCVS = repmat({'elec234'},size(lVCVS_prestim{1, 1},1),1);
elec567_lVCVS = repmat({'elec567'},size(lVCVS_prestim{1, 2},1),1);

elec234_rVCVS = repmat({'elec234'},size(rVCVS_prestim{1, 1},1),1);
elec567_rVCVS = repmat({'elec567'},size(rVCVS_prestim{1, 2},1),1);


VCVS_eleclabels = vertcat(elec234_lVCVS,elec567_lVCVS,elec234_rVCVS,elec567_rVCVS);

all_VCVS_eleclabels = vertcat(VCVS_eleclabels,VCVS_eleclabels,VCVS_eleclabels,...
    VCVS_eleclabels,VCVS_eleclabels); 

disp('generated and cat elec labels')

%% Concatenate all DBS labels and all stim labels 

all_stim_labels = vertcat(stim_state_SCC_labels,stim_state_VCVS_labels);
all_DBS_labels = vertcat(DBS_SCC_labels,DBS_VCVS_labels); 
all_elec_labels = vertcat(all_SCC_eleclabels,all_VCVS_eleclabels); 

assert(length(all_stim_labels)==length(all_DBS_labels) && length(all_stim_labels)==length(all_elec_labels)); 
disp('woo everything matches up') 
%% Concatenate all data 

all_data = vertcat(all_SCC_data,all_VCVS_data);

% assert
assert(length(all_stim_labels)==size(all_data,1)); 
disp('size of all data')
disp(size(all_data))

%% Split into freqs 

f = lSCC_data{1, 2}.f_poststim; 

freqband = {'delta', 'theta','alpha','beta','lowgamma','highgamma'}; 
freq_cond = freqband;
avg_FOI = {}; 

for j = 1:numel(freqband)
    freqname = freqband{j}; 
    switch freqname
        case 'delta'
        lower_freq = 1; 
        upper_freq = 4;
        case 'theta'
        lower_freq = 4; 
        upper_freq = 8; % Prior to 05/20/22 it was 7 
        case 'alpha'
        lower_freq = 8; 
        upper_freq = 12;
        case 'beta'
        lower_freq = 12; 
        upper_freq = 30; 
        case 'lowgamma'
        lower_freq = 30; 
        upper_freq = 55; 
        case 'highgamma'
        lower_freq = 65; 
        upper_freq = 100; 
    end
    
freq_idx = find(f >lower_freq & f < upper_freq); 
freq_f{j} = f(freq_idx);

avgFOI = @(pow,FOI_idx)(mean(pow(:,:,FOI_idx),3)); 

avg_FOI{j} = avgFOI(all_data,freq_idx); 

end 
disp('extracted FOI')

%% concatenate into one matrix 

all_data_all_f = squeeze(cat(4,avg_FOI{:}));
disp('cat FOI data into matrix') 
disp('size of all data all f matrix')
disp(size(all_data_all_f))
%% append baseline data

bl_S_data = bl_data.S  ; 
% get frequencies that matter: 
bl_FOI = {}; 
f = bl_data.f  ; 
for j = 1:numel(freqband)
    freqname = freqband{j}; 
    switch freqname
        case 'delta'
        lower_freq = 1; 
        upper_freq = 4;
        case 'theta'
        lower_freq = 4; 
        upper_freq = 8; % Prior to 05/20/22 it was 7 
        case 'alpha'
        lower_freq = 8; 
        upper_freq = 12;
        case 'beta'
        lower_freq = 12; 
        upper_freq = 30; 
        case 'lowgamma'
        lower_freq = 30; 
        upper_freq = 55; 
        case 'highgamma'
        lower_freq = 65; 
        upper_freq = 100; 
    end
    
freq_idx = find(f >lower_freq & f < upper_freq); 
freq_f{j} = f(freq_idx);

avgFOI = @(pow,bl_FOI_idx)(mean(pow(:,:,bl_FOI_idx),3)); 

bl_FOI{j} = avgFOI(bl_S_data,freq_idx); 

end 

bl_data_all_f = squeeze(cat(4,bl_FOI{:}));


%% Concatenate bl and stim labels 
%stim state label 
bl_label = repmat({'BaselineFix'},size(bl_S_data,1),1); 


stim_labels = vertcat(all_stim_labels,bl_label); 
elec_labels = vertcat(all_elec_labels,bl_label);
DBS_labels = vertcat(all_DBS_labels,bl_label); 

assert(length(stim_labels)==length(elec_labels) && length(DBS_labels)==length(stim_labels))


%% FINALLY concatenate all data 

data = vertcat(all_data_all_f,bl_data_all_f); 

%% make tbl with all conds in case we need it. 

tbl = array2table(stim_labels); 
tbl.(2) = DBS_labels;
tbl.(3) = elec_labels; 
tbl.(4) = repmat('f130',size(data,1),1);  

tbl.Properties.VariableNames = {'Stim_Labels','DBS_target','Elec_label','Frequency'}; 

%% Save 
the_date = date; 
addpath(genpath('/Users/anushaallawala/Data/')); 
filename = sprintf('/Users/anushaallawala/Data/DBSTRD/%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID) ; 

save (filename,'data','freq_cond','freq_f','stim_labels','DBS_labels','elec_labels', 'tbl','the_date') 

















