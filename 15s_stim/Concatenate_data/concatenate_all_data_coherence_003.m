%% Script: Concatenate_all_data_coherence_003.m

%Function: this script takes all the spectral analysis outputs from the 
% Chronux mtspectrum.c function (saved on the Oscar server) and
% concatenates all of the data from SCC, VCVS stim, plus baseline
% recordings (5 min pre experiment). It also concatenates all of the labels
% for stim state, DBS target, current configuration. This has only been
% written for 130 Hz stim experiments. 

%  Outputs: 
% 1) single 3D matrix that is trials x trials x channels x frequencies. 
%  This output is used for stats and viz analysis. 
% 2) electrode labels, stim_state (prestim, poststim, etc.), DBS target
% label for each trial. 

% Note: I RUN THIS SCRIPT LOCALLY, NOT ON OSCAR. 

% Warning: script is long. Plz be v careful about making any changes. 
% AA 02/07/22

clear 
clc  
% add paths 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD003/Experiments/'));


% ------------------------ START CODE -------------------------------% 
%% Load stim data. 

PatientID = 'DBSTRD003'; 

lSCC_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*lSCC*.mat',PatientID));
rSCC_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*rSCC*.mat',PatientID));
lVCVS_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*lVCVS*.mat',PatientID));
rVCVS_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*rVCVS*.mat',PatientID));

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

bl_dir = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/BaselineFix_date-04-18-2021/Coherence/*.mat',PatientID));
bl_data = load(sprintf('%s/%s',bl_dir.folder,bl_dir.name));


%% take out channels that are in baselinefix and not in 15s stim 

%load ROI labels 

ROI_labels = load('/Users/anushaallawala/Data/DBSTRD003_ROI_labels_new.mat');

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

summary_ch_info = ROI_labels.tbl ; 
%remove channels from summary table 
%summary_ch_info(stim_remove_idx,:) = []; 
%% get coherence data 

for i = 1:numel(lSCC_files) 
    
    lSCC_delta{i} = lSCC_data{1, 1}.delta_indiv_tr  ; 
    rSCC_delta{i} = rSCC_data{1, 1}.delta_indiv_tr  ; 
    
    lSCC_theta{i} = lSCC_data{1, 1}.theta_indiv_tr  ; 
    rSCC_theta{i} = rSCC_data{1, 1}.theta_indiv_tr  ; 
    
    lSCC_alpha{i} = lSCC_data{1, 1}.alpha_indiv_tr  ; 
    rSCC_alpha{i} = rSCC_data{1, 1}.alpha_indiv_tr ; 
    
    lSCC_beta{i} =lSCC_data{1, 1}.beta_indiv_tr  ; 
    rSCC_beta{i}= rSCC_data{1, 1}.beta_indiv_tr  ; 
    
    lSCC_lowgamma{i} = lSCC_data{1, 1}.lowgamma_indiv_tr  ; 
    rSCC_lowgamma{i}= rSCC_data{1, 1}.lowgamma_indiv_tr  ; 
    
    lSCC_highgamma{i} = lSCC_data{1, 1}.highgamma_indiv_tr  ; 
    rSCC_highgamma{i} = rSCC_data{1, 1}.highgamma_indiv_tr  ; 
    
    
end 

disp('got coherence data from files')

metadata.channelinfo.good_ch_labels = good_ch_labels;
metadata.channelinfo.summary_ch_info = summary_ch_info; 
metadata.files.lSCC_files = lSCC_files;
metadata.files.rSCC_files = rSCC_files; 
metadata.files.lVCVS_files = lVCVS_files; 
metadata.files.rVCVS_files = rVCVS_files; 
%% Get baseline data 

bl_data_delta = bl_data.delta_indiv_tr;
bl_data_theta= bl_data.theta_indiv_tr;
bl_data_alpha= bl_data.alpha_indiv_tr;
bl_data_beta= bl_data.beta_indiv_tr;
bl_data_lowgamma= bl_data.lowgamma_indiv_tr;
bl_data_highgamma= bl_data.highgamma_indiv_tr;

%% Combine left SCC data 

lSCC_all_delta = cat(1,lSCC_delta{:});
lSCC_all_theta = cat(1,lSCC_theta{:}); 
lSCC_all_alpha = cat(1,lSCC_alpha{:}); lSCC_all_beta = cat(1,lSCC_beta{:}); 
lSCC_all_lowgamma = cat(1,lSCC_lowgamma{:}); lSCC_all_highgamma = cat(1,lSCC_highgamma{:}); 
%% Combine right SCC data 


rSCC_all_delta = cat(1,rSCC_delta{:});
rSCC_all_theta = cat(1,rSCC_theta{:}); 
rSCC_all_alpha = cat(1,rSCC_alpha{:}); rSCC_all_beta = cat(1,rSCC_beta{:}); 
rSCC_all_lowgamma = cat(1,rSCC_lowgamma{:}); rSCC_all_highgamma = cat(1,rSCC_highgamma{:}); 

%% Combine both SCC data 

SCC_all_delta = vertcat(lSCC_delta{:},rSCC_delta{:}); 
SCC_all_theta = vertcat(lSCC_theta{:},rSCC_theta{:}); 
SCC_all_alpha = vertcat(lSCC_alpha{:}, rSCC_alpha{:}); 
SCC_all_beta = vertcat(lSCC_beta{:}, rSCC_beta{:});
SCC_all_lowgamma = vertcat(lSCC_lowgamma{:}, rSCC_lowgamma{:});
SCC_all_highgamma = vertcat(lSCC_highgamma{:}, rSCC_highgamma{:}); 

%% label SCC data for DBS and frequency band 

lSCC_label = repmat({'lSCC'}, size(lSCC_all_delta,1),1);  
rSCC_label = repmat({'rSCC'}, size(rSCC_all_delta,1),1); 

all_SCC_label = vertcat(lSCC_label,rSCC_label); 
%% Get VCVS data 

for i = 1:numel(lVCVS_files) 
    
    lVCVS_delta{i} = lVCVS_data{1, 1}.delta_indiv_tr  ; 
    rVCVS_delta{i} = rVCVS_data{1, 1}.delta_indiv_tr  ; 
    
    lVCVS_theta{i} = lVCVS_data{1, 1}.theta_indiv_tr  ; 
    rVCVS_theta{i} = rVCVS_data{1, 1}.theta_indiv_tr  ; 
    
    lVCVS_alpha{i} = lVCVS_data{1, 1}.alpha_indiv_tr  ; 
    rVCVS_alpha{i} = rVCVS_data{1, 1}.alpha_indiv_tr ; 
    
    lVCVS_beta{i} =lVCVS_data{1, 1}.beta_indiv_tr  ; 
    rVCVS_beta{i}= rVCVS_data{1, 1}.beta_indiv_tr  ; 
    
    lVCVS_lowgamma{i} = lVCVS_data{1, 1}.lowgamma_indiv_tr  ; 
    rVCVS_lowgamma{i}= rVCVS_data{1, 1}.lowgamma_indiv_tr  ; 
    
    lVCVS_highgamma{i} = lVCVS_data{1, 1}.highgamma_indiv_tr  ; 
    rVCVS_highgamma{i} = rVCVS_data{1, 1}.highgamma_indiv_tr  ; 
    
    
end 

%% Combine left VCVS data 

lVCVS_all_delta = cat(1,lVCVS_delta{:});
lVCVS_all_theta = cat(1,lVCVS_theta{:}); 
lVCVS_all_alpha = cat(1,lVCVS_alpha{:}); lVCVS_all_beta = cat(1,lVCVS_beta{:}); 
lVCVS_all_lowgamma = cat(1,lVCVS_lowgamma{:}); lVCVS_all_highgamma = cat(1,lVCVS_highgamma{:}); 
%% Combine right VCVS data 


rVCVS_all_delta = cat(1,rVCVS_delta{:});
rVCVS_all_theta = cat(1,rVCVS_theta{:}); 
rVCVS_all_alpha = cat(1,rVCVS_alpha{:}); rVCVS_all_beta = cat(1,rVCVS_beta{:}); 
rVCVS_all_lowgamma = cat(1,rVCVS_lowgamma{:}); rVCVS_all_highgamma = cat(1,rVCVS_highgamma{:}); 

%% Combine both hemispheres for VCVS 

VCVS_all_delta = vertcat(lVCVS_delta{:},rVCVS_delta{:}); 
VCVS_all_theta = vertcat(lVCVS_theta{:},rVCVS_theta{:}); 
VCVS_all_alpha = vertcat(lVCVS_alpha{:}, rVCVS_alpha{:}); 
VCVS_all_beta = vertcat(lVCVS_beta{:}, rVCVS_beta{:});
VCVS_all_lowgamma = vertcat(lVCVS_lowgamma{:}, rVCVS_lowgamma{:});
VCVS_all_highgamma = vertcat(lVCVS_highgamma{:}, rVCVS_highgamma{:}); 

%% label VCVS data for DBS and frequency band 

lVCVS_label = repmat({'lVCVS'}, size(lVCVS_all_delta,1),1);  
rVCVS_label = repmat({'rVCVS'}, size(rVCVS_all_delta,1),1); 

all_VCVS_label = vertcat(lVCVS_label,rVCVS_label); 


%% Combine all stim conditions & labels

all_SCC_VCVS_label = vertcat(all_SCC_label, all_VCVS_label); 


all_SCC_VCVS_delta = vertcat(SCC_all_delta, VCVS_all_delta); 
all_SCC_VCVS_theta = vertcat(SCC_all_theta, VCVS_all_theta); 
all_SCC_VCVS_alpha = vertcat(SCC_all_alpha, VCVS_all_alpha);
all_SCC_VCVS_beta = vertcat(SCC_all_beta, VCVS_all_beta);
all_SCC_VCVS_lowgamma = vertcat(SCC_all_lowgamma, VCVS_all_lowgamma); 
all_SCC_VCVS_highgamma = vertcat(SCC_all_highgamma, VCVS_all_highgamma); 

%% Get electrode labels for SCC
elec234_lSCC = repmat({'elec234'},5,1); 
elec567_lSCC = repmat({'elec567'},5,1); 
elec234_rSCC = repmat({'elec234'},5,1); 
elec567_rSCC = repmat({'elec567'},5,1); 

% combine all electrode labels into one. 

SCC_eleclabels = vertcat(elec234_lSCC, elec567_lSCC,...
    elec234_rSCC, elec567_rSCC); 

disp('extracted electrode labels for SCC and combined them all') 

%% Get electrode labels for VCVS 
elec234_lVCVS = repmat({'elec234'},5,1); 
elec567_lVCVS = repmat({'elec567'},5,1); 
elec234_rVCVS = repmat({'elec234'},5,1); 
elec567_rVCVS = repmat({'elec567'},5,1); 

% combine all electrode labels into one. 

VCVS_eleclabels = vertcat(elec234_lVCVS, elec567_lVCVS,...
    elec234_rVCVS, elec567_rVCVS); 

disp('extracted electrode labels for VCVS and combined them all') 


%% combine electrode labels 

all_elec_labels = vertcat(SCC_eleclabels,VCVS_eleclabels); 

disp('electrode config labels combined for both SCC and VCVS') 
%% Stim state label 

poststim_labels = repmat({'poststim'}, size(all_SCC_VCVS_label,1),1);

disp('got all poststim labels') 

%% Append baseline data & labels 

% append data 

all_delta = vertcat(all_SCC_VCVS_delta(:,keep_stim_idx,keep_stim_idx),bl_data_delta(:,keep_bl_idx,keep_bl_idx)); 
all_theta = vertcat(all_SCC_VCVS_theta(:,keep_stim_idx,keep_stim_idx),bl_data_theta(:,keep_bl_idx,keep_bl_idx));
all_alpha = vertcat(all_SCC_VCVS_alpha(:,keep_stim_idx,keep_stim_idx),bl_data_alpha(:,keep_bl_idx,keep_bl_idx));
all_beta = vertcat(all_SCC_VCVS_beta(:,keep_stim_idx,keep_stim_idx),bl_data_beta(:,keep_bl_idx,keep_bl_idx));
all_lowgamma = vertcat(all_SCC_VCVS_lowgamma(:,keep_stim_idx,keep_stim_idx),bl_data_lowgamma(:,keep_bl_idx,keep_bl_idx));
all_highgamma = vertcat(all_SCC_VCVS_highgamma(:,keep_stim_idx,keep_stim_idx),bl_data_highgamma(:,keep_bl_idx,keep_bl_idx));

% append labels 
baseline_state_labels = repmat({'BaselineFix'},size(bl_data_alpha,1),1); 

%append to stim state 
stim_labels = vertcat(poststim_labels, baseline_state_labels); 
% append elec conf labels 
elec_labels = vertcat(all_elec_labels, baseline_state_labels); 
%append dbs target 
DBS_labels = vertcat(all_SCC_VCVS_label,baseline_state_labels); 

%% concatenate all frequency data 

all_data_all_f = cat(4,all_delta,all_theta,all_alpha,all_beta,all_lowgamma,all_highgamma);

disp('cat FOI data into matrix') 
disp('size of all data all f matrix')
disp(size(all_data_all_f))

%FOI labels 
freq_cond = {'delta', 'theta','alpha','beta','lowgamma','highgamma'}; 

data = all_data_all_f; 
%% make tbl with all conds in case we need it. 

tbl = array2table(stim_labels); 
tbl.(2) = DBS_labels;
tbl.(3) = elec_labels; 
tbl.(4) = repmat('f130',size(data,1),1);  

tbl.Properties.VariableNames = {'Stim_Labels','DBS_target','Elec_label','Frequency'}; 

%% Save 
the_date = date; 
addpath(genpath('/Users/anushaallawala/Data/')); 
filename = sprintf('/Users/anushaallawala/Data/DBSTRD/%s_autocorrbl_data_for_coherence_stats.mat',PatientID) ; 

save (filename,'data','freq_cond','stim_labels','DBS_labels','elec_labels', 'tbl','the_date','metadata') 





