%% Script: Concatenate_all_data_coherence.m

%Function: this script takes all the spectral analysis outputs from the 
% Chronux mtspectrum.c function (saved on the Oscar server) and
% concatenates all of the data from SCC, VCVS stim, plus baseline
% recordings (5 min pre experiment). It also concatenates all of the labels
% for stim state, DBS target, current configuration. This has only been
% written for 130 Hz stim experiments. 

%  Outputs: 
% 1) single 4D matrix that is trials x channels x channels x frequencies. 
%  This output is used for stats and viz analysis. 
% 2) electrode labels, stim_state (prestim, poststim, etc.), DBS target
% label for each trial. 

% Note: I RUN THIS SCRIPT LOCALLY, NOT ON OSCAR. 

% Warning: script is long. Plz be v careful about making any changes. 
% AA 02/07/22

clear 
clc  
% add paths 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/Coherence/')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD001/Experiments/BaselineFix/Processed Data/Coherence/')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD002/Experiments/15s_stim/Processed Data/Coherence/')); 
addpath(genpath('/Volumes/TRD_Project/DBSTRD/DBSTRD002/Experiments/BaselineFix/Processed Data/Coherence/')); 


% ------------------------ START CODE -------------------------------% 
%% Load stim data. 

PatientID = 'DBSTRD002'; 

lSCC_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*lSCC*_alltrials.mat',PatientID));
rSCC_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*rSCC*_alltrials.mat',PatientID));
lVCVS_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*lVCVS*_alltrials.mat',PatientID));
rVCVS_files = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/*rVCVS*_alltrials.mat',PatientID));

%% 

    lSCC_data = load(sprintf('%s/%s',lSCC_files.folder, lSCC_files.name)); 
    rSCC_data = load(sprintf('%s/%s',rSCC_files.folder, rSCC_files.name)); 
    lSCC_name = lSCC_files.name;
    rSCC_name = rSCC_files.name;

    lVCVS_data = load(sprintf('%s/%s',lVCVS_files.folder, lVCVS_files.name)); 
    lVCVS_name = lVCVS_files.name;
    rVCVS_data = load(sprintf('%s/%s',rVCVS_files.folder, rVCVS_files.name));  
    rVCVS_name = rVCVS_files.name;


%% Load baseline data. 

switch PatientID 
    case 'DBSTRD001' 
        bl_dir = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/Coherence/*alltrials.mat',PatientID));       
        bl_data = load(sprintf('%s/%s',bl_dir.folder,bl_dir.name)); 
    
    case 'DBSTRD002' 
        bl_dir = dir(sprintf('/Volumes/TRD_Project/DBSTRD/%s/Experiments/BaselineFix/Processed Data/Coherence/*alltrials.mat',PatientID)); 
        bl_data = load(sprintf('%s/%s',bl_dir.folder,bl_dir.name)); 
       
end 

%% get coherence data 

    
    lSCC_delta = lSCC_data.delta_indiv_tr  ; 
    rSCC_delta= rSCC_data.delta_indiv_tr  ; 
    
    lSCC_theta = lSCC_data.theta_indiv_tr  ; 
    rSCC_theta = rSCC_data.theta_indiv_tr  ; 
    
    lSCC_alpha = lSCC_data.alpha_indiv_tr  ; 
    rSCC_alpha = rSCC_data.alpha_indiv_tr ; 
    
    lSCC_beta = lSCC_data.beta_indiv_tr  ; 
    rSCC_beta = rSCC_data.beta_indiv_tr  ; 
    
    lSCC_lowgamma = lSCC_data.lowgamma_indiv_tr  ; 
    rSCC_lowgamma = rSCC_data.lowgamma_indiv_tr  ; 
    
    lSCC_highgamma = lSCC_data.highgamma_indiv_tr  ; 
    rSCC_highgamma = rSCC_data.highgamma_indiv_tr  ; 
    

disp('got SCC coherence data from files')


%% Combine both SCC data 

SCC_all_delta = vertcat(lSCC_delta,rSCC_delta); 
SCC_all_theta = vertcat(lSCC_theta,rSCC_theta); 
SCC_all_alpha = vertcat(lSCC_alpha, rSCC_alpha); 
SCC_all_beta = vertcat(lSCC_beta, rSCC_beta);
SCC_all_lowgamma = vertcat(lSCC_lowgamma, rSCC_lowgamma);
SCC_all_highgamma = vertcat(lSCC_highgamma, rSCC_highgamma); 

%% label SCC data for DBS and frequency band 

lSCC_label = repmat({'lSCC'}, size(lSCC_delta,1),1);  
rSCC_label = repmat({'rSCC'}, size(rSCC_delta,1),1); 

all_SCC_label = vertcat(lSCC_label,rSCC_label); 
%% Get VCVS data 
lVCVS_delta = lVCVS_data.delta_indiv_tr  ;
rVCVS_delta= rVCVS_data.delta_indiv_tr  ;

lVCVS_theta = lVCVS_data.theta_indiv_tr  ;
rVCVS_theta = rVCVS_data.theta_indiv_tr  ;

lVCVS_alpha = lVCVS_data.alpha_indiv_tr  ;
rVCVS_alpha = rVCVS_data.alpha_indiv_tr ;

lVCVS_beta =lVCVS_data.beta_indiv_tr  ;
rVCVS_beta= rVCVS_data.beta_indiv_tr  ;

lVCVS_lowgamma = lVCVS_data.lowgamma_indiv_tr  ;
rVCVS_lowgamma= rVCVS_data.lowgamma_indiv_tr  ;

lVCVS_highgamma = lVCVS_data.highgamma_indiv_tr  ;
rVCVS_highgamma = rVCVS_data.highgamma_indiv_tr  ;

%% Combine both hemispheres for VCVS 

VCVS_all_delta = vertcat(lVCVS_delta,rVCVS_delta); 
VCVS_all_theta = vertcat(lVCVS_theta,rVCVS_theta); 
VCVS_all_alpha = vertcat(lVCVS_alpha, rVCVS_alpha); 
VCVS_all_beta = vertcat(lVCVS_beta, rVCVS_beta);
VCVS_all_lowgamma = vertcat(lVCVS_lowgamma, rVCVS_lowgamma);
VCVS_all_highgamma = vertcat(lVCVS_highgamma, rVCVS_highgamma); 

%% label VCVS data for DBS and frequency band 

lVCVS_label = repmat({'lVCVS'}, size(lVCVS_delta,1),1);  
rVCVS_label = repmat({'rVCVS'}, size(rVCVS_delta,1),1); 

all_VCVS_label = vertcat(lVCVS_label,rVCVS_label); 
%% Get baseline data 

bl_data_delta = bl_data.delta_indiv_tr;
bl_data_theta= bl_data.theta_indiv_tr;
bl_data_alpha= bl_data.alpha_indiv_tr;
bl_data_beta= bl_data.beta_indiv_tr;
bl_data_lowgamma= bl_data.lowgamma_indiv_tr;
bl_data_highgamma= bl_data.highgamma_indiv_tr;

%% Combine all stim conditions & labels

all_SCC_VCVS_label = vertcat(all_SCC_label, all_VCVS_label); 

all_SCC_VCVS_delta = vertcat(SCC_all_delta, VCVS_all_delta); 
all_SCC_VCVS_theta = vertcat(SCC_all_theta, VCVS_all_theta); 
all_SCC_VCVS_alpha = vertcat(SCC_all_alpha, VCVS_all_alpha);
all_SCC_VCVS_beta = vertcat(SCC_all_beta, VCVS_all_beta);
all_SCC_VCVS_lowgamma = vertcat(SCC_all_lowgamma, VCVS_all_lowgamma); 
all_SCC_VCVS_highgamma = vertcat(SCC_all_highgamma, VCVS_all_highgamma); 

%% Concatenate all electrode configurations 
%%%%%%%%%%* CHECK CODE ON VNC VIEWER FOR WHAT ORDER THE ELECTRODE CONTACT
%%%%%%%%%%IS CONCATENATED IN IN COMPUTE COHERENCE FOR INDIV TRIAL CODE 
%% Get electrode labels for SCC (a bit more complicated) 

switch PatientID 
    case 'DBSTRD001'
        elec1_lSCC = repmat({'elec1'},5,1); 
        elec234_lSCC = repmat({'elec234'},5,1); 
        elec25_lSCC = repmat({'elec25'},5,1); 
        elec36_lSCC = repmat({'elec36'},5,1); 
        elec47_lSCC = repmat({'elec47'},5,1); 
        elec567_lSCC = repmat({'elec567'},5,1); 
        elec8_lSCC = repmat({'elec8'},5,1); 
        %right hemisphere is equal to left hemisphere elec configs for SCC

        elec1_lVCVS = repmat({'elec1'},4,1); 
        elec25_lVCVS = repmat({'elec25'},5,1); 
        elec36_lVCVS = repmat({'elec36'},5,1); 
        elec47_lVCVS = repmat({'elec47'}, 5, 1); 
        elec8_lVCVS = repmat({'elec8'},5,1); 

        elec1_rVCVS = repmat({'elec1'},5,1); 
        elec25_rVCVS = repmat({'elec25'},5,1); 
        elec36_rVCVS = repmat({'elec36'}, 5,1); 
        elec47_rVCVS = repmat({'elec47'},5,1); 
        elec8_rVCVS = repmat({'elec8'},5,1); 
        
        
        %concatenate all labels 
        % left SCC, right SCC, left VCvs, right VCVS 
        
        all_elec_labels = vertcat(elec1_lSCC, elec234_lSCC, elec25_lSCC, elec36_lSCC, elec47_lSCC, elec567_lSCC, elec8_lSCC,...
            elec1_lSCC, elec234_lSCC, elec25_lSCC, elec36_lSCC, elec47_lSCC, elec567_lSCC, elec8_lSCC,...
            elec1_lVCVS,elec25_lVCVS, elec36_lVCVS, elec47_lVCVS, elec8_lVCVS,...
            elec1_rVCVS, elec25_rVCVS, elec36_rVCVS, elec47_rVCVS, elec8_rVCVS) ; 
        
        size(all_elec_labels) 
    case 'DBSTRD002'
        
        elec1_lSCC = repmat({'elec1'},5,1); 
        elec234_lSCC = repmat({'elec234'},5,1); 
        elec25_lSCC = repmat({'elec25'},5,1); 
        elec36_lSCC = repmat({'elec36'},10,1); 
        elec47_lSCC = repmat({'elec47'},5,1); 
        elec567_lSCC = repmat({'elec567'},5,1); 
        elec8_lSCC = repmat({'elec8'},5,1); 
        
        elec1_rSCC = repmat({'elec1'},5,1); 
        elec234_rSCC = repmat({'elec234'},5,1); 
        elec25_rSCC = repmat({'elec25'},10,1); 
        elec36_rSCC = repmat({'elec36'},5,1); 
        elec47_rSCC = repmat({'elec47'},5,1); 
        elec567_rSCC = repmat({'elec567'},5,1); 
        elec8_rSCC = repmat({'elec8'},5,1);
        
        %right hemisphere is equal to left hemisphere elec configs for SCC

        elec1_lVCVS = repmat({'elec1'},5,1); 
        elec234_lVCVS = repmat({'elec234'},5,1);
        elec25_lVCVS = repmat({'elec25'},11,1); 
        elec36_lVCVS = repmat({'elec36'},5,1); 
        elec47_lVCVS = repmat({'elec47'}, 5, 1); 
        elec567_lVCVS = repmat({'elec567'},5,1);
        elec8_lVCVS = repmat({'elec8'},5,1); 

        elec1_rVCVS = repmat({'elec1'},5,1); 
        elec234_rVCVS = repmat({'elec234'},5,1); 
        elec25_rVCVS = repmat({'elec25'},5,1); 
        elec36_rVCVS = repmat({'elec36'}, 5,1); 
        elec47_rVCVS = repmat({'elec47'},5,1); 
        elec567_rVCVS = repmat({'elec567'},5,1);
        elec8_rVCVS = repmat({'elec8'},5,1); 
        
        all_elec_labels = vertcat(elec1_lSCC,elec234_lSCC,elec25_lSCC,elec36_lSCC,elec47_lSCC,elec567_lSCC,elec8_lSCC,...
            elec1_rSCC,elec234_rSCC,elec25_rSCC,elec36_rSCC,elec47_rSCC,elec567_rSCC,elec8_rSCC,...
            elec1_lVCVS,elec234_lVCVS,elec25_lVCVS,elec36_lVCVS,elec47_lVCVS,elec567_lVCVS,elec8_lVCVS,...
            elec1_rVCVS,elec234_rVCVS,elec25_rVCVS,elec36_rVCVS,elec47_rVCVS,elec567_rVCVS,elec8_rVCVS); 
        
end 

assert(length(all_elec_labels) == size(all_SCC_VCVS_alpha,1)); 

%% Load ROI data 


ROI_labels = load(sprintf('/Users/anushaallawala/Data/%s_ROI_labels.mat',PatientID));

metadata.channelinfo.good_ch_labels = ROI_labels.good_ch_labels;
metadata.channelinfo.summary_ch_info = ROI_labels.summary_ch_info; 
metadata.files.lSCC_files = lSCC_files;
metadata.files.rSCC_files = rSCC_files; 
metadata.files.lVCVS_files = lVCVS_files; 
metadata.files.rVCVS_files = rVCVS_files; 

%% Stim state label 

poststim_labels = repmat({'poststim'}, size(all_SCC_VCVS_label,1),1);

disp('got all poststim labels') 

%% get indices for cleaning up channels a bit more 

good_ch_remove = ROI_labels.goodch_idx_remove  ; 

switch PatientID 
    case  'DBSTRD001' 
        keep_stim_idx = [1:58,62:131]; 
        keep_bl_idx = keep_stim_idx; 
    case 'DBSTRD002'
        keep_stim_idx =  [1:57,59:66,68:101,103:129]; 
        keep_bl_idx = keep_stim_idx; 

end 

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

assert(length(stim_labels) == size(all_delta,1))
assert(length(elec_labels) == length(stim_labels))
assert(length(DBS_labels) == length(elec_labels))
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
disp('size of table')
size(tbl)
%% Save 
the_date = date; 
addpath(genpath('/Users/anushaallawala/Data/')); 
filename = sprintf('/Users/anushaallawala/Data/DBSTRD/%s_autocorrbl_data_for_coherence_stats.mat',PatientID) ; 

save (filename,'data','freq_cond','stim_labels','DBS_labels','elec_labels', 'tbl','the_date','metadata') 


