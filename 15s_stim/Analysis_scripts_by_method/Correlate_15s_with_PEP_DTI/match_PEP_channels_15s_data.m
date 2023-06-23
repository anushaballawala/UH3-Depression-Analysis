PatientID = 'DBSTRD001';
hemi ='l';
DBS_target ='SCC'; stim_freq = 130; 
contact_config = 'elec1';
PEP_contact_config = '1';
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 

% load my channel info 

data = load(sprintf('%s%s',experiment_name,'theta')); 
Ch_info_15s = data.metadata.preprocessing.GoodChannelLabels ; 

%load pEP matrix 
PEP_matrix = load(sprintf('E:/DBSTRD/%s/Experiments/PEP_data/Threshold/Threshold_PEP-%s%s-%s.mat',...
    PatientID,hemi,DBS_target,PEP_contact_config)); 

%load Josh Montageinfo file 

load(sprintf('E:/DBSTRD/%s/Experiments/PEP_data/MontageInfo.mat',PatientID)); 

%% Align/index my channels with Josh's channels 






