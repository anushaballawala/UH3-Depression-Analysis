%% PREPROCESSING SETUP 
clear all 
%setup files for preprocessing 
subjecttype = 'DBSTRD'; 
PatientID = 'DBSTRD002'; 
experiment_type = '15s_stim'; 
freq = 130; 
DBStarget = 'SCC'; 
hemi = 'l'; 
run = '01' ; 
remove_artifact = 'yes'; % enter 'yes' if using PARRM, 'no' if not 

% ***** Check that correct Patient ID number is in setup script *****
[outputdir,montageInfo,metadata,data,good_ch_idx,good_ch_labels] = Preprocessing_setup_15s_stim_DBSTRD001(subjecttype,PatientID,experiment_type,hemi,DBStarget,run,freq); 

data = data(montageInfo.MacroContactIndices,:); 
%% 
%behav_type = metadata.general.Experiment; 


if strcmp(remove_artifact,'yes') ==  1 % use PARRM to remove artifact 
    epochdatafile = sprintf('E:/DBSTRD/%s/Experiments/%s/Epochs/epoch_%s_%s%s_f%d.mat',...
    PatientID,experiment_type,PatientID,hemi,DBStarget,freq);
    epochdata = load(epochdatafile);
    preprocessing_15s_PARRM(montageInfo,data,metadata,outputdir,good_ch_idx',good_ch_labels,epochdata)
elseif strcmp(remove_artifact, 'no') == 1
    preprocessing_15s_calculations(montageInfo,data,metadata,outputdir,good_ch_idx',good_ch_labels)
end 





