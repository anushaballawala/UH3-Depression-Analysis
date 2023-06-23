%%%% TO-DO:
% 1. add integrity checks like PCCT 
% 2. add filepaths/dependencies 
% 3. add filepath depending on computer 

%% %% get options for processing
%%% PURPOSE: generate epoched trials for each condition for 15-s stim data
%%% for each file %%% 
%%% INPUT: epoch file with timestamps 
%%% INPUT: corresponding neuraldata file that has been preprocessed 
%%% OUTPUT: epoched data, each cell containing trial x ch x time 
%%% created by AA for UH3 Depression project 06/28/21 %%% 

function [output_and_metadata] = conditionbin_time_series_15s_per_file(PatientID,metadata,neuraldata,epochdata,stim_info) 
    

%% set up epoching variables 

srate = metadata.preprocessing.New_SamplingRate;%assuming all files have the same srate 

trial_len_secs = 25; %trial length in seconds ** will change for other pts 
stim_win = 15; % stim window in seconds 
trial_length = trial_len_secs*srate; 
num_channels = size(neuraldata,1); 
length_time = size(neuraldata,2);
tbl = epochdata.tbl; 

%load data info 
[output_and_metadata,current_condition,conditions,tbl_idx,...
    num_conditions] = bin_time_series_data(neuraldata,tbl,srate,length_time,trial_length,stim_win,num_channels);  

%% Assign metadata 
metadata.epoched.num_conditions = num_conditions; 
metadata.epoched.current_condition = current_condition; 
metadata.epoched.table_summary_conditions = tbl_idx; 
metadata.epoched.trial_averaged_conditions = conditions;
metadata.epoched.trial_length = trial_length; 
metadata.epoched.trial_length_secs = trial_len_secs; 
metadata.epoched.total_file_length = length_time; 
metadata.epoched.epoch_table = tbl; 

%% Create dir 
 
outputdir = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/%s/EXP/15s_stim/Epoched Data/%s',PatientID,stim_info);
if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end
%% save indiv trial data all in one file **
thisfile = sprintf('15s_stim_all_currdir_timeseries_singletrial_%s.mat', stim_info);
save(sprintf('%s/%s',outputdir,thisfile),'output_and_metadata','metadata','-v7.3'); 

end 




