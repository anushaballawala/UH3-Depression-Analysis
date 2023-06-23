%%%% TO-DO:
% 1. add integrity checks like PCCT 
% 2. add filepaths/dependencies 
% 3. add filepath depending on computer 


%% %% get options for processing
%%% PURPOSE: generate epoched trials for each condition for 1-s stim data
%%% for each file %%% 
%%% INPUT: epoch file with timestamps 
%%% INPUT: corresponding neuraldata file that has been preprocessed 
%%% OUTPUT: epoched data, each cell containing chxtrialxtimexfreqs
%%% created by AA for UH3 Depression project 05/12/2020 %%% 

function [] = conditionbin_15s_stim_perfile(subjectID,metadata,neuraldata,epochdata,stim_info) 
    

%% set up epoching variables 

srate = metadata.preprocessing.New_SamplingRate;%assuming all files have the same srate 

trial_len_secs = 30; %trial length in seconds ** 
stim_win = 15; % stim window in seconds 
trial_length = trial_len_secs*srate; 
num_channels = size(neuraldata,2); 

length_time = size(neuraldata(1).amplitude,2);
num_freqs = length(metadata.decomp_parameters.freqs); 
measured_data = zeros(num_channels,num_freqs,length_time); %initialize 

for i = 1:num_channels 
    measured_data(i,:,:) = abs(neuraldata(i).amplitude).^2; %get power 
end 


%load data info 
[output_and_metadata,conditions,condition_bytrial,tbl_idx,num_conditions] = bin_data_15s(measured_data,...
    epochdata.tbl,srate,num_freqs,length_time,trial_length,stim_win,num_channels);

% *add save  line, including saving index for conditions 
%% Assign metadata 
metadata.epoched.num_conditions = num_conditions; 
metadata.epoched.table_summary_conditions = tbl_idx; 
metadata.epoched.trial_averaged_conditions = conditions;
metadata.epoched.trial_length = trial_length; 
metadata.epoched.trial_length_secs = trial_len_secs; 
metadata.epoched.total_file_length = length_time; 

%% Average over trials if needed 

%**if options = re-run everything "yes" then run otherwise dont 

% trial_averaged_pw = cell(num_conditions,1); 
% 
% for cond_idx = 1:num_conditions 
%     trial_averaged_pw{cond_idx} = squeeze(mean(output{cond_idx},1)); 
% end
% 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s',subjectID,stim_info);

if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end
%% save indiv trial data all in one file 
thisfile = sprintf('15s_stim_all_currdir_singletrial_%s.mat', stim_info);
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
save(fulldestination,'output_and_metadata','metadata','-v7.3'); 

%% need to clean everything below this point 
% ** enable saving individual files for each condition if needed in future
% **

% %assign trial averaged conditionst o each variable for faster saving nad
% %loading 
% condition1 = trial_averaged_pw{1}; 
% condition2 = trial_averaged_pw{2}; 
% condition3 = trial_averaged_pw{3}; 
% condition4 = trial_averaged_pw{4};
% condition5 = trial_averaged_pw{5};
% condition6 = trial_averaged_pw{6};
% condition7 = trial_averaged_pw{7};
% 
% %assign indiv trial cnditons for each variable for faster saving and
% %laoding 
% indiv1 = output{1}; 
% indiv2 = output{2}; 
% indiv3  = output{3}; 
% indiv4 = output{4}; 
% indiv5 = output{5}; 
% indiv6 = output{6}; 
% indiv7  = output{7}; 
% 
% %% save and clean up once above is cleaned up 
% 
% %******** CLEAN UP THIS SECTION - INDIV TRIALS WONT SAVE?
% thisfile = sprintf('15s_stim_trialaveraged_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'trial_averaged_pw','metadata','-v7.3'); 
% 
% %cond1 
% thisfile = sprintf('15s_stim_cond1_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition1','indiv1','metadata','-v7.3'); 
% 
% %cond2 
% thisfile = sprintf('15s_stim_cond2_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition2','indiv2','metadata','-v7.3'); 
% 
% %cond3 
% thisfile = sprintf('15s_stim_cond3_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition3','indiv3','metadata','-v7.3'); 
% 
% %cond4 
% thisfile = sprintf('15s_stim_cond4_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition4','indiv4','metadata','-v7.3'); 
% 
% %cond5 
% thisfile = sprintf('15s_stim_cond5_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition5','indiv5','metadata','-v7.3'); 
% 
% %cond6 
% thisfile = sprintf('15s_stim_cond6_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition6','indiv6','metadata','-v7.3'); 
% 
% %cond7 
% thisfile = sprintf('15s_stim_cond7_%s.mat', stim_info);
% fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
% save(fulldestination,'condition7','indiv7','metadata','-v7.3'); 
%
% 
%save (fulldestination,'trial_averaged_pw','epochdatafile','neuraldatafile','freqs','srate','trial_length','cond_idx','epochdata','conditions','num_conditions','stim_type');
fprintf('finished saving trial avgs for %s', stim_info)
end 
%% Takes a LONG time to run  
% save (fulldestination,'output','conditions','condition_bytrial','epochdatafile','neuraldatafile','timeAx','freqs','srate','trial_length','cond_idx','tbl','-v7.3');
% fprintf('finished saving trial avgs for %s', stim_type)




