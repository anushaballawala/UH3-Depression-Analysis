% Use with run_generate_average_or_median_table_ROI.m 
% INPUT: Epoched data for each DBS target 
% OUTPUT: Average and median values across trials for each stimulation state in table 
%  i.e. pre-stim, stim, post-stim for one contact configuration 

function generate_average_or_median_table_ROI(decompdatafile,contact_config,PatientID,experiment_name,DBS_target)

%% Assign data for each condition 

switch DBS_target 
    case 'VCVS'
        [elec1, elec25, elec36, elec47, elec8] = assign_VCVS_conditions(decomp_data,hemi); 
    case 'SCC'
        [elec1, elec234, elec25, elec36, elec47,...
            elec567,elec8] = assign_SCC_conditions(decomp_data,hemi); 
end

%% Define the different time windows of interest, in samples.
pre_stim_win = 1:4500; 
stim_win = 5001:20000; 
post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
post_stim_win2 = 25001:30000; 
post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);


%% load ROI 
load('E:/DBSTRD/DBSTRD001/Experiments/ROI_labels_DBSTRD001.mat'); 
ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG','lVPF',...
    'rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rVPF'}; 
%% get fOI data 
prompt = 'Enter frequency band of interest  '; 
FOI = input(prompt,'s'); 

switch FOI
    case 'alpha'
        freqs = 8:12; 
    case 'theta'
        freqs = 4:7; 
    case 'beta' 
        freqs = 13:18; 
    case 'low gamma' 
        freqs = 18:24; 
    case 'high gamma' 
        freqs = 24:27; 
end 
    
switch contact_config 
    case 'elec25' 
        data = elec25; 
    case 'elec36'
        data = elec36; 
    case 'elec47'
        data = elec47; 
    case '234'
        data = elec234; 
    case '567'  
        data = elec567; 
    case 'elec1'
        data = elec1; 
    case 'elec8'
        data = elec8; 
end 
%% 
% Average over freqband 
average_over_freqband = @(x) mean(x(:, :, freqs, :), 3);

data_FOI_avg_freq = average_over_freqband(data);  

if ndims(data_FOI_avg_freq) == ndims(data)
    data_FOI_avg_freq = squeeze(data_FOI_avg_freq); 
end 

%Average over time 
average_over_time = @(x, win) mean(x(:, :, win), 3);

pre_stim_avgtime = average_over_time(data_FOI_avg_freq,pre_stim_win); 
stim_avgtime = average_over_time(data_FOI_avg_freq,stim_win); 
post_stim_avgtime = average_over_time(data_FOI_avg_freq,post_stim_total_win); 

%Extract ROIs 
ROI_pre_stim_indiv_tr = generate_ROI_indiv_tr(pre_stim_avgtime); 
ROI_stim_indiv_tr = generate_ROI_indiv_tr(stim_avgtime); 
ROI_post_stim_indiv_tr = generate_ROI_indiv_tr(post_stim_avgtime); 

% Average over trials 

average_over_trials_ROI = @(x) mean(x,1); 

ROI_pre_stim_mn = average_over_trials_ROI(ROI_pre_stim_indiv_tr); 
ROI_stim_mn = average_over_trials_ROI(ROI_stim_indiv_tr) ; 
ROI_post_stim_mn = average_over_trials_ROI(ROI_post_stim_indiv_tr); 

% Take median over trials 
median_over_trials_ROI = @(x) median(x(:,:),1); 

ROI_pre_stim_med = median_over_trials_ROI(ROI_pre_stim_indiv_tr); 
ROI_stim_med = median_over_trials_ROI(ROI_stim_indiv_tr); 
ROI_post_stim_med = median_over_trials_ROI(ROI_post_stim_indiv_tr); 
%% put stuff into table and get metadata 

ROI_table_mn = vertcat(ROI_pre_stim_mn,ROI_stim_mn,ROI_post_stim_mn); 
ROI_table_med = vertcat(ROI_pre_stim_med, ROI_stim_med, ROI_post_stim_med); 

%metadata 

metadata = decomp_data.metadata; 
metadata.processing.contact_config = contact_config; 
metadata.processing.freqband = FOI; 
metadata.processing.frequencyvalues = freqs; 
metadata.processing.experiment_name = experiment_name; 
metadata.processing.filename = decompdatafile; 
metadata.processing.prestim_win = pre_stim_win; 
metadata.processing.stim_win = stim_win; 
metadata.processing.poststim_win = post_stim_total_win; 
metadata.processing.ROI_labels = ROI_labels; 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s',PatientID, experiment_name); 

if ~exist(outputdir,'dir'), mkdir(outputdir), end

tbl_mean = array2table(ROI_table_mn', 'VariableNames',...
    {'Prestim_lSCC','Stim_lSCC','PostStim_lSCC'}) ; 
tbl_mean.Properties.RowNames = ROI_labels; 

tbl_median = array2table(ROI_table_med', 'VariableNames',...
    {'Prestim_lSCC','Stim_lSCC','PostStim_lSCC'}) ; 

thisfile = sprintf('%s_%s_%s_table_ROI.mat', contact_config,experiment_name,FOI); 
fulldestination = fullfile(outputdir, thisfile); 
save(fulldestination, 'ROI_table_mn', 'ROI_table_med', 'tbl_mean',...
    'tbl_median','metadata') 
end 

%% Subtract to get pre-stim and stim-post

% perform log of data 
% z-score data as alternate 
% 
% save pre-stim and stim-post values 
% do subtraction for both pre-stim average and pre-stim single trial level 
% 
% save individual trials to generate error bars later (like in Chang paper)
% 

%% 



