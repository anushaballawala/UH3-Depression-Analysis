% Goal of this script is to bin 15s data into 1-s bins and generate scatter plots 

clear all 
PatientID = 'DBSTRD001'; 
hemi = 'l' ; 
DBS_target = 'VCVS'; 
freq = 130; 
stim_freq = freq; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
% Load single trial time series data 
decompdatafile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s/15s_stim_all_currdir_singletrial_%s.mat',...
    PatientID,experiment_name,experiment_name); 
decomp_data = load(decompdatafile); 
fs = decomp_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(decomp_data.metadata.preprocessing.GoodChannelLabels); 

%% Assign data for each condition 

switch DBS_target 
    case 'VCVS'
        [elec1, elec25, elec36, elec47, elec8] = assign_VCVS_conditions(decomp_data,hemi); 
    case 'SCC'
        [elec1, elec234, elec25, elec36, elec47,...
            elec567,elec8] = assign_SCC_conditions(decomp_data,hemi); 
end

%% Define the different time windows of interest, in samples.
pre_stim_win = 1:5000; 
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
    
data = elec25; 
contact_config = 'elec25'; %%%%%%% ADD SWITCH STATEMENT TO MAKE MODULAR FOR OTHER CURRENT DIRS 
%% 
% Average over freqband 
average_over_freqband = @(x) mean(x(:, :, freqs, :), 3);

data_FOI_avg_freq = average_over_freqband(data);  

if ndims(data_FOI_avg_freq) == ndims(data)
    data_FOI_avg_freq = squeeze(data_FOI_avg_freq); 
end 

%%  

pre_stim_data = data_FOI_avg_freq(:,:,pre_stim_win); 
stim_data = data_FOI_avg_freq(:,:,stim_win); 
post_stim_data = data_FOI_avg_freq(:,:,post_stim_total_win); 

%% Bin data into 1000 sample bins for pre_stim 

pre_stim_bins = [2500;5000]; 
bin1 = pre_stim_data(:,:,1:tmp(1)); 
bin2 = pre_stim_data(:,:,tmp(1):tmp(2)-1); 
bin3 = pre_stim_data(:,:,tmp(2):tmp(3)-1); 
bin4 = pre_stim_data(:,:,tmp(3):tmp(4)-1); 
bin5 = pre_stim_data(:,:,tmp(4):tmp(5)-1); 
pre_stim_bins = vertcat(bin1,bin2,bin3,bin4,bin5); 

%% Bin data into 1000 sample bins for stim 

stim_win_bins = [1000:1000:15000];

bin1 = stim_data(:,:,1:stim_win_bins(1)); 
bin2 = stim_data(:,:,stim_win_bins(1):stim_win_bins(2)-1); 
bin3 = stim_data(:,:,stim_win_bins(2):stim_win_bins(3)-1); 
bin4 = stim_data(:,:,stim_win_bins(3):stim_win_bins(4)-1); 
bin5 = stim_data(:,:,stim_win_bins(4):stim_win_bins(5)-1); 
bin6 = stim_data(:,:,stim_win_bins(5):stim_win_bins(6)-1); 
bin7 = stim_data(:,:,stim_win_bins(6):stim_win_bins(7)-1); 
bin8 = stim_data(:,:,stim_win_bins(7):stim_win_bins(8)-1); 
bin9 = stim_data(:,:,stim_win_bins(8):stim_win_bins(9)-1); 
bin10 = stim_data(:,:,stim_win_bins(9):stim_win_bins(10)-1); 
bin11 = stim_data(:,:,stim_win_bins(10):stim_win_bins(11)-1); 
bin12 = stim_data(:,:,stim_win_bins(11):stim_win_bins(12)-1); 
bin13 = stim_data(:,:,stim_win_bins(12):stim_win_bins(13)-1); 
bin14 = stim_data(:,:,stim_win_bins(13):stim_win_bins(14)-1); 
bin15 = stim_data(:,:,stim_win_bins(14):stim_win_bins(15)-1); 

stim_bins = vertcat(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,...
    bin11,bin12,bin13,bin14,bin15); 

%% bin post stim data 

post_win_bins = [1000:1000:9000]; 




    







post_stim_bins = [1000:1000:9000];












%% Average each bin over time 





% Extract ROIs 




% generate scatter plot 



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





