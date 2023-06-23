%% 
%% take out any ROIs that werent taken care of with bipolar re-referencing. 


% load ROI labels. 

%% 

% for 001. 
remove_idx = [8,22,25,32,38,42,46,52,65,73,82,96,108,115,116,123];

% for 002.
remove_idx = [9,19,34,40,46,67,80,94,99,104,119]; 


% load data. 
path = '/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/concatenated_data/'; 

% average data across each frequency band. 

% compute t-statistic or load data with t-statistic between on on and off. 

%% take out any bad electrodes, etc. 


%% check that ROI is in the same brain region and any ones from other ROIs are removed. 


%% 
%% make histogram for each frequency band. 


% make sure gamma,beta,etc. are correct and based on what I'd chosen
% before. 

%% Make figure of on vs off (violin plot or scatter plot) for figure 2. 