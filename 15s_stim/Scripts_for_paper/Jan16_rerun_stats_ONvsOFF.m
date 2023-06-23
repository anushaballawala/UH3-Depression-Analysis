%% Goal of this script is to 1) run statistical analysis between poststim and baseline 
% 2) make a table that is input into python to make bar plots for paper. 
clc 
close all 
clear
addpath(genpath('/Users/anushaallawala/Data/')); 
%% Load file that contains data for each stimulation trial and baseline. 
% note: this data must have autocorrelation accounted for already. 

PatientID = 'DBSTRD001'; 

data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl;
clear data_file

%% Load ROI info. 

ROI_info = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/Electrodes Info/ROI_labels_mat/%s_ROI_labels.mat',...
    PatientID)); 

switch PatientID 
    % for 001
    case 'DBSTRD001'
        bipolar_remove_idx = [8,22,25,32,38,42,46,52,65,73,82,96,108,115,116,123];
        table_remove_idx = [59,60,61,132:136]; 
        del_idx = horzcat(bipolar_remove_idx,table_remove_idx); 
        % for 002 
    case 'DBSTRD002' 
        bipolar_remove_idx = [9,19,34,40,46,67,80,94,99,104,119]; 
        table_remove_idx = [57,58,102]; 
        del_idx = horzcat(bipolar_remove_idx,table_remove_idx); 
end 
        

%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_stimvsoff_conds_%s_new',PatientID)).comps;

%% Take out any trials that don't meet poststim or baseline conditions. 
nostim_idx = zeros(num_trials,1); 

for i = 1:num_trials
    has_label = @(label) contains(tbl.Stim_Labels(i), label) == 1; 
    if has_label('stim1')|| has_label('stim2')|| has_label('stim3')
        nostim_idx(i) = 0;
    else
        nostim_idx(i) = 1;
    end
end

data_nostim = data(logical(nostim_idx'),:,:); 
idx_nostim = logical(nostim_idx'); 

%get new table with poststim and bl condition labels. 
poststim_tbl = tbl(idx_nostim,:); 

%% Assign VCVS and SCC labels. 

%remove extra labels for other stim conditions. 
for i = 1:numel(comp_labels)    
    comp_labels{1,i}.label_vector = comp_labels{1,i}.label_vector(logical(nostim_idx'));     
end 

% assign labels for DBS target. 
switch PatientID 
    case 'DBSTRD001'                 
        labels = {}; 
        labels{1} = comp_labels{1,9}.label_vector ; % lSCC
        labels{2} = comp_labels{1,17}.label_vector ; % rSCC
        labels{3} = comp_labels{1,31}.label_vector ; %lVCVS
        labels{4} = comp_labels{1,37}.label_vector ; %rVCVS
    case 'DBSTRD002' 
        labels = {}; 
        labels{1} = comp_labels{1,9}.label_vector ; 
        labels{2} = comp_labels{1,17}.label_vector ; 
        labels{3} = comp_labels{1,33}.label_vector ; 
        labels{4} = comp_labels{1,41}.label_vector ;     
end 

%% Calculate z-score 

data_zscore = zeros(size(data_nostim)); 
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - mean(data_nostim(:,j,k)))./std(data_nostim(:,j,k)); 
        end 
    end 
end

%% throw out remaining channels that are bad & corresponding labels. 

ch_data = data_zscore; 
clear data_zscore
ch_data(:,del_idx,:) = []; 

ch_labels = ROI_info.summary_ch_info; 
ch_labels(bipolar_remove_idx,:) = [];

assert(size(ch_labels,1)==size(ch_data,2))

%% Combine and get ROIs of interest for analysis. 

% remove channels with gray matter and get indices for ROIs of interest. 
is_grey = contains(ch_labels.Matter, 'Grey');
is_grey_roi = @(roi) find(contains(ch_labels.ROI, roi) & is_grey); 
switch PatientID 
    case 'DBSTRD001' 
        acc_idx = is_grey_roi('ACC') ; 
        amy_idx = is_grey_roi( 'AMY');

    case 'DBSTRD002' 
        acc_idx = is_grey_roi( 'aCC') ; 
        amy_idx = is_grey_roi( 'Amy');

end 
        
vpf_idx = is_grey_roi( 'VPF') ; 
mof_idx = is_grey_roi( 'MOF') ; 
lof_idx = is_grey_roi( 'LOF'); 
dpf_idx = is_grey_roi( 'DPF') ;

%% combine all rois. 

acc  = mean(ch_data(:,acc_idx,:),2); 
amy = mean(ch_data(:,amy_idx,:),2); 
mof = mean(ch_data(:,mof_idx,:),2); 
lof = mean(ch_data(:,lof_idx,:),2);
vpf = mean(ch_data(:,vpf_idx,:),2); 

matrix_roi = [acc,amy,mof,lof,vpf]; 
ROI_labels = {'ACC','AMY','MOF','LOF','VPF'}; 

% make foi label

foi_labels = {'delta','theta','alpha','beta','lowg','highg'}; 
%% get t-statistic between off and on for SCC and VCVS respectively. 

num_perm = 1000; % change to a 1000 after verified.                                                                                                                                                                                 
% t_stat structure is num_perm+1 X data_dim X num_freq. 

%shuffled_labels = cell(numel(comp_labels),1); 
shuffled_labels = cell(numel(labels),1);
t_stat = cell(numel(labels),1);
%permutations done separately for each condition being compared, just put
%in one cell. 
for i = 1:numel(labels)
    [t_stat{i},shuffled_labels{i}] = permutations(matrix_roi,labels{i},num_perm,@Tfunc);    
end

%% Get num of observations. 

%% restructure t-stat matrix to run exchangeHT 
% testing data from stimulation experiments in left hemisphere 
p_vals_lSCC = exchange(t_stat{1});
p_vals_rSCC = exchange(t_stat{2}); 
p_vals_lVCVS = exchange(t_stat{3});
p_vals_rVCVS= exchange(t_stat{4});

%% get uncorrected values. 

p_vals_lSCC_unc = exchange_unc(t_stat{1});
p_vals_rSCC_unc = exchange_unc(t_stat{2}); 
p_vals_lVCVS_unc = exchange_unc(t_stat{3});
p_vals_rVCVS_unc = exchange_unc(t_stat{4});
%% Create table for plotting in python. 
names = ["PSD", "ROI", "FOI", "Target", "Condition"];
py_table = table([],[],[],[],[],'VariableNames', names);

for i_trial = 1:size(matrix_roi, 1)
    target = poststim_tbl.DBS_target(i_trial);
    condition = poststim_tbl.Stim_Labels(i_trial);
    for i_roi = 1:size(matrix_roi, 2)
        roi = ROI_labels{i_roi};
        for i_foi = 1:size(matrix_roi, 3)
            foi = foi_labels{i_foi};
            x = matrix_roi(i_trial, i_roi, i_foi);
            py_table = cat(1, py_table, table(x, string(roi), string(foi), string(target), string(condition), 'VariableNames', names));
        end
    end
end

%% take out prestim trials. 

py_table = py_table(py_table.Condition ~= 'prestim',:); 

%% Make table for lSCC, each ROI 

targets = {'lSCC','rSCC','lVCVS','rVCVS'}; 

for i = 1:length(targets) 

    Jan16_make_table_py_viz('ACC',targets{i},py_table,PatientID);
    Jan16_make_table_py_viz('AMY',targets{i},py_table,PatientID);
    
    
    Jan16_make_table_py_viz('MOF',targets{i},py_table,PatientID);
    Jan16_make_table_py_viz('LOF',targets{i},py_table,PatientID);
    Jan16_make_table_py_viz('VPF',targets{i},py_table,PatientID);

end 

%% troubleshoot. 
% crit_idx = (py_table.ROI == 'ACC' & py_table.Target == 'lSCC') | (py_table.ROI == 'ACC' & py_table.Target == 'BaselineFix'); 
% 
% tmp_tbl = py_table((py_table.ROI == 'ACC') & py_table.Target == 'lSCC'| py_table.Target == 'BaselineFix',:); 
% addpath(genpath('/Users/anushaallawala/Data/01-16-23-paper/'))
% writetable(tbl,sprintf('/Users/anushaallawala/Data/01-16-23-paper/%s_%s_%s_offvson.csv',PatientID,ROI_name,DBS_target));

%% 

function p_vals_all = exchange(t_stat_for_roi)
    num_ROI = size(t_stat_for_roi,2); 
    num_FOI = size(t_stat_for_roi,3); % number of spectral features of interest 
    num_obs = num_ROI*num_FOI; 
    t_stat = reshape(permute(t_stat_for_roi, [2,3,1]), [num_obs,1001]);
    [pvals,~,~] = exchangeHT(abs(t_stat),0.05,false); 
    p_vals_all = reshape(pvals(:,4),[num_ROI,num_FOI]); %grab p values from k == 4
end

function p_vals_all_unc = exchange_unc(t_stat_for_roi)
    num_ROI = size(t_stat_for_roi,2); 
    num_FOI = size(t_stat_for_roi,3); % number of spectral features of interest 
    num_obs = num_ROI*num_FOI; 
    t_stat = reshape(permute(t_stat_for_roi, [2,3,1]), [num_obs,1001]);
    [pvals,~,~] = exchangeHT(abs(t_stat),0.05,false); 
    p_vals_all_unc = reshape(pvals(:,1),[num_ROI,num_FOI]); %grab p values from k == 4
end




