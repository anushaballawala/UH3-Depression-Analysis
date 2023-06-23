% load data 

% - look at channel labels that were significant for 001 and 002 
% - look at rVCVS and lVCVS vs SCC for stim on 
% - look at rVCVS and lVCVS vs SCC for stim on vs stim off (only keep VCVS
% e.g.) 

%% Load file that contains data and descriptive labels for each trial. 
clear 
clc 
close all 

PatientID = 'DBSTRD002'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data = data_file.all_data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_state_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 
%Load channels labels that will be used later. 
load(sprintf('%s_goodch.mat',PatientID)) 
%*** add new path for 002 good channels 

%% **** Fix this in the preprocessing code, temporary code for now ***** 
switch PatientID 
    case 'DBSTRD001'
        disp('change nothing')
    case 'DBSTRD002'
        tmp_idx = [1:57,59:66,68:101,103:129];
        data = data(:,tmp_idx,:); 
end 

%% Convert descriptive labels to binary labels in a vector. 

% test baseline data vs a specific DBS target or R/L DBS target 
comparison = 'DBSvsbaseline'; 

%Choose DBS target or hemisphere of interest. 
list = {'VCVS','SCC','lVCVS', 'rVCVS', 'lSCC', 'rSCC', 'r', 'l'};                    
[indx,tf] = listdlg('ListString',list, 'ListSize',[150,250]);
target = list{indx}; 
disp(target)
%your two conditions! 
cond1 = target; 

cond2 = 'Baseline'; % this does not change! 
disp(cond2)

% Choose stim state 
stim_list = {'stimall', 'poststim'}; 
[idx_stim,tf_stim] = listdlg('ListString', stim_list, 'ListSize',[150,250]); 
stim_cond = stim_list{idx_stim}; 
disp(stim_cond)

%Choose individual or all configurations 
elec_type_list = {'all', 'indiv'}; 
[indx_e,tf_e] =  listdlg('ListString', elec_type_list, 'ListSize',[150,250]); 
elec_cond_type = elec_type_list{indx_e};
disp(elec_cond_type)

% Choose electrode configuration or all of them 
elec_config_list = {'elec1', 'elec25', 'elec47', 'elec36', 'elec234', 'elec567', 'elec8','allelec'}; 
[idx_elec,tf_elec] = listdlg('ListString', elec_config_list, 'ListSize',[150,250]); 
elec_config = elec_config_list{idx_elec}; 
disp(elec_config)

%% add to metadata.
cd('/Users/anushaallawala/Documents/MATLAB/Github/UH3-Depression-Analysis/15s_stim/Statistics_scripts_by_method')

[s_git,git_hash_string] = system('git rev-parse HEAD'); 
metadata.comparison = comparison; 
metadata.cond1 = cond1; 
metadata.cond2 = cond2; 
metadata.stim_cond = stim_cond; 
metadata.elec_cond_type = elec_cond_type; 
metadata.elec_config = elec_config; 
metadata.tbl = tbl; 
metadata.gitinfo.s = s_git; 
metadata.gitinfo.git_hash_string = git_hash_string;
metadata.good_ch_labels = good_ch_labels; 
metadata.good_ch_idx = good_ch_idx; 

%% insert switch case statement for which function to run 


switch comparison 
    case 'DBSvsbaseline'

        labels = get_labels_DBSvsbaseline(comparison, num_trials, stim_cond, tbl, elec_cond_type,...
            cond1, cond2); 
    case 'somethingelse'
        disp('code not written yet')
        % case for DBStarget vs DBS target 
        % case for stimstate (baseline or prestim vs DBS stim or poststim)
end 

%% Get mean out for condition 1 

% !!!!!!!!CONDITION 1 IS ALWAYS 1 !!!!!!!!

bin_labels = labels == 1; 
%idx out condition 1 
data_cond1 = data(bin_labels,:,:); 
%take mean for each channel 
data_cond1_mean = squeeze(mean(data_cond1,1)); 

%% Get mean out for condition 2 

% !!!!!!!!CONDITION 2 IS ALWAYS 0 !!!!!!!!

bin_labels = labels == 0; 
%idx out condition 2 
data_cond2 = data(bin_labels,:,:); 
data_cond2_mean = squeeze(mean(data_cond2,1)); 
%take mean for each channel 
%% Define which type of data, Channel or ROI 

data_type_list = {'ROI', 'Ch'}; 
[idx_data_list,tf_data] = listdlg('ListString', data_type_list, 'ListSize',[150,250]); 
data_list_config = data_type_list{idx_data_list}; 
disp(data_list_config)

data_dim = data_list_config; 
 
if strcmp(data_dim,'ROI')==1 
    % get mean for plotting later. 
    [ROI_data_cond1,~,~,~] = generate_ROI_from_ch(data_cond1_mean,PatientID,1); 
    [ROI_data_cond2,~,~,~] = generate_ROI_from_ch(data_cond2_mean,PatientID,1); 
    % get data that will be input into permutations fxn. 
    [matrix_ROI,ROI_labels,...
    ~, ~] = generate_ROI_from_ch(data,PatientID,2);
end 
%% Permutation testing 
           
num_perm = 1000; 
% t_stat structure is num_perm+1 X data_dim X num_freq. 


switch data_dim 
    case 'Ch'
        [t_stat,shuffled_labels_1] = permutations(data,labels,num_perm,@Tfunc); 
    case 'ROI' 
        [t_stat,shuffled_labels_1] = permutations(matrix_ROI,labels,num_perm,@Tfunc); 
end 




%% Hypothesis testing - get p-value 

% calculate p_value by looking at comparison of t-stat values and dividing
% by total number of permuations 

% Taking absolute value to check against both extreme ends of sample
% population. 
p_value = (squeeze(mean(abs(t_stat)>=abs(t_stat(1,:,:))))); 

%% Find channels that are significantly different. 

alpha = 0.05; 
metadata.alpha = alpha; 

p_idx_delta = find(p_value(:,1)<alpha) ; 
p_idx_theta = find(p_value(:,2)<alpha); 
p_idx_alpha = find(p_value(:,3)<alpha); 
p_idx_beta = find(p_value(:,4)<alpha); 
p_idx_lowgamma = find(p_value(:,5)<alpha); 
p_idx_highgamma = find(p_value(:,6)<alpha); 

%% Label channels that are significantly different 
switch data_dim 
    case 'Ch'
        ch_delta = good_ch_labels(p_idx_delta);
        ch_theta = good_ch_labels(p_idx_theta);
        ch_alpha = good_ch_labels(p_idx_alpha);
        ch_beta = good_ch_labels(p_idx_beta);
        ch_lowgamma = good_ch_labels(p_idx_lowgamma);
        ch_highgamma  = good_ch_labels(p_idx_highgamma);

        p_delta = p_value(p_idx_delta,1); 
        p_theta = p_value(p_idx_theta,2); 
        p_alpha = p_value(p_idx_alpha,3); 
        p_beta = p_value(p_idx_beta,4); 
        p_lowgamma = p_value(p_idx_lowgamma,5); 
        p_highgamma = p_value(p_idx_highgamma,6); 
        
    case 'ROI'
        if ~isempty(p_idx_delta)
            ROI_delta = ROI_labels(p_idx_delta);
            
        end
        if ~isempty(p_idx_theta)
            ROI_theta = ROI_labels(p_idx_theta);
        end
        if ~isempty(p_idx_alpha)
            ROI_alpha = ROI_labels(p_idx_alpha);
        end
        if ~isempty(p_idx_beta)
            ROI_beta = ROI_labels(p_idx_beta);
        end
        if ~isempty(p_idx_lowgamma)
            ROI_lowgamma = ROI_labels(p_idx_lowgamma);
        end
        if ~isempty(p_idx_highgamma)
            ROI_highgamma = ROI_labels(p_idx_highgamma);
        end
        
end 

                
%% select which figures to plot 

list_freqs = {'delta','theta','alpha','beta','lowgamma','highgamma'}; 
for i = 1:6
    fig = list_freqs{i}; 
    switch fig
        case 'delta'
            p_idx = p_idx_delta;
            freqband_idx = 1; 
            freqband = list_freqs{1}; 
        case 'theta'
            p_idx = p_idx_theta; 
            freqband_idx = 2; 
            freqband = list_freqs{2}; 
        case 'alpha'
            p_idx = p_idx_alpha; 
            freqband_idx = 3; 
            freqband = list_freqs{3}; 
        case 'beta' 
            p_idx = p_idx_beta; 
            freqband_idx = 4; 
            freqband = list_freqs{4}; 
        case 'lowgamma'
            p_idx = p_idx_lowgamma; 
            freqband_idx = 5; 
            freqband = list_freqs{5}; 
        case 'highgamma'
            p_idx = p_idx_highgamma; 
            freqband_idx = 6; 
            freqband = list_freqs{6}; 
    end 
    switch data_dim 
        case 'Ch' 
            
        plot_scatter_vs(PatientID,metadata,p_idx,freqband_idx,cond1,cond2,freqband,...
            data_cond1_mean,data_cond2_mean,'Ch')
        
        case 'ROI'
        plot_scatter_vs(PatientID,metadata,p_idx,freqband_idx,cond1,cond2,freqband,...
            ROI_data_cond1,ROI_data_cond2,'ROI')
    end 
end 

%% Save data 

switch data_dim 
    case 'Ch'
        cd(sprintf('/Users/anushaallawala/Data/stats_results/two-tail_pos_neg_Tstat/pvalues/%s',data_dim))
        filename = sprintf('%s_%s_%s_%s_%s-elec_%s.mat',PatientID,comparison,target,stim_cond,elec_cond_type,elec_config); 
        save(filename,'num_perm','alpha','p_idx_delta','p_idx_theta',...
            'p_idx_alpha','p_idx_beta','p_idx_lowgamma','p_idx_highgamma',...
            'p_delta','p_theta','p_alpha','p_beta','p_lowgamma','p_highgamma',...
            'p_value','t_stat','metadata','ch_delta','ch_theta','ch_alpha',...
            'ch_beta','ch_lowgamma','ch_highgamma','data_cond1','data_cond1_mean',...
            'data_cond2','data_cond2_mean','shuffled_labels'); 
    case 'ROI'
        cd(sprintf('/Users/anushaallawala/Data/stats_results/two-tail_pos_neg_Tstat/pvalues/%s',data_dim))
        filename = sprintf('%s_%s_%s_%s_%s-elec_%s.mat',PatientID,comparison,target,stim_cond,elec_cond_type,elec_config); 
        save(filename,'num_perm','alpha','p_idx_delta','p_idx_theta',...
            'p_idx_alpha','p_idx_beta','p_idx_lowgamma','p_idx_highgamma',...
            'p_delta','p_theta','p_alpha','p_beta','p_lowgamma','p_highgamma',...
            'p_value','t_stat','metadata','ROI_data_cond1','ROI_data_cond2',...
            'matrix_ROI','ROI_delta','ROI_highgamma','ROI_lowgamma','ROI_theta','ROI_alpha','shuffled_labels'); 

disp('saved')
end 
%% 

%for i = 1:p_idx_theta(1:end)
%     figure() 
%     h1 = histogram(squeeze(data_cond1(:,p_idx_theta(i),2)),'BinWidth',1,'FaceAlpha',0.55); 
%     h1.FaceColor = [0,0.447058823529412,0.741176470588235];
%     hold on 
%     h2 = histogram(squeeze(data_cond2(:,p_idx_theta(i),2)),'BinWidth',1,'FaceAlpha',0.55); 
%     h2.FaceColor = [0.188235294117647,0.729411764705882,0.443137254901961];
% 
% %     h1.Normalization = 'probability';
% %     h1.BinWidth = 1;
% %     h2.Normalization = 'probability';
% %     h2.BinWidth = 1;
%     
%     title(good_ch_labels(p_idx_theta(i)))
% 
% %end 
% 
% %% 
% figure()
% scatter(1,squeeze(data_cond1(:,p_idx_theta(i),2)))
% hold on 
% scatter(1.5, squeeze(data_cond2(:,p_idx_theta(i),2)))
% 
% %% Plot means 
% 
% figure() 
% 
% scatter(data_cond1_mean(p_idx_beta,4),data_cond2_mean(p_idx_beta,4))
% refline
% ylim([-30 30])
% xlim([-30 30])
%% 
%% 
% 
% ROI_cond1_delta = ROI_data_cond1(p_idx_delta,1); 
% ROI_cond2_delta = ROI_data_cond2(p_idx_delta,1);
% delta = ROI_cond1_delta - ROI_cond2_delta; 
% 
% ROI_cond1_theta = ROI_data_cond1(p_idx_theta,2); 
% ROI_cond2_theta = ROI_data_cond2(p_idx_theta,2);
% theta = ROI_cond1_theta - ROI_cond2_theta; 
% 
% ROI_cond1_alpha = ROI_data_cond1(p_idx_delta,3); 
% ROI_cond2_alpha = ROI_data_cond2(p_idx_delta,3);
% alpha = ROI_cond1_alpha - ROI_cond2_alpha; 
% 
% ROI_cond1_beta = ROI_data_cond1(p_idx_beta,4); 
% ROI_cond2_beta = ROI_data_cond2(p_idx_beta,4);
% beta = ROI_cond1_beta - ROI_cond2_beta; 
% 
% ROI_cond1_lowgamma = ROI_data_cond1(p_idx_lowgamma,5); 
% ROI_cond2_lowgamma = ROI_data_cond2(p_idx_lowgamma,5);
% lowgamma = ROI_cond1_lowgamma - ROI_cond2_lowgamma; 
% 
% ROI_cond1_highgamma = ROI_data_cond1(p_idx_highgamma,6); 
% ROI_cond2_highgamma = ROI_data_cond2(p_idx_highgamma,6);
% highgamma = ROI_cond1_highgamma - ROI_cond2_highgamma; 




