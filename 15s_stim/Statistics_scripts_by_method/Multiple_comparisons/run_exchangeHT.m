%% Running multiple comparisons testing for baseline-subtracted SCC vs 
% baseline-subtracted VCVS 
clear all 
close all 
clc 
PatientID = 'DBSTRD001'; 
load(sprintf('%s_PSD_tstats_25-May-2022.mat',PatientID)); 

%% Restructure t-stat matrix to be provided as an input to exchangeHT fxn. 

% %Running multiple comparisons for each FOI separately 
% for i = 1:6
%     t_stat_for_HT{i} = permute(t_stat_left(:,:,i),[2,1]);
%     [pvals{i},acceptPW{i},acceptSIM{i}] = exchangeHT(t_stat_for_HT{i},0.05,false); 
% 
% end  

%% IF Patient is 002, remove hippocampal ROIs 

switch PatientID 
    case 'DBSTRD001' 
        idx_new = [1:2, 4:5, 8]; 
        t_stat_left = t_stat_left(:,idx_new,:); 
        t_stat_right = t_stat_right(:,idx_new,:); 
    case 'DBSTRD002' 
        idx_new = [1:5]; 
        t_stat_left = t_stat_left(:,idx_new,:); 
        t_stat_right = t_stat_right(:,idx_new,:); 
end 

%% Get num of observations. 

%num_ROI = size(t_stat_left,2); % number of regions of interest 
num_ROI = 5; 
num_FOI = 6; % number of spectral features of interest 
num_obs = num_ROI*num_FOI; 

%% restructure t-stat matrix to run exchangeHT 
% testing data from stimulation experiments in left hemisphere 

k = 4; 
t_stat_for_HT_all_left = permute(t_stat_left, [2,3,1]); 
tmp_left = reshape(t_stat_for_HT_all_left,[num_obs,1001]); 

%reshape(t_stat, [48, 1001]); 

[pvals_all_left,~,~] = exchangeHT(abs(tmp_left),0.05,false); 

p_vals_all_FOI_left = reshape(pvals_all_left(:,k),[num_ROI,num_FOI]); %grab p values from k == 4


%% Check that p-value in first k column is equivalent to what I had obtained previously after perm testing. 

%p-value from perm testing: 
disp('p-val left')
p_value_left
disp('p-val from exchangeHT')
reshape(pvals_all_left(:,1),[num_ROI,num_FOI])
%% Repeat for t_stat_right (values for stim experiments in right hemisphere - independent) 

t_stat_for_HT_all_right = permute(t_stat_right, [2,3,1]); 
tmp_right = reshape(t_stat_for_HT_all_right,[num_obs,1001]); 

%reshape(t_stat, [48, 1001]); 

[pvals_all_right,~,~] = exchangeHT(abs(tmp_right),0.05,false); 

p_vals_all_FOI_right = reshape(pvals_all_right(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

%% 

note = 'multiple comparisons for bl-subtracted SCC vs VCVS, using k = 4'; 
filename = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/mult_compare/%s_%s_mult_compare_data.mat',PatientID,date); 
save(filename,'ROI_labels1','ROI_labels','p_vals_all_FOI_left','p_vals_all_FOI_right',...
    'pvals_all_left','pvals_all_right','note'); 





