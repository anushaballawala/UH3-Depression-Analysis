%% Running multiple comparisons testing for baseline-subtracted SCC vs 
% baseline-subtracted VCVS 
clear all 
close all 
clc 
PatientID = 'DBSTRD002'; 
load(sprintf('%s_t_stat_vals_ROI_offon_27-Jun-2022_zscore.mat',PatientID)); 

%% 

t_stat_lSCC = t_stat{1}; 
t_stat_rSCC = t_stat{2}; 
t_stat_lVCVS = t_stat{3}; 
t_stat_rVCVS = t_stat{4}; 

%lSCC 312-121 nan

%% 

switch PatientID 
    case 'DBSTRD001' 
        idx_new = [1:2, 4:5, 8]; 
        t_stat_lSCC = t_stat_lSCC(:,idx_new,:); 
        t_stat_rSCC = t_stat_rSCC(:,idx_new,:); 
        t_stat_lVCVS = t_stat_lVCVS(:,idx_new,:); 
        t_stat_rVCVS = t_stat_rVCVS(:,idx_new,:); 
        
    case 'DBSTRD002' 
        idx_new = [1:5]; 
        t_stat_lSCC = t_stat_lSCC(:,idx_new,:); 
        t_stat_rSCC = t_stat_rSCC(:,idx_new,:); 
        t_stat_lVCVS = t_stat_lVCVS(:,idx_new,:); 
        t_stat_rVCVS = t_stat_rVCVS(:,idx_new,:); 
end 

%% Get num of observations. 

%num_ROI = size(t_stat_left,2); % number of regions of interest 
num_ROI = 5; 
num_FOI = 6; % number of spectral features of interest 
num_obs = num_ROI*num_FOI; 

%% restructure t-stat matrix to run exchangeHT 
% testing data from stimulation experiments in left hemisphere 

%lSCC
k = 4; 
t_stat_for_HT_all_lSCC = permute(t_stat_lSCC, [2,3,1]); 
tmp_left_lSCC= reshape(t_stat_for_HT_all_lSCC,[num_obs,1001]); 
[pvals_lSCC,~,~] = exchangeHT(abs(tmp_left_lSCC),0.05,false); 
p_vals_all_lSCC = reshape(pvals_lSCC(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

%rSCC
t_stat_for_HT_all_rSCC = permute(t_stat_rSCC, [2,3,1]); 
tmp_left_rSCC= reshape(t_stat_for_HT_all_rSCC,[num_obs,1001]); 
[pvals_rSCC,~,~] = exchangeHT(abs(tmp_left_rSCC),0.05,false); 
p_vals_all_rSCC = reshape(pvals_rSCC(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

%% Repeat for VCVS

%lVCVS
k = 4; 
t_stat_for_HT_all_lVCVS = permute(t_stat_lVCVS, [2,3,1]); 
tmp_left_lVCVS= reshape(t_stat_for_HT_all_lVCVS,[num_obs,1001]); 
[pvals_lVCVS,~,~] = exchangeHT(abs(tmp_left_lVCVS),0.05,false); 
p_vals_all_lVCVS = reshape(pvals_lVCVS(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

%rVCVS
t_stat_for_HT_all_rVCVS = permute(t_stat_rVCVS, [2,3,1]); 
tmp_left_rVCVS= reshape(t_stat_for_HT_all_rVCVS,[num_obs,1001]); 
[pvals_rVCVS,~,~] = exchangeHT(abs(tmp_left_rVCVS),0.05,false); 
p_vals_all_rVCVS = reshape(pvals_rVCVS(:,k),[num_ROI,num_FOI]); %grab p values from k == 4

alt_p_rVCVS = (squeeze(mean(abs(t_stat_rVCVS)>=abs(t_stat_rVCVS(1,:,:))))); 


%% 

note = 'multiple comparisons for OFF v ON, using k = 4'; 
filename = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/off_vs_on_multcompare/%s_%s_OFFvON_mult_compare_data.mat',PatientID,date); 
save(filename,'ROI_labels1','p_vals_all_lSCC','p_vals_all_rSCC','p_vals_all_lVCVS','p_vals_all_rVCVS'); 





