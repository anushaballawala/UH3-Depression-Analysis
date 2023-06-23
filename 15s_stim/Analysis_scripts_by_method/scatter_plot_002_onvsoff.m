addpath(genpath('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/')); 

clear 
clc 
close all 
PatientID = 'DBSTRD002'; 

load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/%s_t_stat_vals_ROI_offon_27-Jun-2022_zscore.mat',PatientID));

%% 
lSCC_tstat = squeeze(t_stat{1}(1,:,:));  
rSCC_tstat = squeeze(t_stat{2}(1,:,:));
lVCVS_tstat = squeeze(t_stat{3}(1,:,:));
rVCVS_tstat = squeeze(t_stat{4}(1,:,:));

%% 

scatter(lSCC_tstat, lVCVS_tstat);
hold on 
scatter(rSCC_tstat, rVCVS_tstat);

%% for 002 

%lSCC 
new_lSCC_tstat(1,:) = lSCC_tstat(4,:); 
new_lSCC_tstat(2,:) = lSCC_tstat(5,:); 
new_lSCC_tstat(3,:) = lSCC_tstat(6,:); 
new_lSCC_tstat(4,:) = lSCC_tstat(2,:);
new_lSCC_tstat(5,:) = lSCC_tstat(1,:);
new_lSCC_tstat(6,:) = lSCC_tstat(7,:); 
new_lSCC_tstat(7,:) = lSCC_tstat(8,:); 
new_lSCC_tstat(8,:) = lSCC_tstat(3,:);

%lVCVS 
new_lVCVS_tstat(1,:) = lVCVS_tstat(4,:); 
new_lVCVS_tstat(2,:) = lVCVS_tstat(5,:); 
new_lVCVS_tstat(3,:) = lVCVS_tstat(6,:); 
new_lVCVS_tstat(4,:) = lVCVS_tstat(2,:);
new_lVCVS_tstat(5,:) = lVCVS_tstat(1,:);
new_lVCVS_tstat(6,:) = lVCVS_tstat(7,:); 
new_lVCVS_tstat(7,:) = lVCVS_tstat(8,:); 
new_lVCVS_tstat(8,:) = lVCVS_tstat(3,:);

%rSCC
new_rSCC_tstat(1,:) = rSCC_tstat(4,:); 
new_rSCC_tstat(2,:) = rSCC_tstat(5,:); 
new_rSCC_tstat(3,:) = rSCC_tstat(6,:); 
new_rSCC_tstat(4,:) = rSCC_tstat(2,:);
new_rSCC_tstat(5,:) = rSCC_tstat(1,:);
new_rSCC_tstat(6,:) = rSCC_tstat(7,:); 
new_rSCC_tstat(7,:) = rSCC_tstat(8,:); 
new_rSCC_tstat(8,:) = rSCC_tstat(3,:);

%rVCVS
new_rVCVS_tstat(1,:) = rVCVS_tstat(4,:); 
new_rVCVS_tstat(2,:) = rVCVS_tstat(5,:); 
new_rVCVS_tstat(3,:) = rVCVS_tstat(6,:); 
new_rVCVS_tstat(4,:) = rVCVS_tstat(2,:);
new_rVCVS_tstat(5,:) = rVCVS_tstat(1,:);
new_rVCVS_tstat(6,:) = rVCVS_tstat(7,:); 
new_rVCVS_tstat(7,:) = rVCVS_tstat(8,:); 
new_rVCVS_tstat(8,:) = rVCVS_tstat(3,:);

ROI_labels_new = {'ACC','AMY','DPF','LOF','MOF','STG','SFG','VPF'}; 

%% FOI labels 

num_ROI = size(new_lSCC_tstat,1); 
foi_1 = string(repmat('delta',num_ROI,1)); 
foi_2 = string(repmat('theta',num_ROI,1)); 
foi_3 = string(repmat('alpha',num_ROI,1)); 
foi_4 = string(repmat('beta',num_ROI,1)); 
foi_5 = string(repmat('lowgamma',num_ROI,1)); 
foi_6 = string(repmat('highgamma',num_ROI,1)); 
all_f = vertcat(foi_1,foi_2,foi_3,foi_4,foi_5,foi_6);  
%% 

% ROI_tmp = repmat(ROI_labels1',6,1);
% all_ROI = ROI_tmp(:);

ROI_tmp = repmat(ROI_labels_new',6,1);
all_ROI = ROI_tmp(:);
%% Create table for plotting in seaborn 
%all_DBS = horzcat(lSCC_tstat(:),lVCVS_tstat(:),rSCC_tstat(:),rVCVS_tstat(:));

all_DBS = horzcat(new_lSCC_tstat(:),new_lVCVS_tstat(:),new_rSCC_tstat(:),new_rVCVS_tstat(:));
%% 
tbl = array2table(all_DBS,'VariableNames',{'lSCC','lVCVS','rSCC','rVCVS'});

tbl(:,5) = all_ROI; 
tbl(:,6) = cellstr(all_f);

tbl.Properties.VariableNames = {'lSCC','lVCVS','rSCC','rVCVS','ROI','FOI'}; 

filename = sprintf('/Users/anushaallawala/Python_projects/%s_offvson_tstat.csv',PatientID); 
writetable(tbl,filename);