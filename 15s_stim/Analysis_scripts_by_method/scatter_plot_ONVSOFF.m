
%% Goal of this script is to generate a scatter plot of : 
% SCC t-stat (BL vs poststim) against VCVS t-stat (BL vs poststim) to 
% observe common patterns of stimulation across neural features of
% interest across the brain. 
% Idea is to generate 6 different plots (6 diff frequencies) and color code
% ROIs in scatter plot. 
% Secondary goal: make .csv file for generating jointplot in seaborn. 

addpath(genpath('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/')); 

clear 
clc 
close all 
PatientID = 'DBSTRD001'; 

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


%% FOI labels 

num_ROI = size(lSCC_tstat,1); 
foi_1 = string(repmat('delta',num_ROI,1)); 
foi_2 = string(repmat('theta',num_ROI,1)); 
foi_3 = string(repmat('alpha',num_ROI,1)); 
foi_4 = string(repmat('beta',num_ROI,1)); 
foi_5 = string(repmat('lowgamma',num_ROI,1)); 
foi_6 = string(repmat('highgamma',num_ROI,1)); 
all_f = vertcat(foi_1,foi_2,foi_3,foi_4,foi_5,foi_6);  

%% 

ROI_tmp = repmat(ROI_labels1',6,1);
all_ROI = ROI_tmp(:);

%% Create table for plotting in seaborn 
all_DBS = horzcat(lSCC_tstat(:),lVCVS_tstat(:),rSCC_tstat(:),rVCVS_tstat(:));

%% 
tbl = array2table(all_DBS,'VariableNames',{'lSCC','lVCVS','rSCC','rVCVS'});

tbl(:,5) = all_ROI; 
tbl(:,6) = cellstr(all_f);

tbl.Properties.VariableNames = {'lSCC','lVCVS','rSCC','rVCVS','ROI','FOI'}; 

filename = sprintf('/Users/anushaallawala/Python_projects/%s_offvson_tstat.csv',PatientID); 
writetable(tbl,filename);









