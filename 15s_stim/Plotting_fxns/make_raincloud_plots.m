clear all
clc
addpath(genpath('/Users/anushaallawala/Data/DBSTRD/15s_stim/'))
%% 
PatientID  = 'DBSTRD003'; 

switch PatientID
    case 'DBSTRD001'
        load('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/DBSTRD001_PSD_tstats_25-May-2022.mat');
        roi_idx = 8; 
        foi_idx = 6; % beta 
    case 'DBSTRD002' 
        load('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/DBSTRD002_PSD_tstats_25-May-2022.mat'); 
        roi_idx = 7; 
        foi_idx = 3; % delta 
    case 'DBSTRD003' 
        load('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/DBSTRD003_PSD_tstats_25-May-2022.mat');
        roi_idx = 4;
        foi_idx = 2; % theta 
end 
        
 %% Get difference between SCC and BL, VCVS and BL 
 
 baseline_data = matrix_ROI(labels{1}==0,:,:);
 mean_bl = mean(baseline_data,1);
 
 
 % SCC.
 lSCC_data = matrix_ROI(labels{1}==1,:,:);
 rSCC_data = matrix_ROI(labels{2}==1,:,:);
 
 % VCVS.
 lVCVS_data = matrix_ROI(labels{3}==1,:,:);
 rVCVS_data = matrix_ROI(labels{4}==1,:,:);
 
 diff_lSCC = lSCC_data - mean_bl;
 diff_rSCC = rSCC_data - mean_bl;
 
 diff_lVCVS = lVCVS_data - mean_bl;
 diff_rVCVS = rVCVS_data - mean_bl;
 
%% Make raincloud plots 

addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/RainCloudPlots-master/')); 
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/cbrewer/')); 
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);

%% Plot right. 


scc_d = squeeze(diff_rSCC(:,roi_idx,foi_idx)); % VMPFC high gamma. 
vcvs_d = squeeze(diff_rVCVS(:,roi_idx,foi_idx)); % VMPFC high gamma. 


h1 = raincloud_plot(scc_d, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', 0.15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(vcvs_d, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
     'box_col_match', 1);
legend([h1{1} h2{1}], {'rSCC', 'rVCVS'});
title('rVCVS vs rSCC');
ax = gca;
ax.XAxisLocation = "top";
%set(gca,'XLim', [0 40]);
box off
camroll(90)

%% Plot left. 

scc_d_left = squeeze(diff_lSCC(:,roi_idx,foi_idx)); % VMPFC high gamma. 
vcvs_d_left = squeeze(diff_lVCVS(:,roi_idx,foi_idx)); % VMPFC high gamma. 


h1 = raincloud_plot(scc_d_left, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', 0.15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(vcvs_d_left, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
     'box_col_match', 1);
legend([h1{1} h2{1}], {'lSCC','lVCVS'});
title('lVCVS vs lSCC');
ax = gca;
ax.XAxisLocation = "top";
%set(gca,'XLim', [0 40]);
box off
camroll(90)


