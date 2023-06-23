%% DBStarget_hypothesistest.m 

%%% FXN: Permutation and hypothesis testing for all combinations of 
%%%      15-s experiments in 001 and 002. 

%%% ***Bug-fixes remaining: 
%%% 1. check whether individual current configurations can be tested 
%%% 2. check if lead  hemisphere data can be combined 
%%% 3. test poststim 
%%% 4. test stim on 
%%% 5. Save list of channels that are significant 
%%% 6. save experiment tested in metadata. 

%% Load file that contains data for each trial

clc 
close all                                                                       
clear
PatientID = 'DBSTRD001'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
%data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
switch PatientID 
    case 'DBSTRD001'
        load(sprintf('%s_Jan-04-2023_allroi_tstat_ch.mat', PatientID)); 
    case 'DBSTRD002'
        load(sprintf('%s_06-Jan-2023_allroi_tstat_ch.mat', PatientID)); 
end 
%% get ROI. 

acc = mean(ch_data(:,acc_idx,:),2); 
amy = mean(ch_data(:,amy_idx,:),2); 
lof = mean(ch_data(:,lof_idx,:),2); 
mof = mean(ch_data(:,mof_idx,:),2); 
vpf = mean(ch_data(:,vpf_idx,:),2); 

% concatenate into one. 
all_roi_data = horzcat(acc,amy,lof,mof,vpf); 

roi_labels_new = {'acc','amy','lof','mof','vpf'}; 

%% foi labels 

foi_labels_new = {'delta','theta','alpha','beta','lowg','highg'}; 

%% separate data matrices. 

baseline_idx = find(contains(stim_tbl.Stim_Labels,'Baseline')==1) ; 
rSCC_idx = find(contains(stim_tbl.Stim_Labels,'poststim')==1 & contains(stim_tbl.DBS_target,'rSCC')==1); 
rVCVS_idx = find(contains(stim_tbl.Stim_Labels,'poststim')==1 & contains(stim_tbl.DBS_target,'rVCVS')==1); 
lSCC_idx = find(contains(stim_tbl.Stim_Labels,'poststim')==1 & contains(stim_tbl.DBS_target,'lSCC')==1); 
lVCVS_idx = find(contains(stim_tbl.Stim_Labels,'poststim')==1 & contains(stim_tbl.DBS_target,'lVCVS')==1); 

lSCC_data = all_roi_data(lSCC_idx,:,:); 
rSCC_data = all_roi_data(rSCC_idx,:,:); 

lVCVS_data = all_roi_data(lVCVS_idx,:,:); 
rVCVS_data = all_roi_data(rVCVS_idx,:,:); 

baseline_data = all_roi_data(baseline_idx,:,:); 
%% Get mean of baseline data. 
 
baseline_mean = mean(baseline_data,1); 

% subtract SCC and VCVS by baseline data. 
diff_lSCC = lSCC_data - baseline_mean;
diff_lVCVS = lVCVS_data - baseline_mean; 
diff_rSCC = rSCC_data - baseline_mean;
diff_rVCVS = rVCVS_data - baseline_mean; 

diff_data_left = [diff_lSCC;diff_lVCVS]; 
diff_data_right = [diff_rSCC;diff_rVCVS]; 
%% save 
% 
% filename_left = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/Figure 3/%s_zscore_bl_poststim_left_SCC_VCVS_%s.mat',PatientID,date); 
% save(filename_left,'diff_data_left'); 
% 
% filename_right = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/Figure 3/%s_zscore_bl_poststim_right_SCC_VCVS_%s.mat',PatientID,date); 
% save(filename_right,'diff_data_right');

%% Assign colors 

addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/RainCloudPlots-master/')); 
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/cbrewer/')); 
addpath(genpath(    '/Users/anushaallawala/Documents/MATLAB/Github/Packages/hex_and_rgb_v1.1.1')); 

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);

VCVS_color = hex2rgb('#6b4266'); 
SCC_color = hex2rgb('#64809a'); 
%baseline_color = hex2rgb('#d3d4da'); 
baseline_color = [0.2196    0.2157    0.2157];

%% get data and plot raincloud plots for all fois and rois. 

for roi_idx = 1:2%numel(roi_labels_new)
    
    for foi_idx = 1:2%numel(foi_labels_new)
        figure() 
        % run loop for right.
        scc_d = squeeze(rSCC_data(:,roi_idx,foi_idx));
        vcvs_d = squeeze(rVCVS_data(:,roi_idx,foi_idx));
        baseline_d = squeeze(baseline_data(:,roi_idx,foi_idx));
        % plot for right.
        % BASELINE.
        h1 = raincloud_plot(baseline_d, 'box_on', 1, 'color', baseline_color, 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
            'box_col_match', 1, 'line_width_density', 0.75, 'bxfacecl', baseline_color, 'line_width',1.25,...
            'lwr_bnd', 4);
        % VCVS. 
        h2 = raincloud_plot(vcvs_d, 'box_on', 1, 'color', VCVS_color, 'alpha', 0.6,...
            'box_dodge', 1, 'box_dodge_amount', .47, 'dot_dodge_amount', .47,...
            'box_col_match', 1, 'line_width_density', 0.75, 'bxfacecl',VCVS_color, 'line_width',1.25,...
            'lwr_bnd', 4);
        % SCC. 
        h3 = raincloud_plot(scc_d, 'box_on', 1, 'color', SCC_color, 'alpha', 0.6,...
            'box_dodge', 1, 'box_dodge_amount', 0.79, 'dot_dodge_amount', .79,...
            'box_col_match', 1, 'line_width_density', 0.75, 'bxfacecl', SCC_color, 'line_width',1.25,...
            'lwr_bnd', 4);
                           
        legend([h1{1} h2{1} h3{1}], {'baseline', 'rVCVS','rSCC'});
        title(sprintf('rVCVS vs rSCC %s %s',roi_labels_new{roi_idx},foi_labels_new{foi_idx}));
        ax = gca;
        ax.XAxisLocation = "top";
        box off
        camroll(90)
        %xlim([-0.5 0.5])
        addpath(genpath('/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/'));
        filepath_right = '/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/Figure 3/'; 
        filename_right = sprintf('%s_%s_new_right_raincloud_%s_%s.svg',PatientID,date,roi_labels_new{roi_idx},foi_labels_new{foi_idx});
        full_file_right = sprintf('%s%s', filepath_right, filename_right); 
        % save 
       % saveas(gcf,full_file_right); 
        sprintf('saved %s',filename_right); 
        %close all  
        %-------------------% 
        % get data for left. 
        scc_d_left = squeeze(lSCC_data(:,roi_idx,foi_idx));
        vcvs_d_left = squeeze(lVCVS_data(:,roi_idx,foi_idx));
        baseline_d = squeeze(baseline_data(:,roi_idx,foi_idx));
%         scc_d_left = squeeze(diff_lSCC(:,roi_idx,foi_idx)); % VMPFC high gamma.
%         vcvs_d_left = squeeze(diff_lVCVS(:,roi_idx,foi_idx)); % VMPFC high gamma.
%       plot for left. 
        figure() 
        % BASELINE.
        h1 = raincloud_plot(baseline_d, 'box_on', 1, 'color', baseline_color, 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
            'box_col_match', 1, 'line_width_density', 0.75, 'bxfacecl', baseline_color, 'line_width',1.25,...
             'lwr_bnd', 4);
        % VCVS. 
        h2 = raincloud_plot(vcvs_d_left, 'box_on', 1, 'color', VCVS_color, 'alpha', 0.7,...
            'box_dodge', 1, 'box_dodge_amount', .47, 'dot_dodge_amount', .47,...
            'box_col_match', 1, 'line_width_density', 0.75, 'bxfacecl',VCVS_color, 'line_width',1.25,...
             'lwr_bnd', 4);
        % SCC. 
        h3 = raincloud_plot(scc_d_left, 'box_on', 1, 'color', SCC_color, 'alpha', 0.7,...
            'box_dodge', 1, 'box_dodge_amount', 0.79, 'dot_dodge_amount', .79,...
            'box_col_match', 1, 'line_width_density', 0.75, 'bxfacecl', SCC_color, 'line_width',1.25,...
             'lwr_bnd', 4);
                           
        legend([h1{1} h2{1} h3{1}], {'baseline', 'lVCVS','lSCC'});
        title(sprintf('lVCVS vs lSCC %s %s',roi_labels_new{roi_idx},foi_labels_new{foi_idx}));
        ax = gca;
        ax.XAxisLocation = "top";
        box off
        camroll(90)
        %xlim([-0.5 0.5])
        filepath_left = '/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/Figure 3/'; 
        filename_left = sprintf('%s_%s_new_left_raincloud_%s_%s.svg',PatientID,date,roi_labels_new{roi_idx},foi_labels_new{foi_idx});
        full_file_left = sprintf('%s%s', filepath_left, filename_left); 
        % save . 
        saveas(gcf,full_file_left); 
        sprintf('saved %s',filename_left); 
        
        %close all 
    end
end

%% 



