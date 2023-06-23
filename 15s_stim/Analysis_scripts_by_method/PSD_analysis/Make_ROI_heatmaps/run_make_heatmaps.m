%% run_make_heatmaps.m - Run this to generate heatmaps for visualization 
% of spectral power in each frequency band of interest in ROIs for each
% contact configuration per DBS lead

% Additional info: This script needs a "make_heatmaps" fxn and needs to be
% rewritten as old fxns have been archived as of 10/15/21
% Inputs: N/A
% Outputs: Spectral power from Chronux multitaper fxn with frequencies, 
% and heatmaps for each DBS lead 
% Dependencies: Chronux toolbox, UH3 github repo 
% Sub-functions: N/A

% Anusha Allawala, 9
/2021

%------------ START OF CODE ---------------% %% 

% This script relies on the run_setup_15s script,
% generate_heatmap_mtspectrumc_power script 
%% setup input parameters & run setup file 
clear all 
close all
PatientID = 'DBSTRD001'; DBS_target = 'VCVS'; hemi = 'l'; stim_freq = 130; 
FOI = input('enter FOI ','s'); 
[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names, upper_freq,lower_freq,metadata,ROI_labels,num_contact_configs,...
    experiment_name,elec1,elec25,elec36,elec47,...
    elec8] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, FOI); 

%% get z-mean for recording 

        % get mean across the entire recording for z-scoring purposes. 
        elec_avg = @(elec)(mean(mean(elec,3),1));% fxn handle 

        z_mean_e1 =  elec_avg(elec1); 
        z_mean_e25 = elec_avg(elec25); 
        z_mean_e36 = elec_avg(elec36); 
        z_mean_e47 = elec_avg(elec47); 
        z_mean_e8 = elec_avg(elec8); 

        all_z_e = vertcat(z_mean_e1,z_mean_e25,z_mean_e36,z_mean_e47,z_mean_e8); 
        z_mean = mean(all_z_e,1);

        %get std 
        elec_std = @(elec)(std(std(elec,0,3),1));

        z_std_e1 = elec_std(elec1); 
        z_std_e25 = elec_std(elec25); 
        z_std_e36 = elec_std(elec36); 
        z_std_e47 = elec_std(elec47); 
        z_std_e8 = elec_std(elec8); 

        all_z_std_e = vertcat(z_std_e1,z_std_e25,z_std_e36,z_std_e47,z_std_e8);
        z_std = std(all_z_std_e,0,1); 

%% run through analysis for each configuration 
tic
tbl_all_data_stimpre = [];
tbl_all_data_postpre =[]; 
tbl_allch_alldata_stimpre = [];
tbl_allch_alldata_postpre = [];

for i = 1:num_contact_configs
    contact_config = all_config_names{i}; 
    switch contact_config
        case 'elec25'
            data = elec25;
        case 'elec1'
            data = elec1;
        case 'elec8'
            data = elec8;
        case 'elec36'
            data = elec36;
        case 'elec47'
            data = elec47;

    end 
    disp(all_config_names{i}); 
    names{i} = all_config_names{i}; 

[metadata,tbl_allch_stimpre, tbl_allch_postpre,...
    tbl_ROI_stimpre, tbl_ROI_postpre] = generate_heatmap_mtspectrumc_power(PatientID,metadata,data,z_mean,z_std,...
    lower_freq,upper_freq, ROI_labels,contact_config,ch_labels); 

tbl_all_data_stimpre = [tbl_all_data_stimpre;{tbl_ROI_stimpre}]; 
tbl_all_data_postpre = [tbl_all_data_postpre;{tbl_ROI_postpre}];

tbl_allch_alldata_stimpre = [tbl_allch_alldata_stimpre;{tbl_allch_stimpre}];
tbl_allch_alldata_postpre = [tbl_allch_alldata_postpre;{tbl_allch_postpre}];
end 
toc
%% Combine table into one this is so stupid i hate matlab 

ROI_tbl_all_configs = join(tbl_all_data_stimpre{1,1},tbl_all_data_stimpre{2,1});
tmp = join(ROI_tbl_all_configs,tbl_all_data_stimpre{3,1}); 
tmp = join(tmp,tbl_all_data_stimpre{4,1});
tmp = join(tmp,tbl_all_data_stimpre{5,1});


%% Generate heatmaps 
 
matrix_for_heatmap = horzcat(tmp.stim_pre_elec1,tmp.stim_pre_elec25,...
    tmp.stim_pre_elec36,tmp.stim_pre_elec47,...
    tmp.stim_pre_elec8); 

% matrix_for_heatmap(14,:) = []; %******
% ROI_labels(:,14) = []; %*****
%% 
figure('Position',[500 300 800 600]) 
h = heatmap(matrix_for_heatmap); 

colormap redbluecmap ; 
%co.YDisplayLabels = 'z-scored power: stim minus prestim';
h.YDisplayLabels = char(ROI_labels); 
h.XDisplayLabels = char(all_config_names); 
h.ColorLimits = [-3 3]; 
title(sprintf('stim minus pre %s %s',FOI,experiment_name));

filename = sprintf('%s_%s_%s,stimminuspre1.png',PatientID,experiment_name,FOI); 
saveas(gcf, filename)
% 
% %% group data into table ?
% 
%% post minus pre 
clear tmp ; clear ROI_tbl_all_configs; 
ROI_tbl_all_configs = join(tbl_all_data_postpre{1,1},tbl_all_data_postpre{2,1});
tmp = join(ROI_tbl_all_configs,tbl_all_data_postpre{3,1}); 
tmp = join(tmp,tbl_all_data_postpre{4,1});
tmp = join(tmp,tbl_all_data_postpre{5,1});


%% Generate heatmaps 
 
matrix_for_heatmap_postpre = horzcat(tmp.post_pre_elec1,tmp.post_pre_elec25,...
    tmp.post_pre_elec36,tmp.post_pre_elec47,...
    tmp.post_pre_elec8); 

% matrix_for_heatmap_postpre(14,:) = []; 
%% 

figure('Position',[500 300 800 600]) 
h = heatmap(matrix_for_heatmap_postpre); 

colormap redbluecmap ; 
h.YDisplayLabels = char(ROI_labels); 
h.XDisplayLabels = char(all_config_names); 
h.ColorLimits = [-3 3]; 
title(sprintf('post minus pre %s %s',FOI,experiment_name));

filename = sprintf('%s_%s_%s,postminuspre.png',PatientID,experiment_name,FOI); 
saveas(gcf, filename)

%% Save data 

outputdir = sprintf('E:\DBSTRD\%s\Experiments\15s_stim\Processed Data\',PatientID); 
filename = sprintf('%s\ROI_and_allch_%s_%s.mat',outputdir,experiment_name,FOI); 
save(filename,'metadata','ROI_labels','z_mean','z_std','tbl_all_data_stimpre','tbl_all_data_postpre',...
    'tbl_allch_alldata_stimpre','tbl_allch_alldata_postpre');


% %% Run stats? 