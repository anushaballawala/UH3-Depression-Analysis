%% 

% This script relies on the run_setup_15s script,
% generate_heatmap_mtspectrumc_power script 
%% setup input parameters & run setup file 
clear all 
close all
PatientID = 'DBSTRD002'; DBS_target = 'VCVS'; hemi = 'l'; stim_freq = 130; 

[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,ROI_labels,...
    experiment_name,elec1,elec25,elec36,elec47,...
    elec8,elec234,elec567] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq);

%% get z-mean for recording 

        % get mean across the entire recording for z-scoring purposes. 
        elec_avg = @(elec)(mean(mean(elec,3),1));% fxn handle 

        z_mean_e1 =  elec_avg(elec1); 
        z_mean_e25 = elec_avg(elec25); 
        z_mean_e36 = elec_avg(elec36); 
        z_mean_e47 = elec_avg(elec47); 
        z_mean_e8 = elec_avg(elec8); 
        z_mean_e234 = elec_avg(elec234);
        z_mean_e567 = elec_avg(elec567); 
        all_z_e = vertcat(z_mean_e1,z_mean_e25,z_mean_e36,z_mean_e47,z_mean_e8,...
            z_mean_e234,z_mean_e567); 
        z_mean = mean(all_z_e,1);

        %get std 
        elec_std = @(elec)(std(std(elec,0,3),1));

        z_std_e1 = elec_std(elec1); 
        z_std_e25 = elec_std(elec25); 
        z_std_e36 = elec_std(elec36); 
        z_std_e47 = elec_std(elec47); 
        z_std_e8 = elec_std(elec8); 
        z_std_e234 = elec_std(elec234); 
        z_std_e567 = elec_std(elec567); 
        all_z_std_e = vertcat(z_std_e1,z_std_e25,z_std_e36,z_std_e47,z_std_e8,...
            z_std_e234,z_std_e567);
        z_std = std(all_z_std_e,0,1); 

%% run through analysis for each configuration 
tic
tbl_allch_alldata_stimpre_theta = []; 
tbl_allch_alldata_postpre_theta = []; 

tbl_allch_alldata_stimpre_alpha = []; 
tbl_allch_alldata_postpre_alpha = []; 

tbl_allch_alldata_stimpre_beta = []; 
tbl_allch_alldata_postpre_beta = []; 

tbl_allch_alldata_stimpre_delta = []; 
tbl_allch_alldata_postpre_delta = []; 

num_contact_configs = 7; 
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
        case 'elec234'
            data = elec234;
        case 'elec567'
            data = elec567;
    end 
    disp(all_config_names{i}); 
    names{i} = all_config_names{i}; 

    %add theta beta delta alpha to output 
[metadata,tbl_allch_stimpre_theta, tbl_allch_postpre_theta,...
    tbl_allch_stimpre_alpha, tbl_allch_postpre_alpha,...
    tbl_allch_stimpre_beta, tbl_allch_postpre_beta,...
    tbl_allch_stimpre_delta, tbl_allch_postpre_delta] = generate_heatmap_mtspectrumc_power(PatientID,metadata,data,z_mean,z_std,...
    ROI_labels,contact_config,ch_labels); 

%theta
% tbl_all_data_stimpre_theta = [tbl_all_data_stimpre_theta;{tbl_ROI_stimpre_theta}]; 
% tbl_all_data_postpre_theta = [tbl_all_data_postpre_theta;{tbl_ROI_postpre_theta}];
%theta 
tbl_allch_alldata_stimpre_theta = [tbl_allch_alldata_stimpre_theta;{tbl_allch_stimpre_theta}];
tbl_allch_alldata_postpre_theta = [tbl_allch_alldata_postpre_theta;{tbl_allch_postpre_theta}];
%beta 
% tbl_all_data_stimpre_beta = [tbl_all_data_stimpre_beta;{tbl_ROI_stimpre_beta}]; 
% tbl_all_data_postpre_beta = [tbl_all_data_postpre_beta;{tbl_ROI_postpre_beta}];
tbl_allch_alldata_stimpre_beta = [tbl_allch_alldata_stimpre_beta;{tbl_allch_stimpre_beta}];
tbl_allch_alldata_postpre_beta = [tbl_allch_alldata_postpre_beta;{tbl_allch_postpre_beta}];

%alpha 
% tbl_all_data_stimpre_alpha = [tbl_all_data_stimpre_alpha;{tbl_ROI_stimpre_alpha}]; 
% tbl_all_data_postpre_alpha = [tbl_all_data_postpre_alpha;{tbl_ROI_postpre_alpha}];
tbl_allch_alldata_stimpre_alpha = [tbl_allch_alldata_stimpre_alpha;{tbl_allch_stimpre_alpha}];
tbl_allch_alldata_postpre_alpha = [tbl_allch_alldata_postpre_alpha;{tbl_allch_postpre_alpha}];



%delta 
% tbl_all_data_stimpre_delta = [tbl_all_data_stimpre_delta;{tbl_ROI_stimpre_delta}]; 
% tbl_all_data_postpre_delta = [tbl_all_data_postpre_delta;{tbl_ROI_postpre_delta}];
tbl_allch_alldata_stimpre_delta = [tbl_allch_alldata_stimpre_delta;{tbl_allch_stimpre_delta}];
tbl_allch_alldata_postpre_delta = [tbl_allch_alldata_postpre_delta;{tbl_allch_postpre_delta}];

end 
toc

%% Combine table into one this is so stupid i hate matlab 
% allch_theta = join(tbl_allch_data_stimpre_theta{1,1},tbl_allch_data_stimpre_theta{2,1});
% tmp = join(allch_theta,tbl_allch_data_stimpre_theta{3,1}); 
% tmp = join(tmp,tbl_allch_data_stimpre_theta{4,1});
% tmp = join(tmp,tbl_allch_data_stimpre_theta{5,1});
% tmp = join(tmp,tbl_allch_data_stimpre_theta{6,1});
% tmp = join(tmp,tbl_allch_data_stimpre_theta{7,1});

% add indiv trial as output and save file 
%% Save data
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/data for Logit',PatientID);
if ~exist(outputdir)
    mkdir(outputdir) 
end   
filename = sprintf('%s/AllCh_allFOI_%s.mat',outputdir,experiment_name);
save(filename,'metadata','tbl_allch_alldata_stimpre_theta','tbl_allch_alldata_postpre_theta',...
    'tbl_allch_alldata_stimpre_alpha','tbl_allch_alldata_postpre_alpha',...
    'tbl_allch_alldata_stimpre_beta','tbl_allch_alldata_postpre_beta',...
    'tbl_allch_alldata_stimpre_delta','tbl_allch_alldata_postpre_delta','timeseriesfile','ch_labels','z_mean',...
    'z_std','names');



%% 
% %% Generate heatmaps 
%  
% matrix_for_heatmap = horzcat(tmp.stim_pre_elec1,tmp.stim_pre_elec25,...
%     tmp.stim_pre_elec234,tmp.stim_pre_elec36,tmp.stim_pre_elec47,...
%     tmp.stim_pre_elec567,tmp.stim_pre_elec8); 
% 
% % matrix_for_heatmap(14,:) = []; %******
% % ROI_labels(:,14) = []; %*****
% %% 
% figure('Position',[500 300 800 600]) 
% h = heatmap(matrix_for_heatmap); 
% 
% colormap redbluecmap ; 
% %co.YDisplayLabels = 'z-scored power: stim minus prestim';
% h.YDisplayLabels = char(ROI_labels); 
% h.XDisplayLabels = char(all_config_names); 
% h.ColorLimits = [-3 3]; 
% title(sprintf('stim minus pre %s %s',FOI,experiment_name));
% 
% filename = sprintf('%s_%s_%s,stimminuspre1.png',PatientID,experiment_name,FOI); 
% saveas(gcf, filename)
% % 
% % %% group data into table ?
% % 
% %% post minus pre 
% clear tmp ; clear ROI_tbl_all_configs; 
% ROI_tbl_all_configs = join(tbl_all_data_postpre{1,1},tbl_all_data_postpre{2,1});
% tmp = join(ROI_tbl_all_configs,tbl_all_data_postpre{3,1}); 
% tmp = join(tmp,tbl_all_data_postpre{4,1});
% tmp = join(tmp,tbl_all_data_postpre{5,1});
% tmp = join(tmp,tbl_all_data_postpre{6,1});
% tmp = join(tmp,tbl_all_data_postpre{7,1});
% 
% %% Generate heatmaps 
%  
% matrix_for_heatmap_postpre = horzcat(tmp.post_pre_elec1,tmp.post_pre_elec25,...
%     tmp.post_pre_elec234,tmp.post_pre_elec36,tmp.post_pre_elec47,...
%     tmp.post_pre_elec567,tmp.post_pre_elec8); 
% 
% % matrix_for_heatmap_postpre(14,:) = []; 
% %% 
% 
% figure('Position',[500 300 800 600]) 
% h = heatmap(matrix_for_heatmap_postpre); 
% 
% colormap redbluecmap ; 
% h.YDisplayLabels = char(ROI_labels); 
% h.XDisplayLabels = char(all_config_names); 
% h.ColorLimits = [-3 3]; 
% title(sprintf('post minus pre %s %s',FOI,experiment_name));
% 
% filename = sprintf('%s_%s_%s,postminuspre.png',PatientID,experiment_name,FOI); 
% saveas(gcf, filename)
% 
% %% Save data 
% 
% outputdir = sprintf('E:\DBSTRD\%s\Experiments\15s_stim\Processed Data\',PatientID); 
% filename = sprintf('%s\ROI_and_allch_%s_%s.mat',outputdir,experiment_name,FOI); 
% save(filename,'metadata','ROI_labels','z_mean','z_std','tbl_all_data_stimpre','tbl_all_data_postpre',...
%     'tbl_allch_alldata_stimpre','tbl_allch_alldata_postpre');
% 

% %% Run stats? 