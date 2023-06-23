%% 

% This script relies on the run_setup_15s script,
% generate_heatmap_mtspectrumc_power script 
%% setup input parameters & run setup file 
clear all 
close all
PatientID = 'DBSTRD002'; DBS_target = 'SCC'; hemi = 'r'; stim_freq = 130; 
FOI = input('enter FOI ','s'); 
[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names, upper_freq,lower_freq,metadata,ROI_labels,num_contact_configs,...
    experiment_name,elec1,elec25,elec36,elec47,...
    elec8,elec234,elec567] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, FOI); 

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

ROI_stim_minus_pre = {}; 
ROI_post_minus_pre = {}; 
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

[ROI_labels,ROI_stim_minus_pre{i},...
            ROI_post_minus_pre{i}] = generate_mtspectrumc_power_withstats(PatientID,metadata,data,z_mean,z_std,...
     lower_freq,upper_freq, ROI_labels,contact_config); 
ROI_stim_minus_pre{i} = [ROI_stim_minus_pre{i},{contact_config}];
ROI_post_minus_pre{i} = [ROI_post_minus_pre{i},{contact_config}]; 

 
end 

toc

%% run statistics 
%save data 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/',PatientID); 
filename = sprintf('%s/%s_%s_%s_ROIforstats.mat',...
    outputdir,PatientID,experiment_name,FOI); 
save(filename,'ROI_labels','ROI_stim_minus_pre','ROI_post_minus_pre',...
    'timeseriesfile','FOI')





