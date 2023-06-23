%% This script is used to generate .mat files to plot bar plots overlaid with scatter plots in seaborn for Fig 2c. 
% showing ON vs OFF differences 

% output: tbl with PSD value, FOI, condition, with ROI and patientID in
% name. 

%% load data. 
clear 
clc 
close 
% PatientID = 'DBSTRD002'; 
% 
% switch PatientID 
%     case 'DBSTRD001'
%         load(sprintf('%s_Jan-03-2023_allroi_tstat_ch.mat', PatientID)); 
%     case 'DBSTRD002'
%         load(sprintf('%s_06-Jan-2023_allroi_tstat_ch.mat', PatientID)); 
% end 

%% get all labels. 

size(ch_data) 
stim_labels = stim_tbl.Stim_Labels; 
roi_label = ch_labels.ROI; 

%% get out poststim data 
DBS_target = 'rSCC'; 
poststim_idx = find(contains(stim_tbl.Stim_Labels,'poststim')==1 & contains(stim_tbl.DBS_target,DBS_target)==1); 
poststim_data = ch_data(poststim_idx,:,:);

%% get out FOI data. 

foi_idx = 6; 
foi =  {'delta';'theta';'alpha';'beta';'lowg';'highg'}; 

poststim_foi = poststim_data(:,:,foi_idx);

%% get out ROI data. 

stim_foi_acc = mean(poststim_foi(:,acc_idx),2);
stim_foi_amy = mean(poststim_foi(:,amy_idx),2);
stim_foi_dpf = mean(poststim_foi(:,dpf_idx),2);
stim_foi_mof = mean(poststim_foi(:,mof_idx),2);
stim_foi_lof = mean(poststim_foi(:,lof_idx),2);
stim_foi_vpf = mean(poststim_foi(:,vpf_idx),2);

acc_label = string(repmat('ACC',length(stim_foi_acc),1));
amy_label = string(repmat('AMY',length(stim_foi_amy),1));
dpf_label = string(repmat('DPF',length(stim_foi_dpf),1));
mof_label= string(repmat('MOF',length(stim_foi_mof),1));
lof_label= string(repmat('LOF',length(stim_foi_lof),1));
vpf_label= string(repmat('VPF',length(stim_foi_vpf),1));

% stim labels. 
poststim_label = string(repmat('poststim',length(stim_foi_acc),1));

% foi labels. 
foi_poststim_label = string(repmat(foi{foi_idx},length(stim_foi_acc),1));

%% repeat for baseline data. 

bl_idx = find(contains(stim_tbl.Stim_Labels,'Baseline')==1); 
bl_data = ch_data(bl_idx,:,:);

%% get out FOI data. 

bl_foi = bl_data(:,:,foi_idx);

%% get out ROI data. 

bl_foi_acc = mean(bl_foi(:,acc_idx),2);
bl_foi_amy = mean(bl_foi(:,amy_idx),2);
bl_foi_dpf = mean(bl_foi(:,dpf_idx),2);
bl_foi_mof = mean(bl_foi(:,mof_idx),2);
bl_foi_lof = mean(bl_foi(:,lof_idx),2);
bl_foi_vpf = mean(bl_foi(:,vpf_idx),2);

acc_label = string(repmat('ACC',length(bl_foi_acc),1));
amy_label = string(repmat('AMY',length(bl_foi_amy),1));
dpf_label = string(repmat('DPF',length(bl_foi_dpf),1));
mof_label= string(repmat('MOF',length(bl_foi_mof),1));
lof_label= string(repmat('LOF',length(bl_foi_lof),1));
vpf_label= string(repmat('VPF',length(bl_foi_vpf),1));

% bl labels. 
bl_label = string(repmat('baseline',length(bl_foi_acc),1));

% foi labels. 
foi_bl_label = string(repmat(foi{foi_idx},length(bl_foi_acc),1));

%% concatenate baseline and poststim. 

acc = [bl_foi_acc;stim_foi_acc]; 
amy = [bl_foi_amy;stim_foi_amy]; 
dpf = [bl_foi_dpf;stim_foi_dpf]; 
mof = [bl_foi_mof;stim_foi_mof]; 
lof = [bl_foi_lof;stim_foi_lof]; 
vpf = [bl_foi_vpf;stim_foi_vpf]; 

condition = [bl_label;poststim_label]; 
foi = [foi_bl_label;foi_poststim_label];

%% 

% save ACC 
save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/acc_%s_%s_%s.mat',...
    PatientID, foi{foi_idx}, DBS_target), 'acc','condition','foi');  
% save AMY
save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/amy_%s_%s_%s.mat',...
    PatientID, foi{foi_idx}, DBS_target), 'amy','condition','foi'); 
% save DPF 
save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/dpf_%s_%s_%s.mat',...
    PatientID, foi{foi_idx}, DBS_target), 'dpf','condition','foi'); 
% save MOF 
save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/mof_%s_%s_%s.mat',...
    PatientID, foi{foi_idx}, DBS_target), 'mof','condition','foi'); 
% save LOF 
save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/lof_%s_%s_%s.mat',...
    PatientID, foi{foi_idx}, DBS_target), 'lof','condition','foi'); 
% save VPF 
save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/vpf_%s_%s_%s.mat',...
    PatientID, foi{foi_idx}, DBS_target), 'vpf','condition','foi'); 
