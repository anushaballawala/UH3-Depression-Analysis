%% Load data. 

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

% adjust caxis values and re-run for 001 and 002 
%% clean up and organize data. 

clearvars -except PatientID ch_data ch_labels stim_labels stim_tbl data_file t_stat tbl ROI_info comp_labels acc_idx amy_idx lof_idx mof_idx vpf_idx
 
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

%% Compute variance for each DBS lead, foi and roi. 

num_ROI = size(all_roi_data,2); 
num_foi = size(all_roi_data,3); 

for i = 1:num_ROI 
    for j = 1:num_foi 
        
        var_all_lSCC(i,j) = var(lSCC_data(:,i,j)); 
        var_all_rSCC(i,j) = var(rSCC_data(:,i,j)); 
        var_all_lVCVS(i,j) = var(lVCVS_data(:,i,j)); 
        var_all_rVCVS(i,j) = var(rVCVS_data(:,i,j)); 

    end 
end 

% find highest variance

% for 001 
% beta LOF for lSCC 

% beta mof for lVCVS or beta LOF 

% mof alpha for rSCC or mof theta 

% mof delta 

%% Get stim labels for each DBS lead. 

lSCC_cond = stim_tbl(lSCC_idx,:); 
rSCC_cond = stim_tbl(rSCC_idx,:); 
lVCVS_cond = stim_tbl(lVCVS_idx,:); 
rVCVS_cond = stim_tbl(rVCVS_idx,:); 

%% Get current steering labels.

group_diff_lSCC = lSCC_cond.Elec_label;
group_diff_rSCC = rSCC_cond.Elec_label; 
group_diff_lVCVS = lVCVS_cond.Elec_label; 
group_diff_rVCVS = rVCVS_cond.Elec_label; 

%% Compute anova. 

% lSCC
for i = 1:size(lSCC_data,2)
    for j = 1:size(lSCC_data,3)
        y_lSCC = lSCC_data(:,i,j);
        [p_lSCC(i,j),tbl_lSCC{i,j},stats_lSCC] = anova1(y_lSCC,group_diff_lSCC);
        close
        close all hidden
        means_lSCC{i,j} = stats_lSCC.means;      
        [c_lSCC,~,h_lSCC,gnames_lSCC] = multcompare(stats_lSCC);
    end
end

%rSCC

for i = 1:size(rSCC_data,2)
    for j = 1:size(rSCC_data,3)
        y_rSCC = rSCC_data(:,i,j);
        [p_rSCC(i,j),tbl_rSCC{i,j},stats_rSCC] = anova1(y_rSCC,group_diff_rSCC);
        close
        close all hidden
        means_rSCC{i,j} = stats_rSCC.means;      
        [c_rSCC,~,h_rSCC,gnames_rSCC] = multcompare(stats_rSCC);
    end
end

%lVCVS 
for i = 1:size(lVCVS_data,2)
    for j = 1:size(lVCVS_data,3)
        y_lVCVS = lVCVS_data(:,i,j);
        [p_lVCVS(i,j),tbl_lVCVS{i,j},stats_lVCVS] = anova1(y_lVCVS,group_diff_lVCVS);
        close
        close all hidden
        means_lVCVS{i,j} = stats_lVCVS.means;      
        [c_lVCVS,~,h_lVCVS,gnames_lVCVS] = multcompare(stats_lVCVS);
    end
end


%rVCVS 
for i = 1:size(rVCVS_data,2)
    for j = 1:size(rVCVS_data,3)
        y_rVCVS = rVCVS_data(:,i,j);
        [p_rVCVS(i,j),tbl_rVCVS{i,j},stats_rVCVS] = anova1(y_rVCVS,group_diff_rVCVS);
        close
        close all hidden
        means_rVCVS{i,j} = stats_rVCVS.means;      
        [c_rVCVS,~,h_rVCVS,gnames_rVCVS] = multcompare(stats_rVCVS);
    end
end

%% save 
clear f h_lSCC h_lVCVS h_rSCC h_rVCVS

filename = sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/current_steer/%s_%s_currsteer.mat',...
    PatientID,date); 

save(filename); 
%% Make heatmap 

addpath(genpath('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/')); 
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/customcolormap/')); 
%make colormap. 
% mycolormap = customcolormap(linspace(0,1,9), {'#c3553a', '#ce7963', '#dba598',...
%     '#e3bfb6','#f2f1f1', '#b9cfd5','#9cbcc6','#71a0ae','#3f7f93'}); 

mycolormap = customcolormap(linspace(0,1,9), {'#a9383c', '#b14a4b', '#ba5e5d',...
    '#d39795','#faf5f4', '#d5d9e2','#9eb0cb','#6389bd','#2d6dbc'}); 

%foi_n = 2; 
lSCC_all_ROI = {}; rSCC_all_ROI = {}; lVCVS_all_ROI = {}; rVCVS_all_ROI = {}; 
for k = 1:num_foi
    lSCC_all_ROI{k} = cat(1,means_lSCC{:,k});
    rSCC_all_ROI{k} = cat(1,means_rSCC{:,k});
    lVCVS_all_ROI{k} = cat(1,means_lVCVS{:,k});
    rVCVS_all_ROI{k} = cat(1,means_rVCVS{:,k});
end 

for k = 1:num_foi 
%lSCC
figure() 
f = heatmap(roi_labels_new,gnames_lSCC,lSCC_all_ROI{k}', 'FontName', 'Helvetica Neue');
f.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/Current Steer/%s_%s_lSCC_%s.svg',PatientID,foi_labels_new{k}, date)); 
title(sprintf('lSCC %s',foi_labels_new{k}))

%rSCC
figure() 
f = heatmap(roi_labels_new,gnames_rSCC,rSCC_all_ROI{k}', 'FontName', 'Helvetica Neue');
f.CellLabelFormat = '%.2f';
colormap(mycolormap);

caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/Current Steer/%s_%s_rSCC_%s.svg',PatientID,foi_labels_new{k}, date)); 
title(sprintf('rSCC %s',foi_labels_new{k}))

%lVCVS 
figure() 
f = heatmap(roi_labels_new,gnames_lVCVS,lVCVS_all_ROI{k}', 'FontName', 'Helvetica Neue');
f.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/Current Steer/%s_%s_lVCVS_%s.svg',PatientID,foi_labels_new{k}, date)); 
title(sprintf('lVCVS %s',foi_labels_new{k}))


%rVCVS
figure() 
f = heatmap(roi_labels_new,gnames_rVCVS,rVCVS_all_ROI{k}', 'FontName', 'Helvetica Neue');
f.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/Frontotemporal_DBS_TRD_2022/Current Biology Draft/Figures/Current Steer/%s_%s_rVCVS_%s.svg',PatientID,foi_labels_new{k}, date)); 
title(sprintf('rVCVS %s',foi_labels_new{k}))

end 


    
%% Find highest variance. 

mean_var_lSCC = mean(var_all_lSCC,1); 
mean_var_rSCC = mean(var_all_rSCC,1); 
mean_var_lVCVS = mean(var_all_lVCVS,1);
mean_var_rVCVS = mean(var_all_rVCVS,1);

% find median 
med_var_lSCC = median(var_all_lSCC,1); 
med_var_rSCC = median(var_all_rSCC,1); 
med_var_lVCVS = median(var_all_lVCVS,1);
med_var_rVCVS = median(var_all_rVCVS,1);

%% sort data. 
[sd_lSCC,r_lSCC]=sort(mean_var_lSCC,'descend'); 
[sd_rSCC,r_rSCC]=sort(mean_var_rSCC,'descend'); 
[sd_lVCVS,r_lVCVS]=sort(mean_var_lVCVS,'descend'); 
[sd_rVCVS,r_rVCVS]=sort(mean_var_rVCVS,'descend'); 

%% 
for i = 1:5
    for j = 1:6 
        new_var_lSCC(i,j) = var(means_lSCC{i,j}); 
        new_var_rSCC(i,j) = var(means_rSCC{i,j}); 
        new_var_lVCVS(i,j) = var(means_lVCVS{i,j}); 
        new_var_rVCVS(i,j) = var(means_rVCVS{i,j});       
    end
end

mean_new_var_lSCC = mean(new_var_lSCC,1); 
mean_new_var_rSCC = mean(new_var_rSCC,1); 
mean_new_var_lVCVS = mean(new_var_lVCVS,1); 
mean_new_var_rVCVS = mean(new_var_rVCVS,1); 

%% rank frequency bands by percentage of ROIs that show the greatest variance. 


for i = 1:5
    for j  = 1:6 
        f_lSCC(i,j) = tbl_lSCC{i,j}{2, 5}; 
        f_lVCVS(i,j) = tbl_lVCVS{i,j}{2, 5};
        f_rSCC(i,j) = tbl_rSCC{i,j}{2, 5};
        f_rVCVS(i,j) = tbl_rVCVS{i,j}{2, 5}  ;
    end 
end 

foi_f_lSCC = sum(f_lSCC,1); 
foi_f_lVCVS = sum(f_lVCVS,1); 
foi_f_rSCC = sum(f_rSCC,1); 
foi_f_rVCVS = sum(f_rVCVS,1); 

%% 

perc_p = p_lSCC<0.05; 



