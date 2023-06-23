
%% Get data (CODE FROM DBSTARGET_HYPOTHESISTEST.M) 


clc 
close all 
clear
PatientID = 'DBSTRD002'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
%data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 

%Load channels labels that will be used later. 
ROI_labels = load(sprintf('%s_ROI_labels.mat',PatientID));
ch_labels = ROI_labels.good_ch_labels;
label_tbl = ROI_labels.summary_ch_info;
%data(:,ROI_labels.goodch_idx_remove,:) = [];

%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
comp_labels = comp_labels.comps;

switch PatientID
    case 'DBSTRD001'
        stim_idx_lSCC = find(comp_labels{1, 9}.label_vector==1);
        stim_idx_rSCC = find(comp_labels{1, 17}.label_vector == 1);
        stim_idx_lVCVS = find(comp_labels{1, 31}.label_vector == 1);
        stim_idx_rVCVS = find(comp_labels{1,37}.label_vector == 1);
        baseline_idx = find(comp_labels{1, 1}.label_vector==0);        
        % Get vector of labels
        
        labels = {};
        labels{1} = comp_labels{9}.label_vector;
        labels{2} = comp_labels{17}.label_vector;
        labels{3} = comp_labels{31}.label_vector;
        labels{4} = comp_labels{37}.label_vector;
    case 'DBSTRD002'
        
        stim_idx_lSCC = find(comp_labels{1, 9}.label_vector==1);
        stim_idx_rSCC = find(comp_labels{1, 17}.label_vector == 1);
        stim_idx_lVCVS = find(comp_labels{1, 33}.label_vector == 1);
        stim_idx_rVCVS = find(comp_labels{1,41}.label_vector == 1); 
        baseline_idx = find(comp_labels{1,9}.label_vector == 0); 
        labels = {};
        labels{1} = comp_labels{9}.label_vector;
        labels{2} = comp_labels{17}.label_vector;
        labels{3} = comp_labels{33}.label_vector;
        labels{4} = comp_labels{41}.label_vector;
          
    case 'DBSTRD003'
        
        stim_idx_lSCC = find(comp_labels{1, 4}.label_vector==1);
        stim_idx_rSCC = find(comp_labels{1, 7}.label_vector == 1);
        stim_idx_lVCVS = find(comp_labels{1, 13}.label_vector == 1);
        stim_idx_rVCVS = find(comp_labels{1, 16}.label_vector == 1); 
        baseline_idx = find(comp_labels{1, 4}.label_vector == 0); 
        labels = {};
        labels{1} = comp_labels{4}.label_vector;
        labels{2} = comp_labels{7}.label_vector;
        labels{3} = comp_labels{13}.label_vector;
        labels{4} = comp_labels{16}.label_vector;
        
end
%% Get data that excludes stim (artifact) 

for i = 1:num_trials
    
    is_stim = contains(tbl.Stim_Labels(i), 'stim1') == 1 || contains(tbl.Stim_Labels(i), 'stim2') == 1 ||contains(tbl.Stim_Labels(i), 'stim3') == 1;
    if is_stim == 1
        nostim_idx(i) = 0;
    elseif is_stim == 0
        nostim_idx(i) = 1;
    end
end
data_nostim = data(logical(nostim_idx'),:,:);

idx = logical(nostim_idx'); 

new_tbl = tbl(idx,:); 
    
%% relabel data after stim trials thrown out from zscoring 

    for i = 1:numel(labels)
        labels{i} = labels{i}(logical(nostim_idx'));  
    end 

%% Calculate z-score 

data_zscore = zeros(size(data_nostim)); 
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - nanmean(data_nostim(:,j,k)))./nanstd(data_nostim(:,j,k)); 
        end 
    end 
end 


%% Define which type of data, Channel or ROI 

data_type_list = {'ROI', 'Ch'}; 
% [idx_data_list,tf_data] = listdlg('ListString', data_type_list, 'ListSize',[150,250]); 
% data_list_config = data_type_list{idx_data_list}; 
% disp(data_list_config)
% data_dim = data_list_config; 
 
data_dim = 'ROI'; 
if strcmp(PatientID,'DBSTRD003')==1
    %load('DBSTRD003_referencedsig_SCC_f130_234.mat','metadata') %** CHECK
    %THIS?? 
end 
if strcmp(data_dim,'ROI')==1 
    % get data that will be input into permutations fxn. 
    [matrix_ROI,ROI_labels1] = generate_ROI_from_ch(data_zscore,PatientID,2);
end 

%% Subtract from baseline. 

baseline_data = matrix_ROI(labels{1}==0,:,:); 
baseline_cond = new_tbl(labels{1}==0,:); 
mean_bl = mean(baseline_data,1); 

% SCC. 
lSCC_data = matrix_ROI(labels{1}==1,:,:); 
rSCC_data = matrix_ROI(labels{2}==1,:,:); 
%SCC conds labels. 
lSCC_cond = new_tbl(labels{1}==1,:); 
rSCC_cond = new_tbl(labels{2}==1,:); 

% VCVS. 
lVCVS_data = matrix_ROI(labels{3}==1,:,:); 
rVCVS_data = matrix_ROI(labels{4}==1,:,:); 
% VCVS conds labels. 
lVCVS_cond = new_tbl(labels{3} ==1,:); 
rVCVS_cond = new_tbl(labels{4} ==1,:); 

%% Get subtracted power values. 

diff_lSCC = lSCC_data - mean_bl; 
diff_rSCC = rSCC_data - mean_bl;

diff_lVCVS = lVCVS_data - mean_bl;
diff_rVCVS = rVCVS_data - mean_bl; 

%% New matrix with just SCC and VCVS (for comparison) . 

% get number of trials for each DBS target. 
num_lSCC_trials = size(diff_lSCC,1); 
num_lVCVS_trials = size(diff_lVCVS,1); 
num_rSCC_trials = size(diff_rSCC,1); 
num_rVCVS_trials = size(diff_rVCVS,1); 

%% choose roi and foi of interest 

%lSCC 
for i = 1:size(diff_lSCC,2)
    for j = 1:size(diff_lSCC,3)
        y_diff_lSCC = diff_lSCC(:,i,j);
        group_diff_lSCC = lSCC_cond.Elec_label;
        [p_diff_lSCC(i,j),~,stats_diff_lSCC] = anova1(y_diff_lSCC,group_diff_lSCC);
        close
        close all hidden
        means_lSCC{i,j} = stats_diff_lSCC.means;      
        [c_diff_lSCC,~,h_diff_lSCC,gnames_diff_lSCC] = multcompare(stats_diff_lSCC);
    end
end
    




%% Run anova for FOI and ROI of interest 

%lVCVS theta OFC 
lVCVS_new = diff_lVCVS(:,2,2); 
%rOFC theta OFC 
rVCVS_new = diff_rVCVS(:,2,2);
%lSCC high gamma amygdala 
lSCC_new = diff_lSCC(:,5,6);
%rSCC high gamma amygdala 
rSCC_new = diff_rSCC(:,5,6);


%anova 
[p_lSCC,~,stats_diff_lSCC] = anova1(lSCC_new,group_diff_lSCC);
[p_rSCC,~,stats_diff_rSCC] = anova1(rSCC_new,group_diff_rSCC);
[p_lVCVS,~,stats_diff_lVCVS] = anova1(lVCVS_new,group_diff_lVCVS);
[p_rVCVS,~,stats_diff_rVCVS] = anova1(rVCVS_new,group_diff_rVCVS);
close all 

%%
[c_lSCC,~,h_lSCC,gnames_lSCC] = multcompare(stats_diff_lSCC);
[c_rSCC,~,h_rSCC,gnames_rSCC] = multcompare(stats_diff_rSCC);
[c_lVCVS,~,h_lVCVS,gnames_lVCVS] = multcompare(stats_diff_lVCVS);
[c_rVCVS,~,h_rVCVS,gnames_rVCVS] = multcompare(stats_diff_rVCVS);



    %% rSCC
    
    for i = 1:size(diff_rSCC,2)
        for j = 1:size(diff_rSCC,3)
            y_diff_rSCC = diff_rSCC(:,i,j);
            group_diff_rSCC = rSCC_cond.Elec_label;
            [p_diff_rSCC(i,j),~,stats_diff_rSCC] = anova1(y_diff_rSCC,group_diff_rSCC);
            close
            close all hidden
            means_rSCC{i,j} = stats_diff_rSCC.means;
            [c_diff_rSCC,~,h_diff_rSCC,gnames_diff_rSCC] = multcompare(stats_diff_rSCC);
        end
    end

    
   
    %% lVCVS 
    for i = 1:size(diff_lVCVS,2)
        for j = 1:size(diff_lVCVS,3)
            y_diff_lVCVS = diff_lVCVS(:,i,j);
            group_diff_lVCVS = lVCVS_cond.Elec_label;
            [p_diff_lVCVS(i,j),~,stats_diff_lVCVS] = anova1(y_diff_lVCVS,group_diff_lVCVS);
            close
            close all hidden
            means_lVCVS{i,j} = stats_diff_lVCVS.means;
            [c_diff_lVCVS,~,h_diff_lVCVS,gnames_diff_lVCVS] = multcompare(stats_diff_lVCVS);
        end
    end
    
   
    %% rVCVS
    for i = 1:size(diff_rVCVS,2)
        for j = 1:size(diff_rVCVS,3)
            y_diff_rVCVS = diff_rVCVS(:,i,j);
            group_diff_rVCVS = rVCVS_cond.Elec_label;
            [p_diff_rVCVS(i,j),~,stats_diff_rVCVS] = anova1(y_diff_rVCVS,group_diff_rVCVS);
            close
            close all hidden
            means_rVCVS{i,j} = stats_diff_rVCVS.means;
            [c_diff_rVCVS,~,h_diff_rVCVS,gnames_diff_rVCVS] = multcompare(stats_diff_rVCVS);
        end
    end
    

%% Make table containing all labels and PSD data. 


roi_n =  1; 
foi_n = 6; 
clear lSCC_tbl rSCC_tbl lVCVS_tbl rVCVS_tbl
%lSCC
lSCC_tbl = lSCC_cond; 
lSCC_tbl{:,5} = diff_lSCC(:,roi_n,foi_n); 
lSCC_tbl.Properties.VariableNames{5} = sprintf('PSD_%s',data_file.freq_cond{foi_n}); 
%rSCC
rSCC_tbl = rSCC_cond; 
rSCC_tbl{:,5} = diff_rSCC(:,roi_n,foi_n); 
rSCC_tbl.Properties.VariableNames{5} = sprintf('PSD_%s',data_file.freq_cond{foi_n}); 
%lVCVS
lVCVS_tbl = lVCVS_cond; 
lVCVS_tbl{:,5} = diff_lVCVS(:,roi_n,foi_n); 
lVCVS_tbl.Properties.VariableNames{5} = sprintf('PSD_%s',data_file.freq_cond{foi_n}); 
%rVCVS
rVCVS_tbl = rVCVS_cond; 
rVCVS_tbl{:,5} = diff_rVCVS(:,roi_n,foi_n); 
rVCVS_tbl.Properties.VariableNames{5} = sprintf('PSD_%s',data_file.freq_cond{foi_n}); 

%% Save with label of ROI, FOI, DBS target (left or right) for sb viz. 

%SCC
writetable(lSCC_tbl,sprintf('/Users/anushaallawala/Python_projects/PSD_%s_foi%d_roi%d_currstr_lSCC.csv',PatientID,foi_n,roi_n)); 
writetable(rSCC_tbl,sprintf('/Users/anushaallawala/Python_projects/PSD_%s_foi%d_roi%d_currstr_rSCC.csv',PatientID,foi_n,roi_n));
%VCVS 
writetable(lVCVS_tbl,sprintf('/Users/anushaallawala/Python_projects/PSD_%s_foi%d_roi%d_currstr_lVCVS.csv',PatientID,foi_n,roi_n)); 
writetable(rVCVS_tbl,sprintf('/Users/anushaallawala/Python_projects/PSD_%s_foi%d_roi%d_currstr_rVCVS.csv',PatientID,foi_n,roi_n));


%% Get avgs for ROI across trials conditions. 
foi_idx = foi_n; %high gamma 
% get FOI of interest. 

lSCC_all_ROI = cat(1,means_lSCC{:,foi_idx});
rSCC_all_ROI = cat(1,means_rSCC{:,foi_idx});
rVCVS_all_ROI = cat(1,means_rVCVS{:,foi_idx}); 
lVCVS_all_ROI = cat(1,means_lVCVS{:,foi_idx}); 

%% for 002 

%lSCC 
new_lSCC_all_ROI(1,:) = lSCC_all_ROI(4,:); 
new_lSCC_all_ROI(2,:) = lSCC_all_ROI(5,:); 
new_lSCC_all_ROI(3,:) = lSCC_all_ROI(6,:); 
new_lSCC_all_ROI(4,:) = lSCC_all_ROI(2,:);
new_lSCC_all_ROI(5,:) = lSCC_all_ROI(1,:);
new_lSCC_all_ROI(6,:) = lSCC_all_ROI(7,:); 
new_lSCC_all_ROI(7,:) = lSCC_all_ROI(8,:); 
new_lSCC_all_ROI(8,:) = lSCC_all_ROI(3,:);

%lVCVS 
new_lVCVS_all_ROI(1,:) = lVCVS_all_ROI(4,:); 
new_lVCVS_all_ROI(2,:) = lVCVS_all_ROI(5,:); 
new_lVCVS_all_ROI(3,:) = lVCVS_all_ROI(6,:); 
new_lVCVS_all_ROI(4,:) = lVCVS_all_ROI(2,:);
new_lVCVS_all_ROI(5,:) = lVCVS_all_ROI(1,:);
new_lVCVS_all_ROI(6,:) = lVCVS_all_ROI(7,:); 
new_lVCVS_all_ROI(7,:) = lVCVS_all_ROI(8,:); 
new_lVCVS_all_ROI(8,:) = lVCVS_all_ROI(3,:);

%rSCC
new_rSCC_all_ROI(1,:) = rSCC_all_ROI(4,:); 
new_rSCC_all_ROI(2,:) = rSCC_all_ROI(5,:); 
new_rSCC_all_ROI(3,:) = rSCC_all_ROI(6,:); 
new_rSCC_all_ROI(4,:) = rSCC_all_ROI(2,:);
new_rSCC_all_ROI(5,:) = rSCC_all_ROI(1,:);
new_rSCC_all_ROI(6,:) = rSCC_all_ROI(7,:); 
new_rSCC_all_ROI(7,:) = rSCC_all_ROI(8,:); 
new_rSCC_all_ROI(8,:) = rSCC_all_ROI(3,:);

%rVCVS
new_rVCVS_all_ROI(1,:) = rVCVS_all_ROI(4,:); 
new_rVCVS_all_ROI(2,:) = rVCVS_all_ROI(5,:); 
new_rVCVS_all_ROI(3,:) = rVCVS_all_ROI(6,:); 
new_rVCVS_all_ROI(4,:) = rVCVS_all_ROI(2,:);
new_rVCVS_all_ROI(5,:) = rVCVS_all_ROI(1,:);
new_rVCVS_all_ROI(6,:) = rVCVS_all_ROI(7,:); 
new_rVCVS_all_ROI(7,:) = rVCVS_all_ROI(8,:); 
new_rVCVS_all_ROI(8,:) = rVCVS_all_ROI(3,:);

ROI_labels_new = {'ACC','AMY','DPF','LOF','MOF','STG','SFG','VPF'}; 
%% 002 MAKE HEATMAP OF AVG VALUES 

addpath(genpath('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/')); 
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/customcolormap/')); 
%make colormap. 
mycolormap = customcolormap(linspace(0,1,9), {'#c3553a', '#ce7963', '#dba598',...
    '#e3bfb6','#f2f1f1', '#b9cfd5','#9cbcc6','#71a0ae','#3f7f93'}); 

%lSCC
%save 
figure() 
f = heatmap(ROI_labels_new,gnames_diff_lSCC,new_lSCC_all_ROI', 'FontName', 'Helvetica Neue')
f.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_lSCC_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 


%rSCC
%save 
figure() 
f1 = heatmap(ROI_labels_new,gnames_diff_rSCC,new_rSCC_all_ROI','FontName', 'Helvetica Neue')
f1.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])

%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_rSCC_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 


%lVCVS 
figure() 
f2 = heatmap(ROI_labels_new,gnames_diff_lVCVS,new_lVCVS_all_ROI','FontName', 'Helvetica Neue')
f2.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])

%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_lVCVS_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 


%rVCVS
figure() 
f3 = heatmap(ROI_labels_new,gnames_diff_rVCVS,new_rVCVS_all_ROI','FontName', 'Helvetica Neue')
f3.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])

%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_rVCVS_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 





%% make heatmap of avg values 

addpath(genpath('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/')); 
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/customcolormap/')); 
%make colormap. 
mycolormap = customcolormap(linspace(0,1,9), {'#c3553a', '#ce7963', '#dba598',...
    '#e3bfb6','#f2f1f1', '#b9cfd5','#9cbcc6','#71a0ae','#3f7f93'}); 

%lSCC
%save 
figure() 
f = heatmap(ROI_labels1,gnames_diff_lSCC,lSCC_all_ROI', 'FontName', 'Helvetica Neue')
f.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_lSCC_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 


%rSCC
%save 
figure() 
f1 = heatmap(ROI_labels1,gnames_diff_rSCC,rSCC_all_ROI','FontName', 'Helvetica Neue')
f1.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_rSCC_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 


%lVCVS 
figure() 
f2 = heatmap(ROI_labels1,gnames_diff_lVCVS,lVCVS_all_ROI','FontName', 'Helvetica Neue')
f2.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_lVCVS_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 


%rVCVS
figure() 
f3 = heatmap(ROI_labels1,gnames_diff_rVCVS,rVCVS_all_ROI','FontName', 'Helvetica Neue')
f3.CellLabelFormat = '%.2f';
colormap(mycolormap);
caxis([-0.8 0.8])
%save 
saveas(gcf, sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Current_steer/%s_%s_rVCVS_%s_%s.svg',PatientID,data_file.freq_cond{foi_idx}, date, ROI_labels1{roi_n})); 

%% make heatmap of FOI for all DBS targets. 



%tmp_mean(i,j) = varfun(@mean,lSCC_cond(:,[3,5]),'GroupingVariables',{'Elec_label'} );



