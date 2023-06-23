%% Make RGB VALUE FILES FOR PSD CHANGES FOR VCVS AND SCC MMVT VIZ
%clear

clc
close all
clear
PatientID = 'DBSTRD001';
addpath(genpath('/Users/anushaallawala/Data/'));

data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats.mat',PatientID)) ;
data = data_file.data;

DBS_labels = string(data_file.DBS_labels);
stim_labels = string(data_file.stim_labels);
num_trials = size(data,1);
tbl = data_file.tbl;

ROI_labels = load(sprintf('%s_ROI_labels.mat',PatientID));
ch_labels = ROI_labels.good_ch_labels;
label_tbl = ROI_labels.summary_ch_info;
data(:,ROI_labels.goodch_idx_remove,:) = [];

%% Convert descriptive labels to binary labels in a vector for each combo being tested.

comp_labels = load(sprintf('labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
comp_labels = comp_labels.comps;

%% z-score data

% calculate zscore for all values except for stim1,2,3

%get data that excludes stim (artifact)
for i = 1:num_trials
    
    is_stim = contains(tbl.Stim_Labels(i), 'stim1') == 1 || contains(tbl.Stim_Labels(i), 'stim2') == 1 ||contains(tbl.Stim_Labels(i), 'stim3') == 1;  ;
    if is_stim == 1
        nostim_idx(i) = 0;
    elseif is_stim == 0
        nostim_idx(i) = 1;
    end
end

data_nostim = data(logical(nostim_idx'),:,:);
disp('size of data_nostim matrix: trials x channels x freqbands')
size(data_nostim)

%% relabel data after stim trials thrown out from zscoring

comp_labels_old = comp_labels;

for i = 1:numel(comp_labels)
    
    comp_labels{1,i}.label_vector = comp_labels{1,i}.label_vector(logical(nostim_idx'));
    
end

labels = {};
for i = 1:numel(comp_labels)
    labels{i} = comp_labels{i}.label_vector;
end
%% Calculate z-score (across windows with no stimulation artifact).

data_zscore = zeros(size(data_nostim));
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - mean(data_nostim(:,j,k)))./std(data_nostim(:,j,k));
        end
    end
end

%% Get labels for poststim for std.

switch PatientID
    case 'DBSTRD001'
        stim_idx_lSCC = find(labels{1, 9}==1);
        stim_idx_rSCC = find(labels{1, 17} == 1);
        stim_idx_lVCVS = find(labels{1, 31} == 1);
        stim_idx_rVCVS = find(labels{1,37} == 1);
        baseline_idx = find(labels{1, 1}==0);
        % Get vector of labels
        
        condition_labels = {};
        condition_labels{1} = labels{9};
        condition_labels{2} = labels{17};
        condition_labels{3} = labels{31};
        condition_labels{4} = labels{37};
    case 'DBSTRD002'
        
        stim_idx_lSCC = find(labels{1, 9}==1);
        stim_idx_rSCC = find(labels{1, 17}== 1);
        stim_idx_lVCVS = find(labels{1, 33} == 1);
        stim_idx_rVCVS = find(labels{1,41} == 1);
        baseline_idx = find(labels{1,9} == 0);
        
        condition_labels = {};
        condition_labels{1} = labels{9};
        condition_labels{2} = labels{17};
        condition_labels{3} = labels{33};
        condition_labels{4} = labels{41};
        
    case 'DBSTRD003'
        
        stim_idx_lSCC = find(labels{1, 4}==1);
        stim_idx_rSCC = find(labels{1, 7} == 1);
        stim_idx_lVCVS = find(labels{1, 13} == 1);
        stim_idx_rVCVS = find(labels{1, 16} == 1);
        baseline_idx = find(labels{1, 4}== 0);
        
        condition_labels = {};
        condition_labels{1} = labels{4};
        condition_labels{2} = labels{7};
        condition_labels{3} = labels{13};
        condition_labels{4} = labels{4};
        
end


%% Get data for each condition - change this for zscored data.

data_lSCC_diff = squeeze(mean(data_zscore(condition_labels{1}==1,:,:),1)) - squeeze(mean(data_zscore(condition_labels{1}==0,:,:),1));
data_rSCC_diff = squeeze(mean(data_zscore(condition_labels{2}==1,:,:),1)) - squeeze(mean(data_zscore(condition_labels{2}==0,:,:),1));
data_lVCVS_diff = squeeze(mean(data_zscore(condition_labels{3}==1,:,:),1)) - squeeze(mean(data_zscore(condition_labels{3}==0,:,:),1));
data_rVCVS_diff = squeeze(mean(data_zscore(condition_labels{4}==1,:,:),1)) - squeeze(mean(data_zscore(condition_labels{4}==0,:,:),1));


%% Generate 3D scatter plot to get RGB values.

experiment = 'rSCC'; 

% load mmvt tables
tmp = load(sprintf('%s_MMVT_labels.mat',PatientID));

freqband = {'delta','theta','alpha','beta','lowgamma','highgamma'};
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/redblue'));


x = ROI_labels.summary_ch_info.x;
y = ROI_labels.summary_ch_info.y;
z = ROI_labels.summary_ch_info.z;

disp(experiment)
condition = {};
for i = 1:length(freqband)
    
    switch experiment  
        case 'lSCC'
            condition{i} = squeeze(data_lSCC_diff(:,i));

        case 'rSCC'
            condition{i} = squeeze(data_rSCC_diff(:,i));
        case 'lVCVS' 
            condition{i} = squeeze(data_lVCVS_diff(:,i));
        case 'rVCVS' 
            condition{i} = squeeze(data_rVCVS_diff(:,i)); 
    end 
    
    num_ch = length(ch_labels);
    labels = ch_labels;
    
    
    figure()
    h = scatter3(x,y,z,70,'filled','CData',condition{i},'MarkerEdgeColor','k');
    title(freqband{i});
    
   %mycolormap = customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
   %mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
    
   %%%%%-------- SWAP OUT LINES FROM CUSTOM COLOR MAP HERE %%%%%%-------------
   
    mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
colorbar('southoutside');
colormap(mycolormap);
axis off;
   
    
    cmap = colormap(mycolormap); %can change 'bluewhitered' to 'redblue'
    
    
    
    %%%%%%% --------------------------
    Cdata = h.CData;    

    cmin = min(Cdata(:)); %can change this to -100, -50, etc.
    cmax = max(Cdata(:)); %can change this to 100,50, etc.
    %%%%%%%%%
    m = length(cmap);
    idx = fix((Cdata-cmin)/(cmax-cmin)*m)+1;
    RGB = squeeze(ind2rgb(idx,cmap));
    tmp = 1:num_ch; tmp = tmp' ;
    RGB = [tmp,RGB];

    colorbar
    caxis([-0.8 0.8])
    
    % Append table with labels
    mmvt_tbl = array2table(RGB,'VariableNames',{'Electrodes','RGB1','RGB2','RGB3'});
    mmvt_tbl.Electrodes = deblank(string(mmvt_labels));
        
    % Save table for MMVT viz.
    filename = sprintf('/Users/anushaallawala/Data/MMVT_PSD/auto%s_%s_%s.csv',PatientID, freqband{i}, experiment);
    
    %export as csv    
   %writetable(mmvt_tbl,filename,'WriteVariableNames',0);
    
    
end

%% 

% tmp1 = [1:num_ch]; 
% tmp_thetatbl = [tmp1',condition{2}]; 
% thetatbl = array2table(tmp_thetatbl,'VariableNames',{'Electrodes','theta'});
% thetatbl.Electrodes = deblank(string(mmvt_labels));
% 
% filename_mmvt_th = sprintf('/Users/anushaallawala/Data/MMVT_PSD/%s_theta_%s.csv',PatientID,experiment);
% writetable(thetatbl,filename_mmvt_th,'WriteVariableNames',0);
% 
% 
% tmp_gammatbl = [tmp1',condition{6}]; 
% highgammatbl = array2table(tmp_gammatbl,'VariableNames',{'Electrodes','highgamma'}); 
% highgammatbl.Electrodes = deblank(string(mmvt_labels)); 
% 
% filename_mmvt_hg = sprintf('/Users/anushaallawala/Data/MMVT_PSD/%s_highgamma_%s.csv',PatientID,experiment);
% writetable(highgammatbl,filename_mmvt_hg,'WriteVariableNames',0);


