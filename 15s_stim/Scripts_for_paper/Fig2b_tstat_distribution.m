%% permandhypothesistest.m 

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
%clear 

%%%%% DOUBLE CHECK # OF TRIALS AND LABELS MATCH UP FOR 002 AND THEN RERUN
%%%%% THIS SCRIPT 

clc 
close all 
clear
PatientID = 'DBSTRD001'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
%data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 
%Load channels labels that will be used later. 
load(sprintf('%s_goodch.mat',PatientID)) 

%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_stimvsoff_conds_%s_new',PatientID));
%comp_labels = load(sprintf('labels_for_poststim_elecconfig_conds_%s_new',PatientID));
comp_labels = comp_labels.comps;

%% Get labels for poststim for std. 

stim_idx_SCC = find(comp_labels{1,1}.label_vector==1);
stim_idx_VCVS = find(comp_labels{1,2}.label_vector==1); 
baseline_idx = find(comp_labels{1,1}.label_vector==0);  

SCC = data(stim_idx_SCC, :,:);  % SCC elec all, 
VCVS = data(stim_idx_VCVS, :,:); % VCVS elec all

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


%% relabel data after stim trials thrown out from zscoring 

comp_labels_old = comp_labels; 

for i = 1:numel(comp_labels)
    
    comp_labels{1,i}.label_vector = comp_labels{1,i}.label_vector(logical(nostim_idx'));  
    
end 

labels = {}; 
for i = 1:numel(comp_labels)
    labels{i} = comp_labels{i}.label_vector;
end 
%% Calculate z-score 

data_zscore = zeros(size(data_nostim)); 
for i = 1:size(data_nostim,1)
    for j = 1:size(data_nostim,2)
        for k = 1:size(data_nostim,3)
            data_zscore(i,j,k) = (data_nostim(i,j,k) - mean(data_nostim(:,j,k)))./std(data_nostim(:,j,k)); 
        end 
    end 
end 

%% Define which type of data, Channel or ROI 

data_type_list = {'ROI', 'Ch'}; 
[idx_data_list,tf_data] = listdlg('ListString', data_type_list, 'ListSize',[150,250]); 
data_list_config = data_type_list{idx_data_list}; 
disp(data_list_config)
data_dim = data_list_config; 
 
if strcmp(data_dim,'ROI')==1 
    % get data that will be input into permutations fxn. 
    [matrix_ROI,ROI_labels,...
    ~, ~] = generate_ROI_from_ch(data_zscore,PatientID,2);
end 


%% Get data mean for plotting later 

%create 3D matrix with data for each stim exp corresponding to labels. 
matrix_ROI_alldata_label1 = cell(1,numel(comp_labels)); 
matrix_ROI_alldata_label0 = cell(1,numel(comp_labels));

for i = 1:numel(comp_labels)
         
    matrix_ROI_alldata_label1{i} = matrix_ROI(labels{i}==1,:,:); 
    matrix_ROI_alldata_label0{i} = matrix_ROI(labels{i}==0,:,:); 
         
end 

%matrix_ROI_alldata_label1 contains both SCC and VCVS data at this point.

%% Get labels for stim experiments for SCC and VCVS  

SCC_idx = 24;
VCVS_idx = 25; 


comp_labels_VCVS = {}; 
comp_labels_SCC = {}; 

for i = 1:SCC_idx 
    comp_labels_SCC{i} = comp_labels{i}; 
end 

for j = VCVS_idx:numel(comp_labels)
    comp_labels_VCVS{j} = comp_labels{j}; 
end 
% delete empties because VCVS has fewer experiments than SCC for 001. 
comp_labels_VCVS = comp_labels_VCVS(~cellfun('isempty',comp_labels_VCVS)); 




%% Get VCVS and SCC data out from zscored data 

%Get VCVS and SCC data out from ROI_matrix_alldata_label1 var

alldata_SCC= {}; 
for i = 1:SCC_idx
alldata_SCC{i} = matrix_ROI_alldata_label1{i}; 
end 

alldata_VCVS = {};
for j = VCVS_idx:numel(comp_labels)
    alldata_VCVS{j} = matrix_ROI_alldata_label1{j};  
end 
% delete empties because VCVS has fewer experiments than SCC for 001. 
alldata_VCVS = alldata_VCVS(~cellfun('isempty',alldata_VCVS)); 

alldata_baseline = matrix_ROI_alldata_label0{1}; 

%% Compute mean power values for SCC and VCVS 

% matrices for each condition are current configuration X frequency band. 
for i = 1:numel(alldata_SCC)
mean_data_SCC{i} = squeeze(mean(alldata_SCC{1,i},1)); 
end 

for i = 1:numel(alldata_VCVS)
mean_data_VCVS{i} = squeeze(mean(alldata_VCVS{1,i},1)); 
end 

mean_data_baseline = squeeze(mean(alldata_baseline,1)); 

%% Run permutations and shuffle labels 

num_perm = 1000;                                                                                                                                                                                 
% t_stat structure is num_perm+1 X data_dim X num_freq. 

%shuffled_labels = cell(numel(comp_labels),1); 
shuffled_labels = cell(numel(labels),1);

%permutations done separately for each condition being compared, just put
%in one cell. 
for i = 1:numel(labels)
    switch data_dim 
        case 'Ch'
            [t_stat{i},shuffled_labels{i}] = permutations(data_zscore,labels{i},num_perm,@Tfunc); 
        case 'ROI' 
            [t_stat{i}, shuffled_labels{i}] = permutations(matrix_ROI,labels{i},num_perm,@Tfunc); 
    end 

end 

%% Save output from permutations 

metadata.stats.labels = labels; 
metadata.stats.num_perm = num_perm; 
metadata.stats.tbl = tbl;
metadata.stats.num_trials = num_trials; 
metadata.stats.DBS_labels = DBS_labels; 

metadata.PatientID = PatientID; 
% metadata.good_ch_idx = good_ch_idx; 
% metadata.good_ch_labels = good_ch_labels; 
%metadata.good_ch_note = note; 
metadata.ROI.ROI_labels = ROI_labels;

%add hard drive to path. 
addpath(genpath('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data')); 
filepath = '/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data'; 
filename = sprintf('%s/%s_t_stat_vals_ROI_offon_%s_zscore.mat',filepath,PatientID,date); 
save(filename,'data','data_zscore','data_nostim','t_stat','shuffled_labels', 'ROI_labels1','labels','shuffled_labels','comp_labels'); 

disp('t stat saved')
%% Convert t_stat into matrix of conditions 

% concatenate and make into matrix that is perm X ROI X freq X SCC conds X
% VCVS conds. 

t_stat_SCC_cell = {}; 
for i = 1:SCC_idx
t_stat_SCC_cell{i} = t_stat{i}; 
end 

t_stat_VCVS_cell = {};
for j = VCVS_idx:numel(comp_labels)
    t_stat_VCVS_cell{j} = t_stat{j}; 
end 
t_stat_VCVS_cell = t_stat_VCVS_cell(~cellfun('isempty',t_stat_VCVS_cell)); 

%% convert to matrix and concatenate

switch PatientID
    case 'DBSTRD001'
        t_stat_SCC = cat(4,t_stat_SCC_cell{:}); 
        t_stat_VCVS = cat(4,t_stat_VCVS_cell{:}); 
        tmp = zeros(size(t_stat_SCC)); 
        tmp(:,:,:,1:18) = t_stat_VCVS; 
        t_stat_VCVS = tmp; 
        %pad VCVS so that it can be concatenated 
        t_stat_mat = cat(5,t_stat_SCC,t_stat_VCVS); 
               
    case 'DBSTRD002'
        t_stat_SCC = cat(4,t_stat_SCC_cell{:}); 
        t_stat_VCVS = cat(4,t_stat_VCVS_cell{:});
        t_stat_mat = cat(5,t_stat_SCC,t_stat_VCVS); 
end 

% compute p values X ROI X freq X SCC conds X VCVS conds  


%% Hypothesis testing - get p-value 
% p values X ROI X freq X SCC conds X VCVS conds
% calculate p_value by looking at comparison of t-stat values and dividing
% by total number of permuations 

% Taking absolute value to check against both extreme ends of sample
% population. 
p_value = (squeeze(mean(abs(t_stat_mat)>=abs(t_stat_mat(1,:,:,:,:))))); 

%% save data for jupyter 

alldata_SCC_1 = alldata_SCC{1}; 
alldata_VCVS_1 = alldata_VCVS{1}; 

diff_SCC = alldata_SCC_1 - mean(alldata_baseline,1); 
diff_VCVS = alldata_VCVS_1 - mean(alldata_baseline,1); 
%% 
diff_data = cat(1,diff_SCC,diff_VCVS); 

py_filename = sprintf('/Users/anushaallawala/Python_projects/%s_zscore_bl_poststim_SCC_VCVS.mat',PatientID); 
save(py_filename,'diff_data'); 


%% ------ Put together data for MMVT -------- %% 

t_lSCC_theta = t_stat{1}(1,:,2); 
t_lSCC_theta(ROI_labels.goodch_idx_remove) = []; 

t_rSCC_theta = t_stat{2}(1,:,2); 
t_rSCC_theta(ROI_labels.goodch_idx_remove) = []; 

t_lVCVS_theta = t_stat{3}(1,:,2); 
t_lVCVS_theta(ROI_labels.goodch_idx_remove) = []; 

t_rVCVS_theta = t_stat{4}(1,:,2); 
t_rVCVS_theta(ROI_labels.goodch_idx_remove) = []; 

%%%%%%%%%%
t_lSCC_hgamma = t_stat{1}(1,:,6); 
t_lSCC_hgamma(ROI_labels.goodch_idx_remove) = []; 

t_rSCC_hgamma = t_stat{2}(1,:,6); 
t_rSCC_hgamma(ROI_labels.goodch_idx_remove) = []; 

t_lVCVS_hgamma = t_stat{3}(1,:,6); 
t_lVCVS_hgamma(ROI_labels.goodch_idx_remove) = []; 

t_rVCVS_hgamma = t_stat{4}(1,:,6); 
t_rVCVS_hgamma(ROI_labels.goodch_idx_remove) = []; 




%% 



% load mmvt tables
load(sprintf('%s_MMVT_labels.mat',PatientID));

freqband = {'delta','theta','alpha','beta','lowgamma','highgamma'};
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/redblue'));


x = ROI_labels.summary_ch_info.x;
y = ROI_labels.summary_ch_info.y;
z = ROI_labels.summary_ch_info.z;
    
    num_ch = length(ch_labels);
    labels = ch_labels;
    
    
    figure()
    h = scatter3(x,y,z,70,'filled','CData',t_lSCC_theta,'MarkerEdgeColor','k');
    %title(freqband{i});
    
   %mycolormap = customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
    %mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
    %colormap(mycolormap);
    
    cmap = colormap(redblue); %can change 'bluewhitered' to 'redblue'
    
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
    caxis([-2 2])
    
    % Append table with labels
    mmvt_tbl = array2table(RGB,'VariableNames',{'Electrodes','RGB1','RGB2','RGB3'});
    mmvt_tbl.Electrodes = deblank(string(mmvt_labels));
        
    % Save table for MMVT viz.
    filename = sprintf('/Users/anushaallawala/Data/DBSTRD/MMVT_PSD/%s_%s_%s_%s.csv',date,PatientID, 'theta', 'lSCC');
    
    %export as csv    
   writetable(mmvt_tbl,filename,'WriteVariableNames',0);
    
    

