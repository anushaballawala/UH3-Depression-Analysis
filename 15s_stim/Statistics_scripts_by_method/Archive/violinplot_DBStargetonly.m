%clear 
clc 
close all 
clear 
PatientID = 'DBSTRD002'; 
addpath(genpath('/Users/anushaallawala/Data/')); 
data_file = load(sprintf('%s_data_for_PSD_stats',PatientID)) ; 
data = data_file.all_data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_state_labels); 
num_trials = size(data,1); 
tbl = data_file.tbl; 
%Load channels labels that will be used later. 
load(sprintf('%s_goodch.mat',PatientID)) 
%*** add new path for 002 good channels 

%% **** Fix this in the preprocessing code, temporary code for now ***** 
switch PatientID 
    case 'DBSTRD001'
        disp('change nothing')
    case 'DBSTRD002'
        tmp_idx = [1:57,59:66,68:101,103:129];
        data = data(:,tmp_idx,:); 
end 

%% Convert descriptive labels to binary labels in a vector for each combo being tested. 

comp_labels = load(sprintf('labels_for_poststim_conds_%s_new',PatientID)); 
comp_labels = comp_labels.comps;

labels = {}; 
for i = 1:numel(comp_labels)
    labels{i} = comp_labels{i}.label_vector;
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
    ~, ~] = generate_ROI_from_ch(data,PatientID,2);
end 
% load matrix ROI data 
%% Get data mean for plotting later 

%create 3D matrix with data for each stim exp corresponding to labels. 
matrix_ROI_alldata_label1 = cell(1,numel(comp_labels)); 
matrix_ROI_alldata_label0 = cell(1,numel(comp_labels));

for i = 1:numel(comp_labels)
         
    matrix_ROI_alldata_label1{i} = matrix_ROI(labels{i}==1,:,:); 
    matrix_ROI_alldata_label0{i} = matrix_ROI(labels{i}==0,:,:); 
         
end 

%matrix_ROI_alldata_label1 contains both SCC and VCVS data at this point. 
%% What is this for ???** 
% switch PatientID
%     case 'DBSTRD001' 
        SCC_idx = 24; 
        VCVS_idx = 25;
%     case 'DBSTRD002'
%         SCC_idx = 24; 
%         VCVS_idx = 25; 
% end 
%% Get labels for stim experiments for SCC and VCVS  

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

%% Get VCVS and SCC data out from ROI_matrix_alldata_label1 var

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
% 
SCC_data = alldata_SCC{1, 1}  ; 
VCVS_data = alldata_VCVS{1, 1}  ;
baseline_data = alldata_baseline; 

%% 
data_labels_SCC = repmat('SCC ',size(SCC_data(:,1,1))); 
data_labels_VCVS = repmat('VCVS', size(VCVS_data(:,1,1))); 
data_labels_baseline = repmat('bl  ', size(baseline_data(:,1,1))); 

% concatenate 
complete_data = vertcat(SCC_data, VCVS_data, baseline_data); 
complete_labels = vertcat(data_labels_SCC,data_labels_VCVS,data_labels_baseline); 
assert(length(complete_data)==length(complete_labels)); 

%% Make violin plot for DBS target vs target 

addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/Violinplot-Matlab-master'))

% for i = 1:numel(ROI_labels)
% figure() 
% violinplot(complete_data(:,i,1), cellstr(complete_labels)); 
% title(ROI_labels{i})
% end 

%% Plot difference data 

% subtract data 

SCC_diff = SCC_data - mean(baseline_data,1); 
VCVS_diff = VCVS_data - mean(baseline_data,1); 

% concatenate 

diff_data = vertcat(SCC_diff, VCVS_diff); 
diff_labels_DBStarget = vertcat(data_labels_SCC,data_labels_VCVS); 

diff_data_contactconfig = vertcat(SCC_diff,SCC_diff,VCVS_diff,VCVS_diff); 
%% Labels for contact config 
SCC_all = repmat('SCC          ',size(alldata_SCC{1, 1}(:,1,1))); 
SCC_1 = repmat('SCC elec1    ',size(alldata_SCC{1, 2}(:,1,1))); 
SCC_25 = repmat('SCC elec25   ',size(alldata_SCC{1, 3}(:,1,1))); 
SCC_36 = repmat('SCC elec36   ',size(alldata_SCC{1, 4}(:,1,1))); 
SCC_47= repmat('SCC elec47   ',size(alldata_SCC{1, 5}(:,1,1))); 
SCC_8= repmat('SCC elec8    ',size(alldata_SCC{1, 6}(:,1,1))); 
SCC_234= repmat('SCC elec234  ',size(alldata_SCC{1, 7}(:,1,1))); 
SCC_567= repmat('SCC elec567  ',size(alldata_SCC{1, 8}(:,1,1))); 

VCVS_all = repmat('VCVS         ',size(alldata_VCVS{1, 1}(:,1,1))); 
VCVS_1= repmat('VCVS elec1   ',size(alldata_VCVS{1, 2}(:,1,1))); 
VCVS_25= repmat('VCVS elec25  ',size(alldata_VCVS{1, 3}(:,1,1))); 
VCVS_36= repmat('VCVS elec36  ',size(alldata_VCVS{1, 4}(:,1,1))); 
VCVS_47= repmat('VCVS elec47  ',size(alldata_VCVS{1, 5}(:,1,1))); 
VCVS_8= repmat('VCVS elec8   ',size(alldata_VCVS{1, 6}(:,1,1))); 


if numel(alldata_VCVS)>18
    VCVS_234= repmat('VCVS elec234 ',size(alldata_VCVS{1, 7}(:,1,1))); 
    VCVS_567= repmat('VCVS elec567 ',size(alldata_VCVS{1, 8}(:,1,1))); 
diff_labels_contactconfig = vertcat(SCC_all,SCC_1,SCC_25,SCC_36,SCC_47,SCC_8,SCC_234,SCC_567,...
    VCVS_all,VCVS_1,VCVS_25,VCVS_36,VCVS_47,VCVS_8,VCVS_234,VCVS_567);    

else 
diff_labels_contactconfig = vertcat(SCC_all,SCC_1,SCC_25,SCC_36,SCC_47,SCC_8,SCC_234,SCC_567,...
    VCVS_all,VCVS_1,VCVS_25,VCVS_36,VCVS_47,VCVS_8); 
end 


%% Make plot for all contact configurations 

freqband = {'delta','theta','alpha','beta','lowgamma','highgamma'}; 

for j = 1:numel(freqband)
    for i = 1:numel(ROI_labels)
    figure() 
    violinplot(diff_data_contactconfig(:,i,j), cellstr(diff_labels_contactconfig)); 
    title(ROI_labels{i})
    ylabel('Power difference from baseline')
    end
end 
saveas(gcf, sprintf('%s_%s_%s.png',PatientID,freqband{j},ROI_labels{i}))



