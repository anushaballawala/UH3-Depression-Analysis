%% Generate_labels_all_comps.m

% The goal of this script is to generate labels for all combinations so 
% that they can be loaded altogether when running permutations for all 
% conditions tested against each other 


%% %% Load file that contains data and descriptive labels for each trial. 
clear 
clc 
close all 

PatientID = 'DBSTRD003'; 

addpath(genpath('/Users/anushaallawala/Data/')); 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats_052022.mat',PatientID)) ; 
data = data_file.data; 
DBS_labels = string(data_file.DBS_labels); 
stim_labels = string(data_file.stim_labels); 
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

%% Labels 

comp = [];
% comp.DBS = ;
% comp.elec = '';
% comp.Stim_state ='';
comp_id = 0;
comps = {};
%DBS = {'both_hemi_SCC_VCVS', 'left_SCC_VCVS', 'right_SCC_VCVS'};
DBS = {'SCC','VCVS','left_VCVS','right_VCVS','left_SCC','right_SCC'}; 
Stim_state = 'poststim';

for dbs_idx = 1:numel(DBS)
    comp_id = comp_id + 1;
    
    comp.DBS = DBS{dbs_idx};
    comp.Stim_state = Stim_state;
    
    label_vector = zeros(num_trials, 1);
    for i = 1:num_trials
        
        is_baseline = contains(tbl.DBS_target(i),'Baseline'); 
        dbs_matches_SCC = contains(tbl.DBS_target(i), 'SCC');
        dbs_matches_VCVS = contains(tbl.DBS_target(i), 'VCVS');
        
        dbs_matches_rSCC =  contains(tbl.DBS_target(i), 'rSCC');
        dbs_matches_rVCVS = contains(tbl.DBS_target(i), 'rVCVS');
        
        dbs_matches_lSCC = contains(tbl.DBS_target(i), 'lSCC');
        dbs_matches_lVCVS =contains(tbl.DBS_target(i), 'lVCVS');  
        
        Stim_state_matches = strcmp(tbl.Stim_Labels(i), comp.Stim_state) == 1;
        
        switch dbs_idx 
            case 1
                if dbs_matches_SCC == 1 && Stim_state_matches
                    label_vector(i) = 1; 
                elseif is_baseline == 1
                    label_vector(i) = 0; 
                else
                    label_vector(i) = nan; 
                end 
            case 2 
                if dbs_matches_lSCC == 1 && Stim_state_matches
                    label_vector(i) = 1; 
                elseif dbs_matches_lVCVS ==1 && Stim_state_matches
                    label_vector(i) = 2; 
                elseif is_baseline == 1
                    label_vector(i) = 0; 
                else
                    label_vector(i) = nan; 
                end 
                
            case 3 
                if dbs_matches_rSCC == 1 && Stim_state_matches
                    label_vector(i) = 1; 
                elseif dbs_matches_rVCVS ==1 && Stim_state_matches
                    label_vector(i) = 2; 
                elseif is_baseline == 1
                    label_vector(i) = 0; 
                else
                    label_vector(i) = nan; 
                end               
        end 
    end
    
    comp.label_vector = label_vector;
    comps{comp_id} = comp;
    
end



%% Save 


save(sprintf('/Users/anushaallawala/Data/DBSTRD/labels_for_poststim_DBStargetvstarget_conds_%s_new',...
    PatientID), 'comps');  

