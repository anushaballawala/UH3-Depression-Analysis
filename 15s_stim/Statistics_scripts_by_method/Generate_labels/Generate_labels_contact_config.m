%% 

% The goal of this script is to generate labels for all combinations so 
% that they can be loaded altogether when running permutations for all 
% conditions tested against each other 


%% %% Load file that contains data and descriptive labels for each trial. 
clear 
clc 
close all 

PatientID = 'DBSTRD001'; 

addpath(genpath('/Users/anushaallawala/Data/')); 
data_file = load(sprintf('%s_autocorrbl_data_for_PSD_stats.mat',PatientID)) ; 
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

% generate labels from 1-7 or 1-5 for contact configurations so that anova
% can be run within a lead to compare conditions 

comp = [];  
% comp.DBS = ;
% comp.elec = '';
% comp.Stim_state ='';
comp_id = 0;
comps = {}; 
% DBS = {'SCC', 'VCVS'}; 


DBS = {'both_hemi_SCC_VCVS', 'left_SCC_VCVS', 'right_SCC_VCVS'};

Stim_state = {'poststim'}; 

for dbs_idx = 1:numel(DBS)

        comp_id = comp_id + 1;
        
        comp.DBS = DBS{dbs_idx};
        comp.Stim_state = Stim_state;
        
        label_vector = zeros(num_trials, 1);
        for i = 1:num_trials
            
            dbs_matches = contains(tbl.DBS_target(i), comp.DBS);
            Stim_state_matches = strcmp(tbl.Stim_Labels(i), comp.Stim_state) == 1;
            is_baseline = contains(tbl.DBS_target(i), 'Baseline') == 1 ; 

            if dbs_matches && strcmp(tbl.Elec_label(i),'elec1') == 1 && Stim_state_matches
                label_vector(i) = 1;
            elseif dbs_matches && strcmp(tbl.Elec_label(i),'elec25') == 1 && Stim_state_matches
                label_vector(i) = 2;
            elseif dbs_matches && strcmp(tbl.Elec_label(i),'elec36') == 1 && Stim_state_matches
                label_vector(i) = 3;
            elseif dbs_matches && strcmp(tbl.Elec_label(i),'elec47') == 1 && Stim_state_matches
                label_vector(i) = 4;
            elseif dbs_matches && strcmp(tbl.Elec_label(i),'elec8') == 1 && Stim_state_matches
                label_vector(i) = 5;
            elseif dbs_matches && strcmp(tbl.Elec_label(i),'elec234') == 1 && Stim_state_matches
                label_vector(i) = 6;
            elseif dbs_matches && strcmp(tbl.Elec_label(i),'elec567') == 1 && Stim_state_matches
                label_vector(i) = 7;
            elseif is_baseline == 1
                label_vector(i) = 0; 
            else
                label_vector(i) = NaN;
            end
            
        end
               
            comp.label_vector = label_vector;
            comps{comp_id} = comp;
        end
        

%% Save 


save(sprintf('/Users/anushaallawala/Data/DBSTRD/Stats_data/labels_for_poststim_elecconfig_conds_%s_new',...
    PatientID), 'comps'); 






