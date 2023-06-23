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
%load(sprintf('%s_goodch.mat',PatientID)) 

%% **** Fix this in the preprocessing code, temporary code for now ***** 
% switch PatientID 
%     case 'DBSTRD001'
%         disp('change nothing')
%     case 'DBSTRD002'
%         tmp_idx = [1:57,59:66,68:101,103:129];
%         data = data(:,tmp_idx,:); 
%         
% end 

%% Labels 

comp = [];  
% comp.DBS = ;
% comp.elec = '';
% comp.Stim_state ='';
comp_id = 0;
comps = {}; 
DBS = {'SCC', 'lSCC', 'rSCC', 'VCVS', 'lVCVS', 'rVCVS'}; 
elec = {'all','elec1', 'elec25', 'elec36', 'elec47', 'elec8', 'elec234', 'elec567'}; 
Stim_state = {'poststim'}; 

for dbs_idx = 1:numel(DBS)
    for elec_idx = 1:numel(elec)
        comp_id = comp_id + 1;
        
        comp.DBS = DBS{dbs_idx};
        comp.elec = elec{elec_idx};
        comp.Stim_state = Stim_state;
        
        %comp.DBS = DBS{1};
        %comp.elec = elec{1};
        %comp.Stim_state = Stim_state{1};
        
        label_vector = zeros(num_trials, 1);
        for i = 1:num_trials
            is_baseline = contains(tbl.DBS_target(i), 'Baseline') == 1 ; 
            dbs_matches = contains(tbl.DBS_target(i), comp.DBS);
            elec_matches = strcmp(comp.elec, 'all') == 1 || strcmp(tbl.Elec_label(i), comp.elec) == 1;
            Stim_state_matches = strcmp(tbl.Stim_Labels(i), comp.Stim_state) == 1;
            
            if is_baseline
                label_vector(i) = 0;
            elseif dbs_matches && elec_matches && Stim_state_matches
                label_vector(i) = 1;
            else
                label_vector(i) = nan;
            end
        end
        % TODO: maybe this skipping is a bad idea.
        % Could check if label_vector has no ones, in that case there were
        % no trials matching the comparison conditions and so probably we
        % specified some invalid combination -- let's just skip.
        if ~any(label_vector == 1)
            disp('Skipping this combination -- EITHER ERROR OR COMBINATION WAS NOT TESTED');
            disp(comp.DBS)
            disp(comp.elec)
            % Reset comp id
            comp_id = comp_id - 1;
            % skip this comp
            continue;
        else
            comp.label_vector = label_vector;
            comps{comp_id} = comp;
        end
    end
end

%% Save 


save(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_labels/labels_for_poststim_stimvsoff_conds_%s_new',...
    PatientID), 'comps');  

