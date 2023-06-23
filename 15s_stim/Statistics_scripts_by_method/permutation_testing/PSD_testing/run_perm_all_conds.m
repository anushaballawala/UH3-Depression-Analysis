%% Permutation testing for all conditions at once 

%output is t-statistic that is num_perm X num_ch X num_freq X num_SCCcond X
%num_VCVScond 

%% ***** NEEDS additional code to test for each current config 

%% %% Load file that contains data and descriptive labels for each trial. 
clear 
clc 
close all 

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


%% Labels 

cond1 = 'SCC'; 
cond2 = 'Baseline'; 


for i = 1:num_trials
    switch comparison 
        case 'SCC_bothhemi_allelec'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 
                labels_comp1(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp1(i,1) = 0;
            else
                labels_comp1(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec1'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec1') == 1 
                labels_comp2(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp2(i,1) = 0;
            else
                labels_comp2(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec25'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec25') == 1
                labels_comp3(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp3(i,1) = 0;
            else
                labels_comp3(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec36'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec36') == 1
                labels_comp4(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp4(i,1) = 0;
            else
                labels_comp4(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec47'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec47') == 1
                labels_comp5(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp5(i,1) = 0;
            else
                labels_comp5(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec8'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec8') == 1
                labels_comp6(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp6(i,1) = 0;
            else
                labels_comp6(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec234'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec234') == 1
                labels_comp7(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp7(i,1) = 0;
            else
                labels_comp7(i,1) = NaN;
            end
        case 'SCC_bothhemi_elec567'
            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec567') == 1
                labels_comp8(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp8(i,1) = 0;
            else
                labels_comp8(i,1) = NaN;
            end
            
            
        case 'SCC_left_elec1'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 
                labels_comp9(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp9(i,1) = 0;
            else
                labels_comp9(i,1) = NaN;
            end
            
            
            %** fix 
        case 'SCC_left_elec25'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec25') == 1
                labels_comp10(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp10(i,1) = 0;
            else
                labels_comp10(i,1) = NaN;
            end
        case 'SCC_left_elec36'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec36') == 1
                labels_comp11(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp11(i,1) = 0;
            else
                labels_comp11(i,1) = NaN;
            end
        case 'SCC_left_elec47'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec47') == 1
                labels_comp12(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp12(i,1) = 0;
            else
                labels_comp12(i,1) = NaN;
            end
        case 'SCC_left_elec8'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec8') == 1
                labels_comp13(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp13(i,1) = 0;
            else
                labels_comp13(i,1) = NaN;
            end
        case 'SCC_left_elec234'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec234') == 1
                labels_comp14(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp14(i,1) = 0;
            else
                labels_comp14(i,1) = NaN;
            end
        case 'SCC_left_elec567'
            if contains(tbl.DBS(i),'lSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec567') == 1
                labels_comp15(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp15(i,1) = 0;
            else
                labels_comp15(i,1) = NaN;
            end
        case 'SCC_right_elec1'
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec1') == 1
                labels_comp16(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp16(i,1) = 0;
            else
                labels_comp16(i,1) = NaN;
            end
        case 'SCC_right_elec25'
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec25') == 1
                labels_comp17(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp17(i,1) = 0;
            else
                labels_comp17(i,1) = NaN;
            end
        case 'SCC_right_elec36'
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec36') == 1
                labels_comp18(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp18(i,1) = 0;
            else
                labels_comp18(i,1) = NaN;
            end
        case 'SCC_right_elec47'
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec47') == 1
                labels_comp19(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp19(i,1) = 0;
            else
                labels_comp19(i,1) = NaN;
            end
        case 'SCC_right_elec8'
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec8') == 1
                labels_comp20(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp20(i,1) = 0;
            else
                labels_comp20(i,1) = NaN;
            end
        case 'SCC_right_elec234'
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec234') == 1
                labels_comp21(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp21(i,1) = 0;
            else
                labels_comp21(i,1) = NaN;
            end
        case 'SCC_right_elec567'    
            if contains(tbl.DBS(i),'rSCC')==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),'elec567') == 1
                labels_comp22(i,1) = 1;
            elseif contains(tbl.DBS(i),cond2)== 1
                labels_comp22(i,1) = 0;
            else
                labels_comp22(i,1) = NaN;
            end

            
            
        case 'VCVS_bothhemi_elec1'
        case 'VCVS_bothhemi_elec25'
        case 'VCVS_bothhemi_elec36'
        case 'VCVS_bothhemi_elec47'
        case 'VCVS_bothhemi_elec8'
        case 'VCVS_bothhemi_elec234'
        case 'VCVS_bothhemi_elec567'
            
        case 'VCVS_left_elec1'
        case 'VCVS_left_elec25'
        case 'VCVS_left_elec36'
        case 'VCVS_left_elec47'
        case 'VCVS_left_elec8'
        case 'VCVS_left_elec234'
        case 'VCVS_left_elec567'
            
        case 'VCVS_right_elec1'
        case 'VCVS_right_elec25'
        case 'VCVS_right_elec36'
        case 'VCVS_right_elec47'
        case 'VCVS_right_elec8'
        case 'VCVS_right_elec234'
        case 'VCVS_right_elec567'
            
end







%labels for comparison 1, labels for comparison 2, etc. 