%%% The purpose of this script is to generate line plots for each ROI
%%% showing pre, during, and post stim effects to assess variability in the
%%% data between stimulation ON and OFF states 
%%% INPUT: Epoched single trial data, tr x ch x freqs x time 
%%% OUTPUT: vector for each ROI *** 


% set up directory 
files = dir('E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Epoched Data/*/*15s_stim_all_currdir_singletrial_*.mat');
[datafile, data] = load_files_from_dir(files); 
%% 
num_files = length(datafile); 
% enter FOI 
freqs = 4:7; 
freqband = 'theta'; 

%% 
%find name of the DBS target for each file and assign var w/corresponding name 
for i = 1:num_files 
    if find(contains(data{i, 1}.metadata.files.RawDataFile,'lSCC')) == 1
        disp('this file is lSCC')
        lSCC = data{i,1}; 
    elseif find(contains(data{i, 1}.metadata.files.RawDataFile,'rSCC')) == 1
        disp('this file is rSCC') 
        rSCC = data{i,1}; 
    elseif find(contains(data{i, 1}.metadata.files.RawDataFile,'lVCVS')) == 1
        disp('this file is lVCVS') 
        lVCVS = data{i,1}; 
    elseif find(contains(data{i, 1}.metadata.files.RawDataFile,'rVCVS')) == 1
        disp('this file is rVCVS') 
        rVCVS = data{i,1};           
    end 
end 
%% Perform 
DBStarget = 'rSCC'; 
data_target = rSCC; 
[rSCC_stim_e1,rSCC_poststim1_e1,rSCC_poststim2_e1,...
    rSCC_stim_e234,rSCC_poststim1_e234,rSCC_poststim2_e234,...
    rSCC_stim_e25,rSCC_poststim1_e25,rSCC_poststim2_e25,...
    rSCC_stim_e36,rSCC_poststim1_e36,rSCC_poststim2_e36,...
    rSCC_stim_e47,rSCC_poststim1_e47,rSCC_poststim2_e47,...
    rSCC_stim_e567,rSCC_poststim1_e567,rSCC_poststim2_e567,...
    rSCC_stim_e8,rSCC_poststim1_e8,rSCC_poststim2_e8] = compute_single_tr_for_each_contact_config(DBStarget,data_target,freqs) ; 

num_ch = size(rSCC_poststim1_e1,2); 


%% concatenate - maybe don't need? 

stim_e = horzcat(rSCC_stim_e1, rSCC_stim_e25, rSCC_stim_e36, rSCC_stim_e47, rSCC_stim_e234, rSCC_stim_e567, rSCC_stim_e8);
poststim1_e = horzcat(rSCC_poststim1_e1, rSCC_poststim1_e25,...
    rSCC_poststim1_e36, rSCC_poststim1_e47, rSCC_poststim1_e234, rSCC_poststim1_e567, rSCC_poststim1_e8);
poststim2_e = horzcat(rSCC_poststim2_e1, rSCC_poststim2_e25,...
    rSCC_poststim2_e36, rSCC_poststim2_e47, rSCC_poststim2_e234, rSCC_poststim2_e567, rSCC_poststim2_e8);




