% load averaged data 

% set up directory 

% files = dir('E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Epoched Data/*/*trialaveraged*.mat');
% [datafile, data] = load_files_from_dir(files); 

%% 
clear all 
target = 'rSCC';
%% 

if strcmp(target,'rVCVS') == 1 
    elec1 = (sprintf('15s_stim_elec9_%s_f130.mat',target));
    elec8 = load(sprintf('15s_stim_elec16_%s_f130.mat',target));
    elec2_5 = load(sprintf('15s_stim_elec10_13_%s_f130.mat',target));
    elec3_6 = load(sprintf('15s_stim_11_14_%s_f130.mat',target));
    elec4_7 = load(sprintf('15s_stim_12_15_%s_f130.mat',target));
    %assign values 
    elec1 = elec1.indiv1; 
    elec8 = elec8.indiv2; 
    elec2_5 = elec2_5.indiv3; 
    elec3_6 = elec3_6.indiv4; 
    elec4_7 = elec4_7.indiv5;
end 

if strcmp(target,'rSCC') == 1 
    elec1 = load(sprintf('15s_stim_cond1_%s_f130.mat',target));
    elec8 = load(sprintf('15s_stim_cond7_%s_f130.mat',target));
    elec2_5 = load(sprintf('15s_stim_cond3_%s_f130.mat',target));
    elec3_6 = load(sprintf('15s_stim_cond4_%s_f130.mat',target));
    elec4_7 = load(sprintf('15s_stim_cond5_%s_f130.mat',target));
    elec2_3_4 = load(sprintf('15s_stim_cond2_%s_f130.mat',target));
    elec5_6_7 = load(sprintf('15s_stim_cond6_%s_f130.mat',target));
    %assign values 
    elec1 = elec1.indiv1; 
    elec8 = elec8.indiv7; 
    elec2_5 = elec2_5.indiv3; 
    elec3_6 = elec3_6.indiv4; 
    elec4_7 = elec4_7.indiv5; 
    elec2_3_4 = elec2_3_4.indiv2;
    elec5_6_7 = elec5_6_7.indiv6; 
end 


freqs = 8:13; 
freqband = 'alpha'; 

[stim_elec1,poststim1_elec1,poststim2_elec1] = perform_single_tr_norm(elec1,freqs); 
[stim_elec8,poststim1_elec8,poststim2_elec8] = perform_single_tr_norm(elec8,freqs); 
[stim_elec2_5,poststim1_elec2_5,poststim2_elec2_5] = perform_single_tr_norm(elec2_5,freqs); 
[stim_elec3_6,poststim1_elec3_6,poststim2_elec3_6] = perform_single_tr_norm(elec3_6,freqs); 
[stim_elec4_7,poststim1_elec4_7,poststim2_elec4_7] = perform_single_tr_norm(elec4_7,freqs); 
[stim_elec2_3_4,poststim1_elec2_3_4,poststim2_elec2_3_4] = perform_single_tr_norm(elec2_3_4,freqs); 
[stim_elec5_6_7,poststim1_elec5_6_7,poststim2_elec5_6_7] = perform_single_tr_norm(elec5_6_7,freqs); 



%% get ROI values 

[ROI_elec1] = generate_ROI(poststim2_elec1);
[ROI_elec8] = generate_ROI(poststim2_elec8);
[ROI_elec2_5] = generate_ROI(poststim2_elec2_5);
[ROI_elec3_6] = generate_ROI(poststim2_elec3_6);
[ROI_elec4_7] = generate_ROI(poststim2_elec4_7);

if strcmp(target,'rSCC') == 1 
    [ROI_elec2_3_4] = generate_ROI(poststim2_elec2_3_4);
    [ROI_elec5_6_7] = generate_ROI(poststim2_elec5_6_7);
end 


%% combine all ROIs into one 

%all_target_ROI = [ROI_elec1;ROI_elec8;ROI_elec2_5;ROI_elec3_6;ROI_elec4_7]; % add 234 and 567 
all_target_ROI = [ROI_elec1;ROI_elec8;ROI_elec2_5;ROI_elec3_6;ROI_elec4_7;ROI_elec2_3_4;ROI_elec5_6_7]; % add 234 and 567 

%% 
labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG','lVPF',...
    'rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rVPF'}; 


h = heatmap(all_target_ROI); 
h.XDisplayLabels = labels; 
h.YDisplayLabels = {'elec1','elec8','elec2 5','elec3 6','elec4 7','elec2 3 4','elec5 6 7'};  
colormap(redblue)
h.ColorLimits = [-25 100]; 
    

    
    
    
    
%%  
    
    
    
    


% fs = file.metadata.preprocessing.New_SamplingRate  ;
% 
% 
% %% Metadata 
% % Get other channel info from metadata
% metadata = file.metadata; 
% channel_label = metadata.preprocessing.GoodChannelLabels; 
% ch_idx = metadata.preprocessing.GoodChannelsIdx;
% freqs = 1:4; 
% freqband = 'delta'; 
% %% average across channels for each ROI? 
% 
% [stim_sig,poststim1_sig,poststim2_sig] = perform_single_tr_norm(data,freqs); 
% 
% %% 
% 
% filename = sprintf('ROI_%s_%s_currdir7.mat',freqband,target); 
% save(filename,'matrix_ROI','metadata'); 
% average percent change across channels for each ROI 

% plot heatmap or square map of percent value for each ROI 
% plot heatmap of each ROI across time (similar to before) before, during
% and after stim 


%% 







