%% Generate_heat_map.m 
%%%%%% The goal of this script is to generate heatmaps (current direction
%%%%%% by ROI) for each of the DBS targets being stimulated %%%%%%
%%% INPUT: trial x frequencies x channel x time 
%%% OUTPUT: Visuals: heatmap with averaged power in freqbands for ROIs
%%%         Data: cell array containing averaged data *, with metadata 

% set up directory 
files = dir('E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Epoched Data/*/*15s_stim_all_currdir_singletrial_*.mat');
[datafile, data] = load_files_from_dir(files); 
%% 
num_files = length(datafile); 
% enter FOI 
freqs = 13:18; 
freqband = 'beta'; 

%% set up variable names 
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

%% Perform averaging and baseline correction for each cell entry for each DBS target 

DBStarget = 'lSCC'; 
data_target = lSCC; 
[lSCC_stim_e1,lSCC_poststim1_e1,lSCC_poststim2_e1,...
    lSCC_stim_e234,lSCC_poststim1_e234,lSCC_poststim2_e234,...
    lSCC_stim_e25,lSCC_poststim1_e25,lSCC_poststim2_e25,...
    lSCC_stim_e36,lSCC_poststim1_e36,lSCC_poststim2_e36,...
    lSCC_stim_e47,lSCC_poststim1_e47,lSCC_poststim2_e47,...
    lSCC_stim_e567,lSCC_poststim1_e567,lSCC_poststim2_e567,...
    lSCC_stim_e8,lSCC_poststim1_e8,lSCC_poststim2_e8,...
    lSCC_stim_e, lSCC_poststim1_e, lSCC_poststim2_e] = compute_for_each_contact_config(DBStarget,data_target,freqs) ; 


%% Repeat for rSCC

DBStarget = 'rSCC'; 
data_target = rSCC; 
[ rSCC_stim_e1, rSCC_poststim1_e1, rSCC_poststim2_e1,...
     rSCC_stim_e234, rSCC_poststim1_e234, rSCC_poststim2_e234,...
     rSCC_stim_e25, rSCC_poststim1_e25, rSCC_poststim2_e25,...
     rSCC_stim_e36, rSCC_poststim1_e36, rSCC_poststim2_e36,...
     rSCC_stim_e47, rSCC_poststim1_e47, rSCC_poststim2_e47,...
     rSCC_stim_e567, rSCC_poststim1_e567, rSCC_poststim2_e567,...
     rSCC_stim_e8, rSCC_poststim1_e8, rSCC_poststim2_e8,...
     rSCC_stim_e, rSCC_poststim1_e, rSCC_poststim2_e] = compute_for_each_contact_config(DBStarget,data_target,freqs); 

%% Repeat for lVCVS 

DBStarget = 'lVCVS'; 
data_target = lVCVS; 
[ lVCVS_stim_e1, lVCVS_poststim1_e1, lVCVS_poststim2_e1,...
     lVCVS_stim_e25, lVCVS_poststim1_e25, lVCVS_poststim2_e25,...
     lVCVS_stim_e36, lVCVS_poststim1_e36, lVCVS_poststim2_e36,...
     lVCVS_stim_e47, lVCVS_poststim1_e47, lVCVS_poststim2_e47,...
     lVCVS_stim_e8, lVCVS_poststim1_e8, lVCVS_poststim2_e8,...
     lVCVS_stim_e, lVCVS_poststim1_e, lVCVS_poststim2_e] = compute_for_each_contact_config_VCVS(DBStarget,data_target,freqs); 

%% Repeat for rVCVS 

DBStarget = 'rVCVS'; 
data_target = rVCVS; 
[rVCVS_stim_e1,rVCVS_poststim1_e1,rVCVS_poststim2_e1,...
    rVCVS_stim_e25,rVCVS_poststim1_e25,rVCVS_poststim2_e25,...
    rVCVS_stim_e36,rVCVS_poststim1_e36,rVCVS_poststim2_e36,...
    rVCVS_stim_e47,rVCVS_poststim1_e47,rVCVS_poststim2_e47,...
    rVCVS_stim_e8,rVCVS_poststim1_e8,rVCVS_poststim2_e8,...
    rVCVS_stim_e, rVCVS_poststim1_e, rVCVS_poststim2_e] = compute_for_each_contact_config_VCVS(DBStarget,data_target,freqs); 

%% load ROI 

load('E:/DBSTRD/DBSTRD001/Experiments/ROI_labels_DBSTRD001.mat'); 
%% extract ROIs

% lSCC
[ROI_lSCC_stim] = generate_ROI(lSCC_stim_e); % gives ROI x current dir matrix 
[ROI_lSCC_poststim1] = generate_ROI(lSCC_poststim1_e); 
[ROI_lSCC_poststim2] = generate_ROI(lSCC_poststim2_e); 
% rSCC
[ROI_rSCC_stim] = generate_ROI(rSCC_stim_e); 
[ROI_rSCC_poststim1] = generate_ROI(rSCC_poststim1_e); 
[ROI_rSCC_poststim2] = generate_ROI(rSCC_poststim2_e); 
%lVCVS
[ROI_lVCVS_stim] = generate_ROI(lVCVS_stim_e); 
[ROI_lVCVS_poststim1] = generate_ROI(lVCVS_poststim1_e); 
[ROI_lVCVS_poststim2] = generate_ROI(lVCVS_poststim2_e); 
%rVCVS
[ROI_rVCVS_stim] = generate_ROI(rVCVS_stim_e); 
[ROI_rVCVS_poststim1] = generate_ROI(rVCVS_poststim1_e); 
[ROI_rVCVS_poststim2] = generate_ROI(rVCVS_poststim2_e); 

%% generate heatmaps 

ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG','lVPF',...
    'rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rVPF'}; 
SCC_labels = {'e1','e25', 'e36', 'e47', 'e234', 'e567', 'e8'};
VCVS_labels = {'e1', 'e25', 'e36', 'e47', 'e8'}; 
%% 
%lSCC during stim 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_lSCC_stim'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = SCC_labels; 
h.ColorLimits = [-50 50]; 
title('lSCC during stim') 
filename = 'lSCC_stim.png'; 
saveas(gcf, filename)

%lSCC - poststim1 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_lSCC_poststim1'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = SCC_labels; 
h.ColorLimits = [-50 50]; 
title('lSCC for 5s post-stim') 
filename = 'lSCC_poststim1.png'; 
saveas(gcf, filename)

%lSCC - poststim2
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_lSCC_poststim2'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = SCC_labels; 
h.ColorLimits = [-50 50]; 
title('lSCC for 5s post-stim2') 
filename = 'lSCC_poststim2.png'; 
saveas(gcf, filename)

%% rSCC 

%lSCC during stim 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_rSCC_stim'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = SCC_labels; 
h.ColorLimits = [-50 50]; 
title('rSCC during stim') 
filename = 'rSCC_stim.png'; 
saveas(gcf, filename)

%rSCC - poststim1 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_rSCC_poststim1'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = SCC_labels; 
h.ColorLimits = [-50 50]; 
title('rSCC for 5s post-stim') 
filename = 'rSCC_poststim1.png'; 
saveas(gcf, filename)

%rSCC - poststim2
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_rSCC_poststim2'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = SCC_labels; 
h.ColorLimits = [-50 50]; 
title('rSCC for 5s post-stim2') 
filename = 'rSCC_poststim2.png'; 
saveas(gcf, filename)

%% lVCVS

%lVCVS during stim 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_lVCVS_stim'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = VCVS_labels; 
h.ColorLimits = [-50 50]; 
title('lVCVS during stim') 
filename = 'lVCVS_stim.png'; 
saveas(gcf, filename)

%lVCVS - poststim1 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_lVCVS_poststim1'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = VCVS_labels; 
h.ColorLimits = [-50 50]; 
title('lVCVS for 5s post-stim') 
filename = 'lVCVS_poststim1.png'; 
saveas(gcf, filename)

%lVCVS - poststim2
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_lVCVS_poststim2'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = VCVS_labels; 
h.ColorLimits = [-50 50]; 
title('lVCVS for 5s post-stim2') 
filename = 'lVCVS_poststim2.png'; 
saveas(gcf, filename)

%% rVCVS

%rVCVS during stim 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_rVCVS_stim'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = VCVS_labels; 
h.ColorLimits = [-50 50]; 
title('rVCVS during stim') 
filename = 'rVCVS_stim.png'; 
saveas(gcf, filename)

%rVCVS - poststim1 
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_rVCVS_poststim1'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = VCVS_labels; 
h.ColorLimits = [-50 50]; 
title('rVCVS for 5s post-stim') 
filename = 'rVCVS_poststim1.png'; 
saveas(gcf, filename)

%rVCVS - poststim2
figure('Position',[500 300 800 600]) 
h = heatmap(ROI_rVCVS_poststim2'); 
colormap redbluecmap 
co.YDisplayLabels = '% change BL';
h.XDisplayLabels = ROI_labels; 
h.YDisplayLabels = VCVS_labels; 
h.ColorLimits = [-50 50]; 
title('rVCVS for 5s post-stim2') 
filename = 'rVCVS_poststim2.png'; 
saveas(gcf, filename)

%% 

filename_data = sprintf('ROI_power_%s',freqband); 
save(filename_data, 'ROI_lSCC_stim', 'ROI_lSCC_poststim1', 'ROI_lSCC_poststim2',...
    'ROI_rSCC_stim', 'ROI_rSCC_poststim1', 'ROI_rSCC_poststim2',...
    'ROI_lVCVS_stim', 'ROI_lVCVS_poststim1', 'ROI_lVCVS_poststim2',...
    'ROI_rVCVS_stim', 'ROI_rVCVS_poststim1', 'ROI_rVCVS_poststim2',...
    'SCC_labels', 'VCVS_labels', 'ROI_labels', 'freqs', 'freqband'); 

% append table with all values in it for stim, poststim1 and poststim2 
% cellarray for elec1 
   

% cellarray for elec2_5 
% lSCC_elec1 = 
% rSCC_elec1 = 
% lVCVS_elec1 = 
% rVCVS elec1 = 
 
% make new variables for each current dir containing all SCC & VCVS data 





