%% Generate heatmaps of coherence data after it has been subtracted by bl.
close all 
i = 1;
addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/customcolormap'));
addpath(genpath('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/'));

DBStarget = 'lSCC';

freqname = {'delta','theta','alpha','beta','low gamma','high gamma'};

%ch_labels = string(ROI_labels.summary_ch_info.Label);
num_ch = numel(ch_labels); 
num_freq = numel(freqname);
%for i = 1:num_freq
freqband = freqname{i};

switch DBStarget
    case 'baseline' 
        figure('Position',[100,100,1200,1000]);
        imagesc(squeeze(data_baseline(:,:,:,i)))
    case 'lSCC'
        figure('Position',[100,100,1200,1000]);
        imagesc(squeeze(data_lSCC_diff(:,:,i)))
        %imagesc(squeeze(data_lSCC(:,:,:,i)));
        
    case 'rSCC'
        figure('Position',[100,100,1200,1000]);
        imagesc(squeeze(data_rSCC_diff(:,:,i)))
        
    case 'lVCVS'
        figure('Position',[100,100,1200,1000]);
        imagesc(squeeze(data_lVCVS_diff(:,:,i)))
        
    case 'rVCVS'
        figure('Position',[100,100,1200,1000]);
        imagesc(squeeze(data_rVCVS_diff(:,:,i)))
        
end

xticks([1:num_ch]);
xticklabels(ch_labels);
yticks([1:num_ch]);
yticklabels(ch_labels);
%define colormap
% mycolormap = customcolormap(linspace(0,1,11), {'#f8f697','#fac023','#f37816','#d64844','#b63457',...
%     '#902567','#731b6e', '#5d126c', '#380862', '#220c4a', '#060418' });
colorbar('eastoutside');
        colormap redbluecmap ; 

%colormap(mycolormap);
%axis off;
%caxis([-0.1 0.1])
caxis([0.1 0.3])
title(sprintf('%s %s %s', PatientID, DBStarget, freqband))
directory = '/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Coherence';
filename = sprintf('%s/%s_%s_%s_bl_subtr_%s.fig',directory,PatientID,DBStarget,freqband,date);
filename_svg = sprintf('%s/%s_%s_%s_bl_subtr_%s.svg',directory,PatientID,DBStarget,freqband,date);
% saveas(gcf,filename);
% saveas(gcf,filename_svg);

%end

%% 

