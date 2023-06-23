
clear 
clc 
SCC = load('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Processed Data/PSD_avg/rSCC_f130_highgamma.mat');
VCVS = load('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Processed Data/PSD_avg/rVCVS_f130_highgamma.mat');
roi_info = SCC.metadata.channel_info; 
roi_names = unique(roi_info.ROI_vis); 
num_ch = height(roi_info); 
FOI = 'highgamma'
%% get ROI names 

clear roi_idx ch_idx SCC_avg VCVS_avg
   j_roi = 6; 
   currentroi = string(roi_names(j_roi)); 
   disp(currentroi)
   
for ichan = 1:num_ch 
   if (strcmp(roi_info.ROI_vis{ichan}, string(roi_names{j_roi})) == 1 && strcmp(roi_info.Matter_vis{ichan}, 'Grey') == 1);
    
   roi_idx(ichan) = ichan; 
   end 
end 

ch_idx = find(roi_idx>0) ; 

SCC_avg = mean(SCC.FOI_poststim_avg(:,ch_idx),2);  
VCVS_avg = mean(VCVS.FOI_poststim_avg(:,ch_idx),2);  


% plot 

x = [ones(1,70) 2*ones(1,70)]; 
y = [SCC_avg', VCVS_avg'];
y_z = zscore(y); 

swarmchart(x(1:70),y_z(1:70),30,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
hold on 
swarmchart(x(71:140),y_z(71:140),30,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
scatter(1,mean(y_z(1:70)),60,'filled','k')
scatter(2,mean(y_z(71:140)),60,'filled','k')

ylabel('zscore power')
%xlim([0.8 2.1])
title(sprintf('%s %s',currentroi,FOI)) 
legend('SCC','VCVS')

% saveas(gcf,'/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Figures/PSD_avg/%s',...
%     figname)
%% 




