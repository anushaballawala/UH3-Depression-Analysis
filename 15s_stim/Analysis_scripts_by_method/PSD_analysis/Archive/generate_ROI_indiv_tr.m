function [matrix_ROI] = generate_ROI_indiv_tr(stim_sig)
% set up ROIs 

lMOF = [1:3,26:31]; 
lLOF = [4:8]; 
lVPF = [9:12];
lDPF = [18:24,32:37,44:45];
lACC = [39:43];
lMTG = [52:58];
lAMY = [47:51]; 

rLOF = [67:70]; 
rVPF = [71:75]; 
rDPF = [86:90,99:105,113:117]; 
rACC = [107:112];
rMTG = [126:131];
rAMY = [120:125];
rMOF = [62:66,77:80,92:98]; 

lACC_stim = mean(stim_sig(:,lACC,:),2); 
lAMY_stim = mean(stim_sig(:,lAMY,:),2); 
lDPF_stim = mean(stim_sig(:,lDPF,:),2); 
lLOF_stim = mean(stim_sig(:,lLOF,:),2); 
lMOF_stim = mean(stim_sig(:,lMOF,:),2); 
lMTG_stim = mean(stim_sig(:,lMTG,:),2); 
lVPF_stim = mean(stim_sig(:,lVPF,:),2);  

rACC_stim = mean(stim_sig(:,rACC,:),2);
rAMY_stim = mean(stim_sig(:,rAMY,:),2); 
rDPF_stim = mean(stim_sig(:,rDPF,:),2); 
rLOF_stim = mean(stim_sig(:,rLOF,:),2); 
rMOF_stim = mean(stim_sig(:,rMOF,:),2); 
rMTG_stim = mean(stim_sig(:,rMTG,:),2); 
rVPF_stim = mean(stim_sig(:,rVPF,:),2);  

matrix_ROI = horzcat(lACC_stim,lAMY_stim,lDPF_stim,lLOF_stim,lMOF_stim,lMTG_stim,...
    lVPF_stim,...
    rACC_stim,rAMY_stim,rDPF_stim,rLOF_stim,rMOF_stim,rMTG_stim,rVPF_stim); 
% matrix_ROI = [lACC_stim,lAMY_stim,lDPF_stim,lLOF_stim,lMOF_stim,lMTG_stim,...
%     lVPF_stim,...
%     rACC_stim,rAMY_stim,rDPF_stim,rLOF_stim,rMOF_stim,rMTG_stim,rVPF_stim];


end 