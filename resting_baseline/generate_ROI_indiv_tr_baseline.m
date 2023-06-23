function [matrix_ROI] = generate_ROI_indiv_tr(data,PatientID)
% set up ROIs 

switch PatientID
    case 'DBSTRD001'
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
        
        lACC_stim = mean(data(:,lACC,:),2);
        lAMY_stim = mean(data(:,lAMY,:),2);
        lDPF_stim = mean(data(:,lDPF,:),2);
        lLOF_stim = mean(data(:,lLOF,:),2);
        lMOF_stim = mean(data(:,lMOF,:),2);
        lMTG_stim = mean(data(:,lMTG,:),2);
        lVPF_stim = mean(data(:,lVPF,:),2);
        
        rACC_stim = mean(data(:,rACC,:),2);
        rAMY_stim = mean(data(:,rAMY,:),2);
        rDPF_stim = mean(data(:,rDPF,:),2);
        rLOF_stim = mean(data(:,rLOF,:),2);
        rMOF_stim = mean(data(:,rMOF,:),2);
        rMTG_stim = mean(data(:,rMTG,:),2);
        rVPF_stim = mean(data(:,rVPF,:),2);
        
        matrix_ROI = horzcat(lACC_stim,lAMY_stim,lDPF_stim,lLOF_stim,lMOF_stim,lMTG_stim,...
            lVPF_stim,...
            rACC_stim,rAMY_stim,rDPF_stim,rLOF_stim,rMOF_stim,rMTG_stim,rVPF_stim);
        % matrix_ROI = [lACC_stim,lAMY_stim,lDPF_stim,lLOF_stim,lMOF_stim,lMTG_stim,...
        %     lVPF_stim,...
        %     rACC_stim,rAMY_stim,rDPF_stim,rLOF_stim,rMOF_stim,rMTG_stim,rVPF_stim];
 
        
    case 'DBSTRD002' 
        
            lMOF = [2:4,8,9,14,15]; %
            lLOF = [11:13]; %
            lDPF = [26,54]; %
            lACC = [27:30]; % 
            lSFG = [33:36];
            lSTG = [43:46]; 
            lEC = [37]; 
            lAMY = [38:39]; %
            lVPF = [47:50]; 

            rMOF = [55:59,69:72]; 
            rLOF = [60:63,66:68]; 
            rDPF = [75:83,113,118]; 
            rACC = [84:88]; 
            rSFG = [89:93]; 
            rPHG = [95]; 
            rAMY = [96:98]; 
            rSTG = [101:104]; 
            rVPF = [106:112]; 
         
            % Group channels into ROI 
                   
          lACC_stim = mean(data(:,lACC,:),1);
          lAMY_stim = mean(data(:,lAMY,:),1); 
          lDPF_stim = mean(data(:,lDPF,:),1); 
          lLOF_stim = mean(data(:,lLOF,:),1); 
          lMOF_stim = mean(data(:,lMOF,:),1); 
          lSFG_stim = mean(data(:,lSFG,:),1); 
          lSTG_stim = mean(data(:,lSTG,:),1); 
          lVPF_stim = mean(data(:,lVPF,:),1); 
          
          rACC_stim = mean(data(:,rACC,:),1); 
          rAMY_stim = mean(data(:,rAMY,:),1); 
          rDPF_stim = mean(data(:,rDPF,:),1); 
          rLOF_stim = mean(data(:,rLOF,:),1); 
          rMOF_stim = mean(data(:,rMOF,:),1); 
          rPHG_stim = mean(data(:,rPHG,:),1); 
          rSFG_stim = mean(data(:,rSFG,:),1); 
          rSTG_stim = mean(data(:,rSTG,:),1); 
          rVPF_stim = mean(data(:,rVPF,:),1); 
          
          %combine into one matrix 
          matrix_ROI = horzcat(lACC_stim,lAMY_stim,lDPF_stim,lLOF_stim,...
              lMOF_stim,lSFG_stim,lSTG_stim,lVPF_stim,rACC_stim,rAMY_stim,...
              rDPF_stim,rLOF_stim,rMOF_stim,rPHG_stim,rSFG_stim,rSTG_stim,...
              rVPF_stim); 
         

        
end 