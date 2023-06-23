function [matrix_ROI,ROI_labels] = generate_ROI_from_ch(data,PatientID,ch_dim)

    %bothhemi_matrix_ROI, bothhemi_ROI_labels] = generate_ROI_from_ch(data,PatientID,ch_dim,metadata)

% set up ROIs 
if ch_dim == 1
    avg_ROI = @(data,data_idx)(nanmean(data(data_idx,:),ch_dim)); 
elseif ch_dim ==2 
    avg_ROI = @(data,data_idx)(nanmean(data(:,data_idx,:),ch_dim)); 
elseif ch_dim ==3
   avg_ROI = @(data,data_idx)(nanmean(data(:,:,data_idx),ch_dim)); 
elseif ch_dim ==4 % avg coherence data 
    avg_ROI = @(data,data_idx)(nanmean(nanmean(data(:,data_idx,data_idx,:),2),3));
end 


    switch PatientID 
        case 'DBSTRD001' 
                load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/Electrodes Info/ROI_labels_mat/%s_ROI_labels.mat',PatientID)); 
          %%%%%% avg hemispheres 
                  ACC = avg_ROI(data,ACC_idx_ch);
                  AMY = avg_ROI(data,AMY_idx_ch); 
                  DPF = avg_ROI(data,DPF_idx_ch); 
                  LOF = avg_ROI(data,LOF_idx_ch); 
                  MOF = avg_ROI(data,MOF_idx_ch); 
                  MTG = avg_ROI(data,MTG_idx_ch); 
                  SFG = avg_ROI(data,SFG_idx_ch); 
                  VPF = avg_ROI(data,VPF_idx_ch); 
                  
            %%%%%%Average 
            
            matrix_ROI = cat(ch_dim,ACC,AMY,DPF,LOF,MOF,MTG,SFG,VPF);
            ROI_labels = {'ACC','AMY','DPF','LOF','MOF','MTG','SFG','VPF'};             

          %%%%% left hemi 
                  lACC = avg_ROI(data,lACC_idx_ch); 
                  lAMY = avg_ROI(data,lAMY_idx_ch); 
                  lDPF = avg_ROI(data,lDPF_idx_ch); 
                  lLOF = avg_ROI(data,lLOF_idx_ch); 
                  lMOF = avg_ROI(data,lMOF_idx_ch); 
                  lMTG = avg_ROI(data,lMTG_idx_ch); 
                  lSFG = avg_ROI(data,lSFG_idx_ch); 
                  lVPF = avg_ROI(data,lVPF_idx_ch); 
                  
          %%%%% right hemi 
                  rACC = avg_ROI(data,rACC_idx_ch); 
                  
                  rAMY = avg_ROI(data,rAMY_idx_ch); 
                  rDPF = avg_ROI(data,rDPF_idx_ch); 
                  rLOF = avg_ROI(data,rLOF_idx_ch); 
                  rMOF = avg_ROI(data,rMOF_idx_ch); 
                  rMTG = avg_ROI(data,rMTG_idx_ch); 
                  rSFG = avg_ROI(data,rSFG_idx_ch); 
                  rVPF = avg_ROI(data,rVPF_idx_ch); 
                  
            %%%%%%% All 
           
            bothhemi_matrix_ROI = cat(ch_dim,lACC,lAMY,lDPF,lLOF,lMOF,lMTG,lSFG,lVPF,...
                rACC,rAMY,rDPF,rLOF,rMOF,rMTG,rSFG,rVPF);
            
            bothhemi_ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG','lSFG','lVPF',...
                'rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rSFG','rVPF'};           
               
        case 'DBSTRD002'
                load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/Electrodes Info/ROI_labels_mat/%s_ROI_labels.mat',PatientID)); 

                %%%%%% avg hemispheres 
%                 tmp_LOF_idx = [1:5,7:11];
%                 tmp_AMY_idx = [1:4,6];
%                 tmp_rAMY_idx = [1,3];
%                 tmp_rLOF_idx = [2:6];
%                 
%                 LOF_idx_ch = LOF_idx_ch(tmp_LOF_idx);
%                 AMY_idx_ch = AMY_idx_ch(tmp_AMY_idx); 
%                 rAMY_idx_ch = rAMY_idx_ch(tmp_rAMY_idx);
%                 rLOF_idx_ch = rLOF_idx_ch(tmp_rLOF_idx);
                
                MOF = avg_ROI(data,MOF_idx_ch); 
                LOF = avg_ROI(data,LOF_idx_ch); 
                VPF = avg_ROI(data,VPF_idx_ch); 
                ACC = avg_ROI(data,ACC_idx_ch); 
                AMY = avg_ROI(data,AMY_idx_ch); 
                DPF = avg_ROI(data,DPF_idx_ch); 
                STG = avg_ROI(data,STG_idx_ch); 
                SFG = avg_ROI(data,SFG_idx_ch); 
                PHG = avg_ROI(data,PHG_idx_ch); 
                EC = avg_ROI(data,EC_idx_ch); 
                
                %58, 67(LOF,rLOF), 102 (deleted from rAMY)
                %Average 
                
                matrix_ROI = cat(ch_dim,MOF,LOF,VPF,ACC,AMY,DPF,STG,SFG,PHG,EC); 
                ROI_labels = {'MOF','LOF','VPF','ACC','AMY','DPF','STG','SFG','PHG','EC'}; 
                        %%%%%% left hemi 

                lMOF = avg_ROI(data,lMOF_idx_ch); 
                lLOF = avg_ROI(data,lLOF_idx_ch);  
                lVPF = avg_ROI(data,lVPF_idx_ch);  
                lACC = avg_ROI(data,lACC_idx_ch); 
                lAMY = avg_ROI(data,lAMY_idx_ch); 
                lDPF = avg_ROI(data,lDPF_idx_ch);   
                lSTG = avg_ROI(data,lSTG_idx_ch); 
                lSFG = avg_ROI(data,lSFG_idx_ch); 
                %lPHG_idx_ch  
                lEC = avg_ROI(data,lEC_idx_ch);   

                        %%%%%%% right hemi  
                rMOF = avg_ROI(data,rMOF_idx_ch); 
                rLOF = avg_ROI(data,rLOF_idx_ch); 
                rVPF = avg_ROI(data,rVPF_idx_ch); 
                rACC = avg_ROI(data,rACC_idx_ch); 
                rAMY = avg_ROI(data,rAMY_idx_ch); 
                rDPF = avg_ROI(data,rDPF_idx_ch);
                rSTG = avg_ROI(data,rSTG_idx_ch); 
                rSFG = avg_ROI(data,rSFG_idx_ch);
                rPHG = avg_ROI(data,rPHG_idx_ch);  
                %rEC_idx_ch
            
                bothhemi_matrix_ROI = cat(ch_dim,lMOF,lLOF,lVPF,lACC,lAMY,lDPF,lSTG,lSFG,lEC,...
                    rMOF,rLOF,rVPF,rACC,rAMY,rDPF,rSTG,rSFG,rPHG); 
                bothhemi_ROI_labels = {'lMOF','lLOF','lVPF','lACC','lAMY','lDPF','lSTG','lSFG','lEC',...
                    'rMOF','rLOF','rVPF','rACC','rAMY','rDPF','rSTG','rSFG','rPHG'};
                
        case 'DBSTRD003' 
            % both hemis 
            
            % Load ROI labels 
                load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/Electrodes Info/ROI_labels_mat/%s_ROI_labels.mat',PatientID)); 
%             ACC_idx_ch = [9:10,47];
%             DPF_idx_ch = [6:8,41:46];
%             LOF_idx_ch = [28:31];
%             MFG_idx_ch = [23:27,59:62];
%             MOF_idx_ch = [15:16,51:53];
%             SFG_idx_ch = [11:14,48:50];
%             VPF_idx_ch = [1:4,32:37];
            
            
            
            
%             ACC_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'aCC') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey'));
%             DPF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'DPF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey')); 
%             LOF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'LOF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey'));
%             MFG_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'MFG') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey')); 
%             MOF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'MOF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey'));
%             SFG_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'SFG') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey'));
%             VPF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'VPF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey'));
%         
%             % right hemi 
%             rACC_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'aCC') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             rDPF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'DPF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             rLOF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'LOF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             rMFG_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'MFG') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             rMOF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'MOF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             rSFG_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'SFG') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             rVPF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'VPF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Right'));
%             % left hemi 
%             lACC_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'aCC') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
%             lDPF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'DPF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
%             lLOF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'LOF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
%             lMFG_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'MFG') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
%             lMOF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'MOF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
%             lSFG_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'SFG') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
%             lVPF_idx_ch = find(contains(metadata.preprocessing.ChannelTbl.ROI, 'VPF') & contains(metadata.preprocessing.ChannelTbl.Matter_vis, 'Grey') & contains(metadata.preprocessing.ChannelTbl.Hemisphere, 'Left'));
            %index both hemis 
            ACC = avg_ROI(data, ACC_idx_ch);
            DPF = avg_ROI(data, DPF_idx_ch); 
            LOF = avg_ROI(data, LOF_idx_ch);
           % MFG = avg_ROI(data, MFG_idx_ch);
            MOF = avg_ROI(data, MOF_idx_ch);
            SFG = avg_ROI(data, SFG_idx_ch);
            VPF = avg_ROI(data, VPF_idx_ch);
%             % index right hemi 
%             rACC = avg_ROI(data, rACC_idx_ch);
%             rDPF = avg_ROI(data, rDPF_idx_ch); 
%             %rLOF = avg_ROI(data, rLOF_idx_ch);
%             rMFG = avg_ROI(data, rMFG_idx_ch);
%             rMOF = avg_ROI(data, rMOF_idx_ch);
%             rSFG = avg_ROI(data, rSFG_idx_ch);
%             rVPF = avg_ROI(data, rVPF_idx_ch);
            % index left hemi 
%             lACC = avg_ROI(data, lACC_idx_ch);
%             lDPF = avg_ROI(data, lDPF_idx_ch); 
%             lLOF = avg_ROI(data, lLOF_idx_ch);
%             lMFG = avg_ROI(data, lMFG_idx_ch);
%             lMOF = avg_ROI(data, lMOF_idx_ch);
%             lSFG = avg_ROI(data, lSFG_idx_ch);
%             lVPF = avg_ROI(data, lVPF_idx_ch);
            
            %matrix_ROI = cat(ch_dim,ACC,DPF,LOF,MOF,MFG,SFG,VPF); 
            %ROI_labels = {'ACC','DPF','LOF','MOF','MFG','SFG','VPF'};
            
            matrix_ROI = cat(ch_dim,ACC,DPF,LOF,MOF,SFG,VPF); 
            ROI_labels = {'ACC','DPF','LOF','MOF','SFG','VPF'};
            
%             
%             bothhemi_matrix_ROI = cat(ch_dim,lACC,lDPF,lLOF,lMOF,lMFG,lSFG,lVPF,...
%                rACC,rDPF,rMOF,rMFG,rSFG,rVPF);
%             bothhemi_ROI_labels = {'lACC','lDPF','lLOF','lMOF','lMFG','lSFG','lVPF',...
%                'rACC','rDPF','rMOF','rMFG','rSFG','rVPF'}; 
            
            
            
                
 
    end 





end 

