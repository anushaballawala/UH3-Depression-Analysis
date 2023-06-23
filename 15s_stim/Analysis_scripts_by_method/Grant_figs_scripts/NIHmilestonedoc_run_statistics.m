clear all 
PatientID = 'DBSTRD002'; 
experiment_name = 'lSCC_f130'; 
FOI = 'theta'; 
load(sprintf('%s_%s_%s_ROIforstats.mat',PatientID,experiment_name,FOI));  

% might be diff if VCVS ** 
contact_configuration = {'elec1','elec234','elec25','elec36','elec47','elec567','elec8'}; 

%% 

switch PatientID 
    case 'DBSTRD001'
for i = 1:7
    ROI_lACC(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,1);
    ROI_lAMY(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,2);
    ROI_lDPF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,3);
    ROI_lLOF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,4);
    ROI_lMOF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,5);
    ROI_lMTG(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,6);
    ROI_lVPF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,7);
    ROI_rACC(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,8);
    ROI_rAMY(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,9);
    ROI_rDPF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,10);
    ROI_rLOF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,11);
    ROI_rMOF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,12);
    ROI_rMTG(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,13);
    ROI_rVPF(:,i) =  ROI_stim_minus_pre{1,i}{1,1}(:,14);
end
    case 'DBSTRD002' 
       for i = 1:7
           ROI_lACC(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,1);
           ROI_lAMY(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,2);
           ROI_lDPF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,3);
           ROI_lLOF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,4);
           ROI_lMOF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,5);
           ROI_lSFG(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,6);
           ROI_lSTG(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,7);
           ROI_lVPF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,8);
           ROI_rACC(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,9);
           ROI_rAMY(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,10);
           ROI_rDPF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,11);
           ROI_rLOF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,12);
           ROI_rMOF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,13);
           ROI_rPHG(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,14);
           ROI_rSFG(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,15);
           ROI_rSTG(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,16);
           ROI_rVPF(:,i) = ROI_stim_minus_pre{1,i}{1,1}(:,17);
       end 
end 
     
%%  convert all to table 
switch PatientID 
    case 'DBSTRD001' 
lACC_tbl = array2table(ROI_lACC,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
lAMY_tbl  = array2table(ROI_lAMY,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
lDPF_tbl = array2table(ROI_lDPF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
lLOF_tbl = array2table(ROI_lLOF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
lMOF_tbl = array2table(ROI_lMOF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
lMTG_tbl = array2table(ROI_lMTG,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
lVPF_tbl = array2table(ROI_lVPF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rACC_tbl = array2table(ROI_rACC,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rAMY_tbl = array2table(ROI_rAMY,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rDPF_tbl = array2table(ROI_rDPF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rLOF_tbl = array2table(ROI_rLOF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rMOF_tbl = array2table(ROI_rMOF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rMTG_tbl = array2table(ROI_rMTG,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
rVPF_tbl = array2table(ROI_rVPF,'VariableNames',{'elec1','elec234','elec25',...
    'elec36','elec47','elec567','elec8'}); 
    case 'DBSTRD002' 
        disp('not switching to table') 
end 
%% run statistics 

switch PatientID 
    case 'DBSTRD001' 
          
        [p,anovatab,stats] = anova1(ROI_lACC); 
        multcompare(stats)
        title('lACC')
        
        [p,anovatab,stats] = anova1(ROI_lAMY); 
        multcompare(stats)
        title('lamy')
        
        [p,anovatab,stats] = anova1(ROI_lDPF); 
        multcompare(stats)
        title('ldpf')
        
        [p,anovatab,stats] = anova1(ROI_lLOF); 
        multcompare(stats)
        title('lLOF')
        [p,anovatab,stats] = anova1(ROI_lMOF); 
        multcompare(stats)
        title('lMOF')
         [p,anovatab,stats] = anova1(ROI_lMTG); 
        multcompare(stats)
        title('lMTG')

        [p,anovatab,stats] = anova1(ROI_lVPF); 
        multcompare(stats)
        title('lVPF')
        
        % right hemi 
        
        [p,anovatab,stats] = anova1(ROI_rACC); 
        multcompare(stats)
        title('racc')
        [p,anovatab,stats] = anova1(ROI_rAMY); 
        multcompare(stats)
        title('rAMY')        
        [p,anovatab,stats] = anova1(ROI_rDPF); 
        multcompare(stats)
        title('rDPF')        
        [p,anovatab,stats] = anova1(ROI_rLOF); 
        multcompare(stats)
        title('rLOF')        
         [p,anovatab,stats] = anova1(ROI_rMOF); 
        multcompare(stats)
        title('rMOF')       

        [p,anovatab,stats] = anova1(ROI_rMTG); 
        multcompare(stats)
        title('rMTG')
        [p,anovatab,stats] = anova1(ROI_rVPF); 
        multcompare(stats)
        title('rVPF')       
        
        
        
        
    case 'DBSTRD002' 
         
        [p,anovatab,stats] = anova1(ROI_lACC); 
        multcompare(stats)
        title('lACC')
        
        [p,anovatab,stats] = anova1(ROI_lAMY); 
        multcompare(stats)
        title('lamy')
        
        [p,anovatab,stats] = anova1(ROI_lDPF); 
        multcompare(stats)
        title('ldpf')
        
        [p,anovatab,stats] = anova1(ROI_lLOF); 
        multcompare(stats)
        title('lLOF')
        [p,anovatab,stats] = anova1(ROI_lMOF); 
        multcompare(stats)
        title('lMOF')
         [p,anovatab,stats] = anova1(ROI_lSTG); 
        multcompare(stats)
        title('lSTG')
         [p,anovatab,stats] = anova1(ROI_lSFG); 
        multcompare(stats)
        title('lSFG')
        [p,anovatab,stats] = anova1(ROI_lVPF); 
        multcompare(stats)
        title('lVPF')
        
        % right hemi 
        
        [p,anovatab,stats] = anova1(ROI_rACC); 
        multcompare(stats)
        title('racc')
        [p,anovatab,stats] = anova1(ROI_rAMY); 
        multcompare(stats)
        title('rAMY')        
        [p,anovatab,stats] = anova1(ROI_rDPF); 
        multcompare(stats)
        title('rDPF')        
        [p,anovatab,stats] = anova1(ROI_rLOF); 
        multcompare(stats)
        title('rLOF')        
         [p,anovatab,stats] = anova1(ROI_rMOF); 
        multcompare(stats)
        title('rMOF')       
        [p,anovatab,stats] = anova1(ROI_rSFG); 
        multcompare(stats)
        title('rSFG')
        [p,anovatab,stats] = anova1(ROI_rSTG); 
        multcompare(stats)
        title('rSTG')
        [p,anovatab,stats] = anova1(ROI_rVPF); 
        multcompare(stats)
        title('rVPF')

end 
        

%% Repeat for post minus pre 

switch PatientID
    case 'DBSTRD001' 
for i = 1:7
    ROI_lACC_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,1);
    ROI_lAMY_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,2);
    ROI_lDPF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,3);
    ROI_lLOF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,4);
    ROI_lMOF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,5);
    ROI_lMTG_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,6);
    ROI_lVPF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,7);
    ROI_rACC_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,8);
    ROI_rAMY_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,9);
    ROI_rDPF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,10);
    ROI_rLOF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,11);
    ROI_rMOF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,12);
    ROI_rMTG_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,13);
    ROI_rVPF_post(:,i) =  ROI_post_minus_pre{1,i}{1,1}(:,14);
end
    case 'DBSTRD002'
        for i = 1:7
           ROI_lACC_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,1);
           ROI_lAMY_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,2);
           ROI_lDPF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,3);
           ROI_lLOF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,4);
           ROI_lMOF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,5);
           ROI_lSFG_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,6);
           ROI_lSTG_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,7);
           ROI_lVPF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,8);
           ROI_rACC_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,9);
           ROI_rAMY_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,10);
           ROI_rDPF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,11);
           ROI_rLOF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,12);
           ROI_rMOF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,13);
           ROI_rPHG_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,14);
           ROI_rSFG_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,15);
           ROI_rSTG_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,16);
           ROI_rVPF_post(:,i) = ROI_post_minus_pre{1,i}{1,1}(:,17);
       end
                 
end 

%% run statistics 
switch PatientID 
    case 'DBSTRD001' 
          
        [p,anovatab,stats] = anova1(ROI_lACC_post); 
        multcompare(stats)
        title('lACC')
        
        [p,anovatab,stats] = anova1(ROI_lAMY_post); 
        multcompare(stats)
        title('lamy')
        
        [p,anovatab,stats] = anova1(ROI_lDPF_post); 
        multcompare(stats)
        title('ldpf')
        
        [p,anovatab,stats] = anova1(ROI_lLOF_post); 
        multcompare(stats)
        title('lLOF')
        [p,anovatab,stats] = anova1(ROI_lMOF_post); 
        multcompare(stats)
        title('lMOF')
         [p,anovatab,stats] = anova1(ROI_lMTG_post); 
        multcompare(stats)
        title('lMTG')

        [p,anovatab,stats] = anova1(ROI_lVPF_post); 
        multcompare(stats)
        title('lVPF')
        
        % right hemi 
        
        [p,anovatab,stats] = anova1(ROI_rACC_post); 
        multcompare(stats)
        title('racc')
        [p,anovatab,stats] = anova1(ROI_rAMY_post); 
        multcompare(stats)
        title('rAMY')        
        [p,anovatab,stats] = anova1(ROI_rDPF_post); 
        multcompare(stats)
        title('rDPF')        
        [p,anovatab,stats] = anova1(ROI_rLOF_post); 
        multcompare(stats)
        title('rLOF')        
         [p,anovatab,stats] = anova1(ROI_rMOF_post); 
        multcompare(stats)
        title('rMOF')       

        [p,anovatab,stats] = anova1(ROI_rMTG_post); 
        multcompare(stats)
        title('rMTG')
        [p,anovatab,stats] = anova1(ROI_rVPF_post); 
        multcompare(stats)
        title('rVPF')       
        
        
        
        
    case 'DBSTRD002' 
         
        [p,anovatab,stats] = anova1(ROI_lACC_post); 
        multcompare(stats)
        title('lACC')
        
        [p,anovatab,stats] = anova1(ROI_lAMY_post); 
        multcompare(stats)
        title('lamy')
        
        [p,anovatab,stats] = anova1(ROI_lDPF_post); 
        multcompare(stats)
        title('ldpf')
        
        [p,anovatab,stats] = anova1(ROI_lLOF_post); 
        multcompare(stats)
        title('lLOF')
        [p,anovatab,stats] = anova1(ROI_lMOF_post); 
        multcompare(stats)
        title('lMOF')
         [p,anovatab,stats] = anova1(ROI_lSTG_post); 
        multcompare(stats)
        title('lSTG')
         [p,anovatab,stats] = anova1(ROI_lSFG_post); 
        multcompare(stats)
        title('lSFG')
        [p,anovatab,stats] = anova1(ROI_lVPF_post); 
        multcompare(stats)
        title('lVPF')
        
        % right hemi 
        
        [p,anovatab,stats] = anova1(ROI_rACC_post); 
        multcompare(stats)
        title('racc')
        [p,anovatab,stats] = anova1(ROI_rAMY_post); 
        multcompare(stats)
        title('rAMY')        
        [p,anovatab,stats] = anova1(ROI_rDPF_post); 
        multcompare(stats)
        title('rDPF')        
        [p,anovatab,stats] = anova1(ROI_rLOF_post); 
        multcompare(stats)
        title('rLOF')        
         [p,anovatab,stats] = anova1(ROI_rMOF_post); 
        multcompare(stats)
        title('rMOF')       
        [p,anovatab,stats] = anova1(ROI_rSFG_post); 
        multcompare(stats)
        title('rSFG')
        [p,anovatab,stats] = anova1(ROI_rSTG_post); 
        multcompare(stats)
        title('rSTG')
        [p,anovatab,stats] = anova1(ROI_rVPF_post); 
        multcompare(stats)
        title('rVPF')

end 



%%