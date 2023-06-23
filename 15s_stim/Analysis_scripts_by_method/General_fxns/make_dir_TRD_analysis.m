%% purpose: creates directories in preprocessing folders depending on workstation being used 
%%% meant to work for preprocessing_task_main_PCCT.m script 
%%% created for UH3 Depression Project by AA 05/2020 

function [outputdir] = make_dir_TRD_analysis(PatientID,exp,exp_name,analysis,sub_analysis) 
%% 

if nargin == 5
    
    outputdir = sprintf('E:/DBSTRD/%s/Experiments/%s/Processed Data/%s/%s/%s',...
        PatientID,exp,analysis,sub_analysis,exp_name);
    
elseif nargin == 4
    
    outputdir = sprintf('E:/DBSTRD/%s/Experiments/%s/Processed Data/%s/%s/',...
        PatientID,exp,analysis,exp_name);
end
    
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end

addpath(genpath(outputdir)); 
disp(outputdir) 
      
end 


