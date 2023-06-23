
function [outputdir] = make_directory(PatientID, type, dir_name, experiment_name, contact_config)

switch type 
    case 'figure'
    outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/%s/%s/%s',...
    PatientID,dir_name,experiment_name,contact_config);
if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end 
   addpath(genpath(outputdir));

    case 'data'
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s/%s/%s',...
    PatientID,dir_name,experiment_name,contact_config);

if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end 
   addpath(genpath(outputdir));
end  
   
   
   
    
end


