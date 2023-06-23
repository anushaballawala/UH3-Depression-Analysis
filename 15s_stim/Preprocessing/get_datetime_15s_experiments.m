
clc 
%clear 
close all 

addpath(genpath('/Volumes/TRD_RAWDATA/DBSTRD001'))
addpath(genpath('/Volumes/TRD_RAWDATA/DBSTRD002'))

% 
% lVCVS = load('/Volumes/TRD_RAWDATA/DBSTRD001/lVCVSstim15s_01_f130/sub-DBSTRD001_task-lVCVSstim15s_run-01_blk-f130_rawData.mat'); 
% 
% rVCVS = load('/Volumes/TRD_RAWDATA/DBSTRD001/rVCVSstim15s_01_f130/sub-DBSTRD001_task-rVCVSstim15s_run-01_blk-f130_rawData.mat'); 
% 
% lSCC = load('/Volumes/TRD_RAWDATA/DBSTRD001/lSCCstim15s_01_f130/sub-DBSTRD001_task-lSCCstim15s_run-01_blk-f130_rawData.mat'); 
% 
% rSCC = load('/Volumes/TRD_RAWDATA/DBSTRD001/rSCCstim15s_08_06f130/sub-DBSTRD001_task-rSCCstim15s_run-08_blk-06f130_rawData.mat'); 

lSCC = load('/Volumes/TRD_RAWDATA/DBSTRD002/Stim15s_lSCC_f130/sub-DBSTRD002_task-Stim15s_run-lSCC_blk-f130_rawData.mat'); 
lVCVS = load('/Volumes/TRD_RAWDATA/DBSTRD002/Stim15s_lVCVS_f130/sub-DBSTRD002_task-Stim15s_run-lVCVS_blk-f130_rawData.mat'); 
rSCC = load('/Volumes/TRD_RAWDATA/DBSTRD002/Stim15s_rSCC_f130/sub-DBSTRD002_task-Stim15s_run-rSCC_blk-f130_rawData.mat'); 
rVCVS = load('/Volumes/TRD_RAWDATA/DBSTRD002/Stim15s_rVCVS_05/sub-DBSTRD002_task-Stim15s_run-rVCVS_blk-05_rawData.mat'); 




%% 

lVCVS_time_001 = lVCVS.NS3.MetaTags.DateTime; 
rVCVS_time_001 = rVCVS.NS3.MetaTags.DateTime; 
lSCC_time_001 = lSCC.NS3.MetaTags.DateTime; 
rSCC_time_001 = rSCC.NS3.MetaTags.DateTime; 

lVCVS_time_002 = lVCVS.NS3.MetaTags.DateTime; 
rVCVS_time_002 = rVCVS.NS3.MetaTags.DateTime; 
lSCC_time_002 = lSCC.NS3.MetaTags.DateTime; 
rSCC_time_002 = rSCC.NS3.MetaTags.DateTime; 





