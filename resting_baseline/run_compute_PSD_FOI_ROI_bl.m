%% run_compute_PSD_FOI_ROI_bl.m - Run this to get avg power in freq bands 
% in ROIs of baseline recordings 

% Additional info: ***** filename will require fixing on Oscar *****
% Additional info: ***** outputdir will require fixing on Oscar *****
% Inputs: N/A
% Outputs: Spectral power in ROIs 
% Dependencies: Chronux toolbox, UH3 github repo 
% Sub-functions: compute_PSD_FOI_ROI_bl, outputdir fxn

% Anusha Allawala, 10/2021

%------------ START OF CODE ---------------% 
%% Load PSD data 

PatientID = 'DBSTRD002'; 
experiment = 'BaselineFix'; 
experiment_name = 'BaselineFix_run-01_blk-01'; 

name_of_file = '11-Oct-2021_BaselineFix_run-01_blk-01_mtspectrum_pow_indivtr.mat'; 
ROI_labels = load(sprintf('ROI_labels_%s.mat',PatientID)); 

%% define freq band of interest 

FOI = 'delta'; 
compute_PSD_FOI_ROI(name_of_file, FOI, ROI_labels, PatientID,...
    experiment_name); 

%% plot ROI 1/f 




