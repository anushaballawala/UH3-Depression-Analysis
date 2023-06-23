%% %% run_compute_PSD_FOI_ROI_bl.m - Run this to get avg power in freq bands 
% in ROIs of baseline recordings 

% Additional info: ***** filename will require fixing on Oscar *****
% Additional info: ***** outputdir will require fixing on Oscar *****
% Inputs: data from running Chronux multitaper fxn 
% Outputs: Avg spectral power in each frequency band, in "fake" prestim,
% stim, poststim states 
% Dependencies: UH3 github repo 
% Sub-functions:  outputdir fxn

% Anusha Allawala, 10/2021

%------------ START OF CODE ---------------% 




function [] = compute_PSD_FOI_ROI_bl(name_of_file, FOI, ROI_labels,PatientID,...
    experiment_name)

load(name_of_file); 

%% Define FOI freq band and input 
switch FOI 
    case 'delta'
        lower_freq = 1; 
        upper_freq = 4; 
    case 'theta'
        lower_freq = 4; 
        upper_freq = 7; 
    case 'alpha'
        lower_freq = 8; 
        upper_freq = 13; 
    case 'beta'
        lower_freq = 13; 
        upper_freq = 30; 
    case 'gamma' 
        lower_freq = 50; 
        upper_freq = 80; 
end

%% add to metadata 

metadata.FOI.upper_freq = upper_freq; 
metadata.FOI.lower_freq = lower_freq; 
metadata.FOI.freqband = FOI; 
%% Get FOI indices 

%prestim 
freq_prestim_idx = find(f_prestim >lower_freq & f_prestim < upper_freq); 
freq_prestim_f = f_prestim(freq_prestim_idx);

% stim1. 
freq_stim1_idx = find(f_stim1 >lower_freq & f_stim1 < upper_freq);
freq_stim1_f = f_stim1(freq_stim1_idx);

%poststim
freq_poststim_idx = find(f_poststim >lower_freq & f_poststim < upper_freq);
freq_poststim_f = f_poststim(freq_poststim_idx);

%% Avg PSD across FOI 

avgFOI = @(pow,FOI_idx)(mean(pow(:,:,FOI_idx),3)); 

avg_prestim_FOI = avgFOI(S_prestim,freq_prestim_idx); 
avg_stim1_FOI = avgFOI(S_stim1,freq_stim1_idx);
avg_stim2_FOI = avgFOI(S_stim2,freq_stim1_idx);
avg_stim3_FOI = avgFOI(S_stim3,freq_stim1_idx); 
avg_poststim_FOI = avgFOI(S_poststim,freq_poststim_idx); 

% save
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD_Analysis/%s',...
    PatientID,experiment_name); 
filename = sprintf('%s/avg_FOI_allch_%s',outputdir,FOI); 

if ~exist(outputdir,'dir'); mkdir(outputdir); end;
disp('saving.....') 
save(filename,'avg_prestim_FOI','avg_stim1_FOI','avg_stim2_FOI',...
    'avg_stim3_FOI','avg_poststim_FOI','freq_prestim_f','freq_stim1_f',...
    'freq_poststim_f','metadata'); 
disp('done saving')


%% Get PSD for ROI 

prestim_ROI = generate_ROI_indiv_tr_baseline(avg_prestim_FOI,PatientID); 
stim1_ROI = generate_ROI_indiv_tr_baseline(avg_stim1_FOI,PatientID); 
stim2_ROI = generate_ROI_indiv_tr_baseline(avg_stim2_FOI,PatientID); 
stim3_ROI = generate_ROI_indiv_tr_baseline(avg_stim3_FOI,PatientID); 
poststim_ROI = generate_ROI_indiv_tr_baseline(avg_poststim_FOI,PatientID); 
% save 
disp('saving...') 
 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD_Analysis/%s',...
    PatientID,experiment_name); 
filename = sprintf('%s/avg_FOI_ROI_%s',outputdir,FOI);
if ~exist(outputdir,'dir'); mkdir(outputdir); end;

save(filename, 'prestim_ROI','stim1_ROI','stim2_ROI','stim3_ROI',...
    'poststim_ROI','freq_prestim_f','freq_stim1_f',...
    'freq_poststim_f','ROI_labels','metadata') 
disp('done saving') 
end 




