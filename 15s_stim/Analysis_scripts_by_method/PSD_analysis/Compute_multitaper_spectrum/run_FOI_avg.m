clear
addpath(genpath('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Processed Data/PSD/'))
filedir = dir('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Processed Data/PSD/DBSTRD006_rVCVS_f130_blk-*_mtspectrum_pow_indivtr.mat'); 
experiment = 'rVCVS_f130'; 
for i = 1:numel(filedir) 
    data{i} = load(filedir(i).name); 
    PSD_data{i} = data{i}.S_poststim; 
end 

%%  concatenate files. 

all_data = cat(1,PSD_data{:}); 
f_poststim = data{1}.f_poststim; 
%% 

FOI = 'delta'; 
switch FOI 
    case 'delta'
        lower_freq = 1; 
        upper_freq = 4; 
    case 'alpha'
        lower_freq = 8; 
        upper_freq = 13; 
    case 'theta'
        lower_freq = 4; 
        upper_freq = 7; 
    case 'beta'       
        lower_freq = 13; 
        upper_freq = 30; 
    case 'lowgamma' 
        lower_freq = 30; 
        upper_freq = 55; 
    case 'highgamma' 
        lower_freq = 65; 
        upper_freq = 100; 
end 

%poststim
freq_poststim_idx = find(f_poststim >lower_freq & f_poststim < upper_freq);
freq_poststim_f = f_poststim(freq_poststim_idx);
%% Get average across FOI 
%fxn handle 
avgFOI = @(pow,FOI_idx)(mean(pow(:,:,FOI_idx),3)); 
FOI_poststim_avg = avgFOI(all_data,freq_poststim_idx); 
disp('......Averaging data in FOI')

%% metadata 

metadata.channel_info = data{1}.summary_ch_info; 
metadata.firstblock = data{1}.metadata; 
metadata.secondblock = data{2}.metadata; 
metadata.thirdblock = data{3}.metadata; 
metadata.FOI_lower_freq = lower_freq; 
metadata.FOI_upper_freq = upper_freq; 
%% save 

filename = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD006/EXP/15s_stim/Processed Data/PSD_avg/%s_%s.mat',...
    experiment,FOI) ; 

save(filename, 'metadata','all_data','FOI_poststim_avg','experiment'); 








