%% Load runsheet for each block 

clear 
disp('cleared')
clc
addpath(genpath('/Users/anushaallawala/Downloads/'))

% files = dir('/Users/anushaallawala/Downloads/Reg_amp_15s/f130/*.xlsx'); 
% 
% file = {}; 
% 
% for i = 1:numel(files) 
%     file{i} = files(i).name; 
% end 

%% 
filename = sprintf('/Users/anushaallawala/Downloads/Reg_amp_15s/f130/DBSTRD008_Stim15s_run-lVCVS_blk-f130_subblk-1.xlsx'); 
filename_new = sprintf('/Users/anushaallawala/Downloads/Low2_amp_15s/f130/DBSTRD008_Stim15s_run-lVCVS_blk-f130_subblk-1.xlsx'); 

old_T = readtable(filename); 

%% Reformat and replace amplitude values. 

% First sweep through values for single contact pw 180 

for i = 1:height(old_T) 
    if strcmp(old_T.Amplitude{i},'4800') == 1 
        disp('replacing') 
        old_T.Amplitude{i} = '2400'; 
        disp('replaced') 
    end 
end 

% then stack configs pw 180 
for i = 1:height(old_T) 
    if strcmp(old_T.Amplitude{i},'2400 2400') == 1 
        disp('replacing') 
        old_T.Amplitude{i} = '1200 1200'; 
        disp('replaced') 
    end 
end 

% then ring configs pw 180 
for i = 1:height(old_T) 
    if strcmp(old_T.Amplitude{i},'1600 1600 1600') == 1 
        disp('replacing') 
       old_T.Amplitude{i} = '800 800 800'; 
        disp('replaced') 
    end 
end 

%% Delete condition column. 


%new_T = removevars(old_T, 'Condition'); 
new_T = old_T; 

%% save file name in file_directory with low amplitude 

writetable(new_T, filename_new); 

