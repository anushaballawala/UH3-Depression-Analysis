%function epochtablestim_15s_008() 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % FUNCTION: epoch neural data file into chunks during 15s stim experiment; 
 % 1. description of stimulation parameters with timestamps from comments is
 % parsed and grouped. 
 % 2. detect start of stim from analog chanel containing waveform from .ns3
 % file 
 % 3. assign timestamps from #2 to parsed comments containing stim
 % description 
 
 % INPUT: .ns5 file containing analog channel with stim waveform, .nev file
 % containing comments with timestamps 
 
 % OUTPUT: table containing description of stimtype for each trial, and
 % corresponding timestamp of when each trial started 
 
 % Dependencies: signal processing toolbox, NPMK repo 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ~
cd data/TRD_Project/DBSTRD/DBSTRD010/
clear 
clc 
close all 
%load .mat file containing comments from Blackrock NSP (originally NEV) 
disp('Select NEV files for 15s STIM')
[filename,filepath]=uigetfile('*.nev','Select NEV File');
addpath(genpath('/users/aallawa1/Documents/MATLAB/Toolboxes/NPMK-master/'))
NEV = openNEV(fullfile(filepath,filename),'nosave');

%filename = 'sub-TRDDBS001_task-shortstim_run_01_blk-rVCVSf130elec12'; %example filename 
%% 
%get timestamps from comment file 
comments = NEV.Data.Comments; 
timestamps = comments.TimeStampSec';

%get text from comment file
comments_text = string(NEV.Data.Comments.Text);

%% parse out stim param information from stimulation comments and remove bad trials 

 header_text = 'StimOn';
 
 %find comments with stimon 
 
 txt_idx = strcmp(comments_text,header_text); 
 comments_text = deblank(comments_text); 
 new_comments = comments_text(txt_idx); 
 
N = length(new_comments); %number of total trials in this file 
 %take out garbage comments (if theres any that say there is packet loss,
 %or error .
 frequency_field = strings(N,1); %preallocating freq string 
trial_text = new_comments; 
  
%% repeat 
%parse out stim param information from stimulation comments and remove bad trials 

 header_text = 'StimOn';
 
 %find comments with stimon 
 
 txt_idx = strcmp(comments_text,header_text); 
 comments_text = deblank(comments_text); 
 new_comments = comments_text(txt_idx); 
 
N = length(new_comments); %number of total trials in this file 
 %take out garbage comments (if theres any that say there is packet loss,
 %or error .
 frequency_field = strings(N,1); %preallocating freq string 
trial_text = new_comments; 
  %%
  
  clean_comments_text = new_comments; 
 %set up parsing out freq info 
 freq_pattern = 'frequency=\d*;';
 %assign freq field to full comments, will parse in next section 
 frequency_field = clean_comments_text; % Cheat to make a string array of length N by copying.
  
 %set up parsing out amplitude info  
 amp_pattern = 'amplitude=\d*(\s*\d*\s*\d*)?;';
 amp_field = clean_comments_text; % Cheat to make a string array of length N by copying.
 
 %set up parsing out pulse width info 

 pw_pattern = 'SetPulseWidth=\d*;' ; 
 pw_field = clean_comments_text;  % Cheat to make a string array of length N by copying.
 
 %set up parsing out electrode info 
 e_pattern = 'electrode=\d*(\s*\d*\s*\d*)?;';
 e_field = clean_comments_text; % Cheat to make a string array of length N by copying.
  
 
 %% Get PatientID info 
 
 %get patientID info 
PatientID = getPatientID(filename); 
metadata.general.PatientID = PatientID; 

fprintf('Patient ID: %s\n', PatientID);
%% create strings with information for each parameter (freq,amp,PW,elec #) 

 for i = 1:N
     frequency_field(i) = '130'
%      %now for amplitude
     amp_field(i) = '3';
     
     %now for pulse width 
     pw_field(i) ='180'; 
     
     %dont forget the electrode number grandma  
     e_field(i) = 'placeholder'; 
     
 end 

 %% set up remaining info for epoch table 

braintarget = '(l|r)(VCVS|SCC)';
[start_idx_br, end_idx_br] = regexp(filename,braintarget);
DBStarget = extractBetween(filename,start_idx_br,end_idx_br);


%create Mx1 variable with condition names for each trial 
condition = cell(N,1);
for i = 1:N
tmp = sprintf('%s_elec%s_amp%s_freq%s_pw%s' ,DBStarget{1}, e_field(i), amp_field(i), frequency_field(i), pw_field(i)); 
condition{i} = tmp; 
end

%get additional information about run/block/trials 
subjectname = getPatientID(filename);  

[subjectname,PatientType,block_file,block_num,...
    freq,hemi, DBStarget] = extract_stim_info(filename) ; 
DBStarget = sprintf('%s%s',hemi,DBStarget); 
%generate epoch table 
table_as_matrix = zeros(N, 10);
table_as_matrix(:, 1) = 1:N;
%table_as_matrix(:, 2) = block_name; can't add string yet 
table_as_matrix(:, 3) = timestamps(txt_idx);
%table_as_matrix(:, 4) = condition{:}; %cant add string yet 
table_as_matrix(:, 5) = frequency_field;
table_as_matrix(:, 6) = pw_field;
table_as_matrix(:, 7) = amp_field;
table_as_matrix(:, 8) = e_field;
%table_as_matrix(:, 9) = DBStarget{1}; %cant add string yet 

tbl = array2table(table_as_matrix, 'VariableNames', ...
{'TrialNum', 'Block', 'Time', 'Condition', 'Frequency', 'PW', 'Amp',...
'Contacts', 'Stimtarget','TimeAnalog'});

% Reformat the block numbers as strings with 3 digits
block_numbers = repmat(block_num(2),1,N);
tbl.Block = block_numbers'; 
tbl.Condition = string(condition);
tbl.Frequency = frequency_field;
tbl.PW = pw_field;
tbl.Amp = amp_field;
tbl.Contacts = e_field;
tbl.Stimtarget = repelem(string(DBStarget),N)';
%% run timestamps fxn 
stim_info = sprintf('%s_%s_%s_f%s', PatientID, DBStarget, block_num, freq); 

[timestamps,cerestim,stimsync] = extract_timestamps_analog(subjectname,stim_info); 

% check edges detected thru analog signal equal to actual # of trials
num_trials = height(tbl); 
num_detected_edges = length(timestamps); 
assert(num_trials == num_detected_edges,'detected edges not equal to num trials');

tbl.TimeAnalog = timestamps'; 

%% save/export epoch table 
 
disp('saving file') 
outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/EXP/15s_stim/Epochs',...
    PatientID);
mkdir(outputdir);   %create the directory
thisfile = sprintf('epoch_%s_%s_%s_f%s.mat', subjectname, DBStarget,...
    block_num,freq); 
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
save(fulldestination, 'tbl');  %save the file there directory 
disp('saved')
%end 



 