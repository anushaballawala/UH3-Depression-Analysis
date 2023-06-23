%function epochtablestim_15s_003() 

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
cd '/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD003/Raw Data/15s stim/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
cd 
%load .mat file containing comments from Blackrock NSP (originally NEV) 
disp('Select NEV files for 15s STIM')
[filename,filepath]=uigetfile('*.nev','Select NEV File');
NEV = openNEV(fullfile(filepath,filename),'nosave');

%filename = 'sub-TRDDBS001_task-shortstim_run_01_blk-rVCVSf130elec12'; %example filename 
%% 
%get timestamps from comment file 
comments = NEV.Data.Comments; 
timestamps = comments.TimeStampSec';

%get text from comment file
comments_text = string(NEV.Data.Comments.Text);

%% parse out stim param information from stimulation comments and remove bad trials 

 header_text = "#StimOn#";
 N = length(comments_text); %number of total trials in this file 
 %take out garbage comments (if theres any that say there is packet loss,
 %or error .
 frequency_field = strings(N,1); %preallocating freq string 

  for i = 1:N
     trial_text = comments_text(i);
     disp(i)
  end 
  
    start_comment_idx = contains(comments_text,"#StimOn#"); 
    start_comment_idx  = find(start_comment_idx); 
    
    end_comment_idx = start_comment_idx + 3 ; 
        
    new_comments = strcat(comments_text(start_comment_idx),comments_text(start_comment_idx+1),...
        comments_text(start_comment_idx+2), comments_text(start_comment_idx+3)); 
  
  %%
  % parse stim param infor from comemnts  
 trashtrials_idx = frequency_field == "NODATA"; %trials without #stimon# header 
 keeptrials_idx = ~trashtrials_idx; 
 clean_comments_text = new_comments; %comments with trash comments thrown out 
 
 
 %set up parsing out electrode info 
 e_pattern = 'electrode=\d*(\s*\d*\s*\d*\s*\d*\s*\d*\s*\d*)?;';
 e_field = clean_comments_text; % Cheat to make a string array of length N by copying.
  
  
 %number of trials after throwing out bad comments 
 num_trials = length(clean_comments_text); 
 
 
 if num_trials == N 
     fprintf ('no comments thrown out') 
 else 
    fprintf('threw out some comments because of error') 
 end
 
 %% Keep good trials 
 
 keeptrials_idx = start_comment_idx;
 frequency_field = frequency_field(keeptrials_idx); 
 %% Get PatientID info 
 
 %get patientID info 
start_idx = strfind(filename, 'DBSTRD'); 
end_idx = start_idx + 8; 
PatientID = extractBetween(filename, start_idx, end_idx);
PatientID = PatientID{1}; 
metadata.general.PatientID = PatientID; 

fprintf('Patient ID: %s\n', PatientID);
%% create strings with information for each parameter (freq,amp,PW,elec #) 

for i = 1:num_trials
    %parse out freq
    frequency_field(i) = 130;
    %now for amplitude
    amp_field(i) = 900;
    % pulsewidth
    pw_field(i) = 180;
    
    %dont forget the electrode number grandma
    [start_idx_e, end_idx_e] = regexp(clean_comments_text(i), e_pattern);
    e_field(i) = extractBetween(clean_comments_text(i), start_idx_e+10, end_idx_e-1);
    e_field(i) = strrep(e_field(i), '  ','_');
    
end

 %% set up remaining info for epoch table 

braintarget = '(VCVS|SCC)';
[start_idx_br, end_idx_br] = regexp(filename,braintarget);
DBStarget = extractBetween(filename,start_idx_br,end_idx_br);


%create Mx1 variable with condition names for each trial 
condition = cell(num_trials,1);
for i = 1:num_trials 
tmp = sprintf('%s_elec%s_amp%s_freq%s_pw%d' ,DBStarget{1}, e_field(i), amp_field(i), frequency_field(i), pw_field(i)); 
condition{i} = tmp; 
end

%get additional information about run/block/trials 
subjectname = PatientID; 

block = 'blk-\d*';
run = 'run-\w*\d';

[start_idx, end_idx] = regexp(filename,  block);
block_num = extractBetween(filename, start_idx, end_idx);
block_num = block_num{1}; 

%extract run num info 
[start_idx, end_idx] = regexp(filename,  run);
run_num = extractBetween(filename, start_idx, end_idx);
run_num = run_num{1}; 

%generate epoch table 
table_as_matrix = zeros(num_trials, 10);
table_as_matrix(:, 1) = 1:num_trials;
%table_as_matrix(:, 2) = block_name; can't add string yet 
table_as_matrix(:, 3) = timestamps(keeptrials_idx);
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
block_numbers = repmat(block_num,num_trials,1);
tbl.Block = block_numbers; 
tbl.Condition = string(condition);
tbl.Frequency = frequency_field;
tbl.PW = pw_field';
tbl.Amp = amp_field';
tbl.Contacts = e_field;
tbl.Stimtarget = repelem(string(DBStarget{1}),num_trials)';
%% run timestamps fxn 

stim_info = strcat(run_num(5:end),'_',block_num); 

[timestamps,cerestim,stimsync] = extract_timestamps_analog(PatientID,stim_info); 

% check edges detected thru analog signal equal to actual # of trials
num_trials = height(tbl); 
num_detected_edges = length(timestamps); 
assert(num_trials == num_detected_edges,'detected edges not equal to num trials');

tbl.TimeAnalog = timestamps'; 

%% save/export epoch table 
 
disp('saving file') 
outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Epochs',PatientID);
mkdir(outputdir);   %create the directory
thisfile = sprintf('epoch_%s_%s.mat', PatientID, stim_info); 
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
save(fulldestination, 'tbl');  %save the file there directory 

%end 