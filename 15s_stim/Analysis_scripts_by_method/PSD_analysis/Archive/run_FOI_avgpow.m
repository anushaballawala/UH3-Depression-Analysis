% turn this into fxn 
clc 
clear
addpath(genpath('/Users/anushaallawala/Data'))
PatientID = 'DBSTRD001'; 
FOI = 'high gamma'; 
%% Get freq indices 

dirs = dir(sprintf('/Users/anushaallawala/Data/DBSTRD/%s/PSD/*.mat',PatientID)); 

files =  {}; 
num_files = numel(dirs); 
for i = 1:num_files
    files{i} = load(sprintf('%s/%s',dirs(i).folder,dirs(i).name)); 
end 

%% get freq band info out from PSD & avg over FOI 
clc 

disp(FOI)
FOI_prestim_avg = cell(1,num_files); 
FOI_stim1_avg = cell(1,num_files); 
FOI_stim2_avg = cell(1,num_files);
FOI_stim3_avg = cell(1,num_files);
FOI_poststim_avg = cell(1,num_files);

for i = 1:num_files

[FOI_prestim_avg{i},FOI_stim1_avg{i},FOI_stim2_avg{i},...
    FOI_stim3_avg{i}, FOI_poststim_avg{i}] = get_FOI_idx_avgpow(FOI,...
    files{1, i}.pow_prestim,...
    files{1, i}.pow_stim1, files{1, i}.pow_stim2, ...
    files{1, i}.pow_stim3, files{1, i}.pow_poststim,...
    files{1, i}.f_prestim, files{1, i}.f_stim1,...
    files{1, i}.f_poststim); 
end 

%% Array of trials X ch. 

prestim_avg = vertcat(FOI_prestim_avg{:});
stim1_avg = vertcat(FOI_stim1_avg{:});
stim2_avg = vertcat(FOI_stim2_avg{:});
stim3_avg = vertcat(FOI_stim3_avg{:});
poststim_avg = vertcat(FOI_poststim_avg{:});


%% Labels 

exp_name = {}; 
for i = 1:num_files 
    exp_name{i} = dirs(i).name; 
    exp_name{i} = char(exp_name{i}); 
    exp_name{i} = exp_name{i}(1:end-32); 
end 

exp_name{1} = exp_name{1}(1:end-9); 

%% REpeat this for first 9 trials and remaining 5 trials except for #9

for i = 1:num_files 
    num_trials(i) = size(FOI_poststim_avg{1, i}, 1); 
end 
    
total_num_trials = sum(num_trials); 

%% Generate experiment label for each trial 
    
        for i = 1:num_files 
       
            exp_label{i,1} = repmat({exp_name{i}},num_trials(i),1); 
        end 

exp_label = vertcat(exp_label{:}); 


%% Electrode label 

elec_label = cell(total_num_trials,1); 

for i = 1:total_num_trials
    elec = 'elec\d*';
    [start_idx_elec, end_idx_elec] = regexp(exp_label{i},elec); 

    if isempty(start_idx_elec)
        disp('do nothing')
        elec_label{i, 1} = 'BaselineFix'; 
    elseif ~isempty(start_idx_elec)
        elec_label{i, 1} = extractBetween(exp_label{i},start_idx_elec,end_idx_elec); 
    end 
end 


%% DBS labels 

for i = 1:length(exp_label)
    
braintarget = '(l|r)(VCVS|SCC)';
[start_idx_br, end_idx_br] = regexp(exp_label{i},braintarget);

    if isempty(start_idx_br) 
        disp('do nothing')
        DBS_target{i,1} = 'BaselineFix'; 
    elseif ~isempty(start_idx_br)
        DBS_target{i,1} = extractBetween(exp_label{i},start_idx_br,end_idx_br);
    end 
end 

%% Frequency labels 

freq_label = {}; 
for i = 1:length(exp_label)
    
frequency = 'f\d*';
[start_idx_f, end_idx_f] = regexp(exp_label{i},frequency);

    if isempty(start_idx_f) 
        disp('do nothing')
        freq_label{i,1} = 'BaselineFix'; 
    elseif ~isempty(start_idx_f)
        freq_label{i,1} = extractBetween(exp_label{i},start_idx_f,end_idx_f);
    end 
end 

%% Save 
filename = sprintf('%s_%s_avg_PSD',PatientID, FOI); 
filepath = sprintf('/Users/anushaallawala/Data/DBSTRD/%s',filename); 

save(filepath,'prestim_avg','stim1_avg','stim2_avg','stim3_avg',...
    'poststim_avg','FOI','exp_label','total_num_trials','DBS_target',...
    'freq_label','elec_label','dirs'); 

