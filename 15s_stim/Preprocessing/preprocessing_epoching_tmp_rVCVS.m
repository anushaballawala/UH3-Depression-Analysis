
% parfor i = 1:5 
%     preprocessing_pipeline_15s_stim(sprintf('0%d',i)); 
%     
% end 

%%  Load data -- this whole script only generated and saved 130 Hz data
% 50 Hz and 6 Hz still need to be generated 
% script also needs to be cleaned up a lot 
% AA - 10/25/21 

rVCVS_dir = dir('E:/DBSTRD/DBSTRD002/Experiments/15s_stim/Epoched Data/rVCVS_*/*.mat'); 
%%  load all datafiles for rVCVS 002 
datafiles = {}; 
for i = 1:numel(rVCVS_dir)
    datafiles{i} = load(rVCVS_dir(i).name); 
    
end 
    
%% 
data_cell_1 = {}; data_cell_2 = {}; 
stim_cond_cell_1 = {}; stim_cond_cell_2 = {}; stim_cond_cell_3 = {}; 
stim_cond_cell_4 = {}; stim_cond_cell_5 = {}; 


for i = 1:27
    
data_cell_1{i} = datafiles{1, 1}.output_and_metadata{1, i}{1, 1}  ;
stim_cond_cell_1{i} = datafiles{1, 1}.output_and_metadata{1, i}{1, 2}  ;

data_cell_2{i} = datafiles{1, 2}.output_and_metadata{1, i}{1, 1}  ;
stim_cond_cell_2{i} = datafiles{1, 2}.output_and_metadata{1, i}{1, 2}  ;

data_cell_3{i} = datafiles{1, 3}.output_and_metadata{1, i}{1, 1}  ;
stim_cond_cell_3{i} = datafiles{1, 3}.output_and_metadata{1, i}{1, 2}  ;

data_cell_4{i} = datafiles{1, 4}.output_and_metadata{1, i}{1, 1}  ;
stim_cond_cell_4{i} = datafiles{1, 4}.output_and_metadata{1, i}{1, 2}  ;

data_cell_5{i} = datafiles{1, 5}.output_and_metadata{1, i}{1, 1}  ;
stim_cond_cell_5{i} = datafiles{1, 5}.output_and_metadata{1, i}{1, 2}  ;

end 
    
%% remove extra trial for cell 3  

data_cell_3{1,1} = data_cell_3{1, 1}(2,:,:); 

%% concatenate 

alldata_cell = vertcat(data_cell_1,...
    data_cell_2,data_cell_3,data_cell_4,data_cell_5); 
allstim_cell = horzcat(stim_cond_cell_1,...
    stim_cond_cell_2,stim_cond_cell_3,stim_cond_cell_4,stim_cond_cell_5);

alldata_cell2 = vertcat(alldata_cell{:}); 
allstim_cell2 = vertcat(allstim_cell{:}); 

%% convert cell to matrix 

%data = vertcat(alldata_cell2{:}); 
data = alldata_cell2; 
stim_cond = allstim_cell2; 
%stim_cond = string(stim_cond_cell);

%% now group by condition 

%first by frequency 

f130_idx = contains(stim_cond,'freq130') == 1; 
f130_cond = stim_cond(f130_idx); %idx conditions for this freq.%*SAVE THIS
f130_data = data(f130_idx,:,:); %*SAVE THIS

f50_idx = contains(stim_cond,'freq50') == 1; 
f50_cond = stim_cond(f50_idx); %SAVE THIS*
f50_data = data(f50_idx,:,:); %SAVE THIS* 

f6_idx = contains(stim_cond,'freq6') == 1; 
f6_cond = stim_cond(f6_idx); %SAVE THIS*
f6_data = data(f6_idx,:,:); %SAVE THIS* 


%% now keep those with pulsewidth in separate data groups 

pw_100_idx_f130 = contains(f130_cond,'pw100'); 
pw_100_data_f130 = f130_data(pw_100_idx_f130,:,:); %* SAVE THIS
pw_100_cond_f130 = f130_cond(pw_100_idx_f130); %*SAVE THIS 

%% 

pw_180_idx_f130 = contains(f130_cond,'pw180'); 
pw_180_data_f130 = f130_data(pw_180_idx_f130,:,:); %*SAVE THIS
pw_180_cond_f130 = f130_cond(pw_180_idx_f130); %*SAVE THIS 

%% now create a cell that matches normally epoched data for other conditions.
    trial_length = size(pw_180_data_f130,3); 
    num_channels = size(pw_180_data_f130,2); 
    conditions = unique(pw_180_cond_f130); % An array of the unique condition summaries in the column.
    condition_bytrial = pw_180_cond_f130;
    num_conditions = length(conditions); %num of stim conditions tested in this file 
    num_trials = 5; 
%% 
    output = cell(num_conditions,1); % each cell in cell array is data for one condition
    output_and_metadata = {}; 
    tbl_idx = {}; 

    %% put conditions in respective cells 
        
    cond1 = [];%zeros(num_trials,129,30000); 
    cond2 = [];%zeros(num_trials,129,30000); 
    cond3 = [];%zeros(num_trials,129,30000); 
    cond4 = [];%zeros(num_trials,129,30000); 
    cond5 = [];%zeros(num_trials,129,30000); 
    cond6 = [];%zeros(num_trials,129,30000); 
    cond7 = [];%zeros(num_trials,129,30000); 
    
    for i = 1:size(pw_180_data_f130,1)
        %for j = 1:num_trials
            switch pw_180_cond_f130(i)
                
                case conditions(1)
                    cond1 = vertcat(cond1, pw_180_data_f130(i,:,:));
                    %cond1(j,:,:) = pw_180_data_f130(i,:,:) ;
                case conditions(2)
                    cond2 = vertcat(cond2, pw_180_data_f130(i,:,:));
                    %cond2(j,:,:) = pw_180_data_f130(i,:,:);
                case conditions(3)
                    cond3 = vertcat(cond3, pw_180_data_f130(i,:,:));
                    %cond3(j,:,:) = pw_180_data_f130(i,:,:);
                case conditions(4)
                    cond4 = vertcat(cond4, pw_180_data_f130(i,:,:));
                    %cond4(j,:,:) = pw_180_data_f130(i,:,:);
                case conditions(5)
                    cond5 = vertcat(cond5, pw_180_data_f130(i,:,:));
                    %cond5(j,:,:) = pw_180_data_f130(i,:,:);
                case conditions(6)
                    cond6 = vertcat(cond6, pw_180_data_f130(i,:,:));
                    %cond6(j,:,:) = pw_180_data_f130(i,:,:);
                case conditions(7)
                    cond7 = vertcat(cond7, pw_180_data_f130(i,:,:));
                    %cond7(j,:,:) = pw_180_data_f130(i,:,:);
            %end
            
        end
        
    end 
    
    output = {cond1;cond2;cond3;cond4;cond5;cond6;cond7}; 
    
    
%% make output cell 

output_and_metadata = cell(1,numel(conditions)); 

for i = 1:numel(conditions) 
    
    output_and_metadata{1,i}{1,1} = output{i};  
    output_and_metadata{1,i}{1,2} = conditions(i); 


end 

    
%% metadata

metadata.metadata_file1 = datafiles{1, 1}.metadata;  
metadata.metadata_file2 = datafiles{1, 2}.metadata ; 
metadata.metadata_file3 = datafiles{1, 3}.metadata  ;
metadata.metadata_file4 = datafiles{1, 4}.metadata ; 
metadata.metadata_file5 = datafiles{1, 5}.metadata  ;

metadata.epoched.tbl_idx = tbl_idx; 
metadata.epoched.conditions = conditions; 

%% save
stim_info = 'rVCVS_f130'; 
PatientID = 'DBSTRD002'; 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s',PatientID,stim_info);
%outputdir = pwd; 
if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end
disp(outputdir)
addpath(genpath(outputdir)); 
%% save indiv trial data all in one file **
thisfile = sprintf('15s_stim_all_currdir_timeseries_singletrial_%s.mat', stim_info);
fulldestination = fullfile(outputdir, thisfile);  %name file relative to that directory
disp('saving.....')
save(fulldestination,'output_and_metadata','metadata','-v7.3'); 
disp('saved')



























    
    



