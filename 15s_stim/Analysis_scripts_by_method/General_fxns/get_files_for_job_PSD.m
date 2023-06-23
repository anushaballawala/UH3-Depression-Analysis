function [files_for_job] = get_files_for_job_PSD(job_idx,num_jobs) 

    %files = load('29-Sep-2021_list_of_PCCT_subjects_ChannelCheck_test.mat');
    files = load('11-Nov-2021_list_of_dirs_PSD_15sstim.mat'); 
    files = files.timeseriesfile; 
    disp('Loaded files')
    %which rows of files will belong to jobidx? 
    total_range = 1:length(files); 
    
    % Assign the rows into num_jobs different buckets, like 1,2,3,1,2,3,1,2,...
    % This is like dealing cards around the circle.
    labels = mod(total_range - 1,num_jobs) + 1; 
    %files_for_job = files(labels == job_idx); % go back to this ** 
    files_for_job = string(files{labels == job_idx}); 
    disp('got files for the batch job')
end 