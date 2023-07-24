function [] = run_15s_preproc(job_idx,num_jobs)   

    dirs_for_job = get_files_for_job(job_idx, num_jobs);
    disp('retrieved directories')
    % Print info
    fprintf("Job %d of %d is working on these %d files:\n", job_idx, num_jobs, length(dirs_for_job));
    for dir_idx = 1:length(dirs_for_job)
        fprintf("\t%s\n", dirs_for_job(dir_idx));
        Preprocessing_setup_15s_stim_DBSTRD008(dirs_for_job(dir_idx))
    end 
    
end
