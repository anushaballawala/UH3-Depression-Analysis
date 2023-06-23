function [dirs_for_job] = run_coherence_fxn(job_idx, num_jobs)

    dirs_for_job = get_files_for_job(job_idx, num_jobs);
    
    % Print info
    fprintf("Job %d of %d is working on these %d dirs:\n", job_idx, num_jobs, length(dirs_for_job)); 
    for dir_idx = 1:length(dirs_for_job)
        fprintf("\t%s\n", dirs_for_job(dir_idx));
        run_compute_coherence_oscar_poststim_DBStarget(dirs_for_job(dir_idx))
    end
  

    
end