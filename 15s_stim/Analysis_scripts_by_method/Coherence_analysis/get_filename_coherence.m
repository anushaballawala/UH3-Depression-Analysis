function [dirs_for_job] = get_filename_coherence(job_idx,num_jobs)

addpath(genpath('/gpfs/data/dborton/TRD_Project/DBSTRD/DBSTRD002/Experiments/15s_stim')); 
addpath(genpath('/gpfs/data/dborton/TRD_Project/DBSTRD/DBSTRD001/Experiments/15s_stim')); 

addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/CODE/UH3-Depression-Preprocessing')); 

    dirs_for_job = get_dirs_for_job(job_idx, num_jobs);
    disp('now performing for DBSTRD001 SCC only') 
    % Print info
    fprintf("Job %d of %d is working on these %d dirs:\n", job_idx, num_jobs, length(dirs_for_job)); 
    for dir_idx = 1:length(dirs_for_job)
        fprintf("\t%s\n", dirs_for_job(dir_idx));
        run_compute_coherence_oscar(dirs_for_job(dir_idx))
        run_compute_coherence_oscar_stim2(dirs_for_job(dir_idx))
        run_compute_coherence_oscar_stim3(dirs_for_job(dir_idx))
        run_compute_coherence_oscar_poststim(dirs_for_job(dir_idx))
        run_compute_coherence_oscar_prestim(dirs_for_job(dir_idx))
        disp(char(dirs_for_job(dir_idx)))
    end

end 

