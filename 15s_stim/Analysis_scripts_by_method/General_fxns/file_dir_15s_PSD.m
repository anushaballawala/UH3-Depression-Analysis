% make files for mtspectrum power 


directories = dir('/gpfs/data/dborton/TRD_Project/DBSTRD/*/Experiments/15s_stim/Epoched Data/*/15s_stim_all_currdir_timeseries_singletrial_*.mat'); 
%% 
for i = 1:8 
    timeseriesfile{i} = sprintf('%s/%s',directories(i).folder,directories(i).name); 
end 


%% 






