[pow,f, S_prestim, Serr_poststim,metadata] = mtspectrum_pow_from_timeseries_baseline_autocorr(PatientID,...
    experiment_name, timeseriesfile, timeseries_data, fs, ch_labels,...
    data,metadata);%% FUNCTION: Generate PSD values from multitaper fxn 
%%% Input: PatientID, experiment_name, filename, data (epoched time series), 
%%% sampling freq, good channel labels and metadata 

%%% Output: Spectral power, frequencies of spectral power for each stim
%%% state


function[pow,f, S_prestim, Serr_poststim,metadata] = mtspectrum_pow_from_timeseries_baseline(PatientID,...
        experiment_name,timeseriesfile, timeseries_data, fs, ch_labels,...
                data,metadata); 
%% 
num_trials = size(data,1); 
num_ch = size(data, 2);
num_samples = size(data, 3); 
%% Define the different time windows of interest, in samples.


 %% Get different stim windows with indiv trial data 


 indivtr_poststim_data= data(:,:,post_stim_win); 
 
 %% compute PSD for each trial 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
S_poststim = []; 
for i=1:num_ch
    for j = 1:num_trials 
        [S_poststim(j,i,:),f_poststim,Serr_poststim(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_poststim_data(j,i,:),[3,1,2])),params);
    end    
end 
%% convert to dB 

pow = 10*log(S_poststim);

disp('converted to log')



end 


