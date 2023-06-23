%% FUNCTION: Generate PSD values from multitaper fxn 
%%% Input: PatientID, experiment_name, filename, data (epoched time series), 
%%% sampling freq, good channel labels and metadata 

%%% Output: Spectral power, frequencies of spectral power for each stim
%%% state


function[pow_prestim,pow_stim1,pow_stim2,pow_stim3,...
        pow_poststim,f_prestim, S_prestim, f_stim1 S_stim1,...
        S_stim2,S_stim3,f_poststim,S_poststim,Serr_prestim, Serr_stim1, Serr_stim2,...
        Serr_stim3, Serr_poststim,metadata] = mtspectrum_pow_from_timeseries_baseline(PatientID,...
        experiment_name,timeseriesfile, timeseries_data, fs, ch_labels,...
                data,metadata); 
%% 
num_trials = size(data,1); 
num_ch = size(data, 2);
num_samples = size(data, 3); 
%% Define the different time windows of interest, in samples.

pre_stim_win = 1:5000; 
stim_win = 5001:20000; 
stim_win1 = 5000:10000; 
stim_win2 = 10000:15000; 
stim_win3 = 15000:20000; 
post_stim_win = 21000:26000; % artifact ends a bit after the 20 second mark 

disp('defined stim windows')

metadata.stim_info.pre_stim_win = pre_stim_win; 
metadata.stim_info.stim_win1 = stim_win1; 
metadata.stim_info.stim_win2 = stim_win2; 
metadata.stim_info.stim_win3 = stim_win3; 
metadata.stim_info.post_stim_win = post_stim_win; 
 %% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data = data(:,:,pre_stim_win); 
 indivtr_stim1_data = data(:,:,stim_win1); 
 indivtr_stim2_data = data(:,:,stim_win2); 
 indivtr_stim3_data = data(:,:,stim_win3); 

 indivtr_poststim_data= data(:,:,post_stim_win); 
 
 %% compute PSD for each trial 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 
prestim_color = [0.6706    0.1882    0.1882]; 
stim_color = [0.6549    0.3765    0.7098]; 
poststim_color = [0.5020    0.5020    0.5020]; 
type = 'data'; 

metadata.PSDparams = params; 
%% 
S_prestim = []; 
S_stim1 = []; S_poststim = []; 
for i=1:num_ch
    for j = 1:num_trials 
    %prestim 
        [S_prestim(j,i,:),f_prestim,Serr_prestim(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(j,i,:),[3,1,2])),params);
        [S_stim1(j,i,:),f_stim1,Serr_stim1(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_stim1_data(j,i,:),[3,1,2])),params);
        [S_stim2(j,i,:),f_stim2,Serr_stim2(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_stim2_data(j,i,:),[3,1,2])),params);
        [S_stim3(j,i,:),f_stim3,Serr_stim3(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_stim3_data(j,i,:),[3,1,2])),params);  
        [S_poststim(j,i,:),f_poststim,Serr_poststim(j,i,:,:)] = mtspectrumc(squeeze(permute(indivtr_poststim_data(j,i,:),[3,1,2])),params);
    end    
end 
%% convert to dB 

pow_prestim = 10*log(S_prestim);
pow_stim1 = 10*log(S_stim1);
pow_stim2 = 10*log(S_stim2); 
pow_stim3 = 10*log(S_stim3); 
pow_poststim = 10*log(S_poststim);

disp('converted to log')



end 


