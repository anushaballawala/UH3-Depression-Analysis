
 
%% Get information about data size & metadata
num_trials = size(data,1); 
num_ch = size(data, 2);
num_samples = size(data, 3); 
metadata.processing.num_trials = num_trials ; 
metadata.processing.num_ch = num_ch; 
metadata.processing.num_samples = num_samples; 
%% Define window of interest. 
prestim = 'full window'; 
stim = 'full window'; 
poststim = 'full window'; 

[pre_stim_win, stim_win, post_stim_win] = define_window_of_interest(PatientID, prestim, stim, poststim); 
 
% add to metadata 
metadata.processing.fake_pre_stim_window = pre_stim_win; 
metadata.processing.fake_stim_window = stim_win; 
metadata.processing.fake_post_stim_window = post_stim_win;
metadata.processing.stim_window_notes = '5 second duration 1-second post stim picked as post-stim window'; 

%% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data = data(:,:,pre_stim_win); 
 indivtr_stim_data = data(:,:,stim_win); 
 indivtr_poststim_total_data = data(:,:,post_stim_win); 
 
 disp('Size of pre-stim window')
 disp(size(indivtr_pre_stim_data))
 disp('Size of stim window') 
 disp(size(indivtr_stim_data))
 disp('size of post-stim window')
 disp(size(indivtr_poststim_total_data))
 %%  Set up parameters for multitaper spectrum fxn. 
params.Fs = 1000;
params.fpass = [1 100]; % 1 hz to 40 Hz 
disp('Bandpass')
disp(params.fpass)
params.err = [1 0.05];
params.trialave = 0; % change to 1 if you want to average over trials. 
params.tapers = [4 7];
disp(params.tapers)

% add to metadata 

metadata.processing.multitaperspec.params = params; 
%% Use Chronux function to run multitaper spectrum 
tic 
S_prestim = []; 
S_stim = []; 
S_poststim = []; 
for i=1:num_ch
    fprintf('....Processing Channel %d \n',i)
    for j = 1:num_trials 
    %prestim 
    fprintf('....Processing Trial %d \n',j)
        [S_prestim(j,i,:),f_prestim,~] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(j,i,:),[3,1,2])),params);
        [S_stim(j,i,:),f_stim,~] = mtspectrumc(squeeze(permute(indivtr_stim_data(j,i,:),[3,1,2])),params);
        [S_poststim(j,i,:),f_poststim,~] = mtspectrumc(squeeze(permute(indivtr_poststim_total_data(j,i,:),[3,1,2])),params);
    end    
end
disp('Multitaper spectrum power computation complete') 
toc 


%% Get z-score for every trial, for every channel, then compute mean z-score. 
% get mean 
mean_bl = @(data)(mean(mean(data,1),3)); 
z_score_prestim = mean_bl(S_prestim); 
z_score_poststim = mean_bl(S_poststim); 
z_score_stim = mean_bl(S_stim); 

full_recording = vertcat(z_score_prestim,z_score_poststim,z_score_stim);%for entire recording 
mean_baseline = mean(full_recording,1); 

% get std 
std_bl = @(data)(std(data,0,3)); 
std_bl_2 = @(data)(std(data,0,1)); 
std_prestim = std_bl(S_prestim); std_prestim = std_bl_2(std_prestim); 
std_poststim = std_bl(S_poststim); std_poststim = std_bl_2(std_poststim);
std_stim = std_bl(S_stim); std_stim = std_bl_2(std_stim);
full_recording_std = vertcat(std_prestim,std_stim,std_poststim); 
std_baseline = std(full_recording_std,0,1); 


disp('......z-scoring data')

%Convert data 
z_score = @(avg,z_mean,z_std)((avg - z_mean)./z_std); 

prestim_z = z_score(S_prestim,mean_baseline,std_baseline); 
stim_z = z_score(S_stim,mean_baseline,std_baseline); 
poststim_z = z_score(S_poststim,mean_baseline,std_baseline); 



%% Convert to dB 

pow_prestim = abs(10*log(prestim_z));
pow_stim = abs(10*log(stim_z));
pow_poststim = abs(10*log(poststim_z));
disp('.....Converting to log dB')

%% Get power in FOI 
% Get FOI indices and avg power for FOI 

%delta 
FOI = 'delta'; 

[FOI_prestim_avg_delta,FOI_stim_avg_delta,FOI_poststim_avg_delta] = get_FOI_idx_power(FOI,...
    pow_prestim,pow_stim,pow_poststim,f_prestim,f_stim,f_poststim) ; 
% theta 
FOI = 'theta'; 
[FOI_prestim_avg_theta,FOI_stim_avg_theta,FOI_poststim_avg_theta] = get_FOI_idx_power(FOI,...
    pow_prestim,pow_stim,pow_poststim,f_prestim,f_stim,f_poststim) ; 
%alpha 
FOI = 'alpha'; 
[FOI_prestim_avg_alpha,FOI_stim_avg_alpha,FOI_poststim_avg_alpha] = get_FOI_idx_power(FOI,...
    pow_prestim,pow_stim,pow_poststim,f_prestim,f_stim,f_poststim) ; 

%beta 
FOI = 'beta'; 
[FOI_prestim_avg_beta,FOI_stim_avg_beta,FOI_poststim_avg_beta] = get_FOI_idx_power(FOI,...
    pow_prestim,pow_stim,pow_poststim,f_prestim,f_stim,f_poststim) ; 

%gamma 
FOI = 'gamma'; 
[FOI_prestim_avg_gamma,FOI_stim_avg_gamma,FOI_poststim_avg_gamma] = get_FOI_idx_power(FOI,...
    pow_prestim,pow_stim,pow_poststim,f_prestim,f_stim,f_poststim) ; 

disp('.....Getting data in FOI')

%% Get ROI 
% for individual trials for each FOI

%delta 
ROI_prestim_delta = generate_ROI_mtspectrum_pow(FOI_prestim_avg_delta,PatientID);
ROI_stim_delta = generate_ROI_mtspectrum_pow(FOI_stim_avg_delta,PatientID);
ROI_poststim_delta = generate_ROI_mtspectrum_pow(FOI_poststim_avg_delta,PatientID); 

%theta 
ROI_prestim_theta = generate_ROI_mtspectrum_pow(FOI_prestim_avg_theta,PatientID);
ROI_stim_theta = generate_ROI_mtspectrum_pow(FOI_stim_avg_theta,PatientID);
ROI_poststim_theta = generate_ROI_mtspectrum_pow(FOI_poststim_avg_theta,PatientID); 

%alpha 
ROI_prestim_alpha = generate_ROI_mtspectrum_pow(FOI_prestim_avg_alpha,PatientID);
ROI_stim_alpha = generate_ROI_mtspectrum_pow(FOI_stim_avg_alpha,PatientID);
ROI_poststim_alpha = generate_ROI_mtspectrum_pow(FOI_poststim_avg_alpha,PatientID); 

%beta 
ROI_prestim_beta = generate_ROI_mtspectrum_pow(FOI_prestim_avg_beta,PatientID);
ROI_stim_beta = generate_ROI_mtspectrum_pow(FOI_stim_avg_beta,PatientID);
ROI_poststim_beta = generate_ROI_mtspectrum_pow(FOI_poststim_avg_beta,PatientID); 

%gamma 
ROI_prestim_gamma = generate_ROI_mtspectrum_pow(FOI_prestim_avg_gamma,PatientID);
ROI_stim_gamma = generate_ROI_mtspectrum_pow(FOI_stim_avg_gamma,PatientID);
ROI_poststim_gamma = generate_ROI_mtspectrum_pow(FOI_poststim_avg_gamma,PatientID); 

% save data 


