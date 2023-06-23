function[metadata,tbl_allch_stimpre_theta, tbl_allch_postpre_theta,... %end theta 
    tbl_allch_stimpre_alpha, tbl_allch_postpre_alpha,...
    tbl_allch_stimpre_beta, tbl_allch_postpre_beta,...
    tbl_allch_stimpre_delta, tbl_allch_postpre_delta] = generate_heatmap_mtspectrumc_power(PatientID,metadata,data,z_mean,z_std,...
     ROI_labels,contact_config,ch_labels)

all_FOI = {'theta','alpha','delta','beta'}; 
%% Get information about data size & metadata
num_trials = size(data,1); 
num_ch = size(data, 2);
num_samples = size(data, 3); 
metadata.processing.num_trials = num_trials ; 
metadata.processing.num_ch = num_ch; 
metadata.processing.num_samples = num_samples; 
metadata.processing.contact_config = contact_config; 
%% Define window of interest. 
prestim = 'full window'; 
stim = 'full window'; 
poststim = 'full window'; 

[pre_stim_win, stim_win, post_stim_win] = define_window_of_interest(PatientID, prestim, stim, poststim); 

% add to metadata 
metadata.processing.pre_stim_window = pre_stim_win; 
metadata.processing.stim_window = stim_win; 
metadata.processing.post_stim_window = post_stim_win;
metadata.processing.stim_window_notes = '5 second duration 1-second post stim picked as post-stim window'; 

%% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data = data(:,:,pre_stim_win); 
 indivtr_stim_data = data(:,:,stim_win); 
 indivtr_poststim_data = data(:,:,post_stim_win); 
 
 %%  Set up parameters for multitaper spectrum fxn. 
params.Fs = 1000;
params.fpass = [1 40];
params.err = [1 0.05];
params.trialave = 0; % change to 1 if you want to average over trials. 
params.tapers = [4 7]; 

% add to metadata 

metadata.processing.multitaperspec.params = params; 
%% Use Chronux function to run multitaper spectrum 
[S_prestim,f_prestim,...
    S_stim,f_stim, S_poststim, f_poststim] = compute_multitaper_spectrum_power(num_ch,params,...
    num_trials,indivtr_pre_stim_data,...
    indivtr_stim_data, indivtr_poststim_data); 

%% Convert to dB 

pow_prestim = 10*log(S_prestim);
pow_stim = 10*log(S_stim);
pow_poststim = 10*log(S_poststim);

%% Get power in FOI 
% Get FOI indices 
for i = 1:numel(all_FOI) 
    FOI = all_FOI{i} ; 
    
    switch FOI 
        case 'theta' 
            lower_freq = 4; 
            upper_freq = 7; 
            [avg_tr_stim_minus_pre_theta,avg_tr_post_minus_pre_theta,...
            tbl_allch_stimpre_theta,tbl_allch_postpre_theta] = calc_norm_POW(lower_freq,upper_freq,f_prestim,...
            f_stim,f_poststim,z_mean,z_std,pow_prestim,...
            pow_stim,pow_poststim,ch_labels,num_trials,PatientID,ROI_labels,contact_config,num_ch);

        case 'alpha'
            lower_freq = 8;
            upper_freq = 13;
            [avg_tr_stim_minus_pre_alpha,avg_tr_post_minus_pre_alpha,...
            tbl_allch_stimpre_alpha,tbl_allch_postpre_alpha] = calc_norm_POW(lower_freq,upper_freq,f_prestim,...
            f_stim,f_poststim,z_mean,z_std,pow_prestim,...
            pow_stim,pow_poststim,ch_labels,num_trials,PatientID,ROI_labels,contact_config,num_ch);
        case 'beta'
            lower_freq = 13;
            upper_freq = 30;
            [avg_tr_stim_minus_pre_beta,avg_tr_post_minus_pre_beta,...
            tbl_allch_stimpre_beta,tbl_allch_postpre_beta] = calc_norm_POW(lower_freq,upper_freq,f_prestim,...
            f_stim,f_poststim,z_mean,z_std,pow_prestim,...
            pow_stim,pow_poststim,ch_labels,num_trials,PatientID,ROI_labels,contact_config,num_ch);
        case 'delta'
            lower_freq = 1;
            upper_freq = 4;
            [avg_tr_stim_minus_pre_delta,avg_tr_post_minus_pre_delta,...
            tbl_allch_stimpre_delta,tbl_allch_postpre_delta] = calc_norm_POW(lower_freq,upper_freq,f_prestim,...
            f_stim,f_poststim,z_mean,z_std,pow_prestim,...
            pow_stim,pow_poststim,ch_labels,num_trials,PatientID,ROI_labels,contact_config,num_ch);
    end
    
end 
   

end 