%% run_mtspectrum_pow_from_timeseries_baseline.m - Run this to get spectral
% power for baseline recordings of TRD participants 

% Additional info: Doesn't incorporate FOI idx fxn; PLZ review function before using
% again as of 10/15/21 ***generate_ROI fxn might break**** 
% Inputs: frequency bounds for freq band of interest, z-score of length of
% recording 
% Outputs: z-scored spectral power, std of power across trials 
% Dependencies: Chronux toolbox, UH3 github repo 
% Sub-functions: 

% Anusha Allawala, 8/2021

%------------ START OF CODE ---------------% 

function [stim_minus_pre,post_minus_pre,...
    avg_tr_stim_minus_pre,avg_tr_post_minus_pre,...
    tbl_allch_stimpre,tbl_allch_postpre] = calc_norm_POW(lower_freq,upper_freq,f_prestim,...
    f_stim,f_poststim,z_mean,z_std,pow_prestim,...
    pow_stim,pow_poststim,ch_labels,num_trials,PatientID,ROI_labels,contact_config,num_ch)


%prestim 
freq_prestim_idx = find(f_prestim >lower_freq & f_prestim < upper_freq); 
freq_prestim_f = f_prestim(freq_prestim_idx); 
%stim 
freq_stim_idx = find(f_stim >lower_freq & f_stim < upper_freq);
freq_stim_f = f_stim(freq_stim_idx); 
%poststim
freq_poststim_idx = find(f_poststim >lower_freq & f_poststim < upper_freq);
freq_poststim_f = f_poststim(freq_poststim_idx);

%% Get average across FOI 
%fxn handle 
avgFOI = @(pow,FOI_idx)(mean(pow(:,:,FOI_idx),3)); 

FOI_prestim_avg = avgFOI(pow_prestim,freq_prestim_idx); 
FOI_stim_avg = avgFOI(pow_stim,freq_stim_idx); 
FOI_poststim_avg = avgFOI(pow_poststim,freq_poststim_idx); 

%% Get z-score for every trial, for every channel, then compute mean z-score. 

z_score = @(FOI_avg,z_mean,z_std)((FOI_avg - z_mean)./z_std); 
%z_score = @(a,b,c) a;  % TODO REPLACE WITH REAL ZSCOE AFAIN
FOI_prestim_avg_z = z_score(FOI_prestim_avg,z_mean,z_std); 
FOI_stim_avg_z = z_score(FOI_stim_avg,z_mean,z_std); 
FOI_poststim_avg_z = z_score(FOI_poststim_avg,z_mean,z_std); 



%% subtract stim minus pre & post minus pre z-scores 

stim_minus_pre = FOI_stim_avg_z - FOI_prestim_avg_z; 
post_minus_pre = FOI_poststim_avg_z - FOI_prestim_avg_z; 

%% Average the data across trials 

avgtrials = @(pow)(mean(pow,1)); 

avg_tr_stim_minus_pre = avgtrials(stim_minus_pre); 
avg_tr_post_minus_pre = avgtrials(post_minus_pre); 

%% Get SEM and std across trials 

SEMtrials = @(pow,num_trials)(squeeze(std(pow,0,1)))./num_trials;

SEM_tr_stim_minus_pre = SEMtrials(stim_minus_pre, num_trials); 
SEM_tr_post_minus_pre = SEMtrials(post_minus_pre, num_trials); 

STDtrials = @(pow,num_trials)(std(pow,0,1)); 

STD_tr_stim_minus_pre = STDtrials(stim_minus_pre); 
STD_tr_post_minus_pre = STDtrials(post_minus_pre); 

%% Get ROI 

%for averages 
ROI_stim_minus_pre = generate_ROI(avg_tr_stim_minus_pre',PatientID); 
ROI_post_minus_pre = generate_ROI(avg_tr_post_minus_pre',PatientID); 


%for SEM 
ROI_SEM_stim_minus_pre =generate_ROI(SEM_tr_stim_minus_pre',PatientID); 
ROI_SEM_post_minus_pre = generate_ROI(SEM_tr_post_minus_pre',PatientID); 

% for STD 
ROI_STD_stim_minus_pre = generate_ROI(STD_tr_stim_minus_pre',PatientID); 
ROI_STD_post_minus_pre = generate_ROI(STD_tr_post_minus_pre',PatientID); 
%% Create table with ROI values and for indiv channel 
% num_ROI = numel(ROI_labels); 
% tmp = 1:num_ROI; 
% tmp = tmp'; 
% %create table with ROI values for stim - pre 
% table_as_matrix_ROI_stimpre = horzcat(tmp,ROI_stim_minus_pre,ROI_SEM_stim_minus_pre,...
%     ROI_STD_stim_minus_pre); 
% 
% tbl_ROI_stimpre = array2table(table_as_matrix_ROI_stimpre,'VariableNames',{'ROI_labels',sprintf('stim_pre_%s',contact_config),sprintf('SEM_%s',contact_config),sprintf('STD_%s',contact_config)}); 
% tbl_ROI_stimpre.ROI_labels = ROI_labels'; 
% 
% %create table for post-pre 
% table_as_matrix_ROI_postpre = horzcat(tmp,ROI_post_minus_pre,ROI_SEM_post_minus_pre,...
%     ROI_STD_post_minus_pre); 
% 
% tbl_ROI_postpre = array2table(table_as_matrix_ROI_postpre,'VariableNames',{'ROI_labels',sprintf('post_pre_%s',contact_config),sprintf('SEM_%s',contact_config),sprintf('STD_%s',contact_config)}); 
% tbl_ROI_postpre.ROI_labels = ROI_labels'; 

%% Create tables for individual channels 

%stim - pre 
tmp = 1:num_ch;  
table_allch_stimpre = horzcat(tmp',avg_tr_stim_minus_pre',SEM_tr_stim_minus_pre',STD_tr_stim_minus_pre'); 
tbl_allch_stimpre = array2table(table_allch_stimpre,'VariableNames',{'Ch_name',sprintf('stim_pre_%s',contact_config),...
    sprintf('SEM_%s',contact_config),sprintf('STD_%s',contact_config)}); 
tbl_allch_stimpre.Ch_name = ch_labels; 

%post - pre 
table_allch_postpre = horzcat(tmp',avg_tr_post_minus_pre',SEM_tr_post_minus_pre',STD_tr_post_minus_pre'); 
tbl_allch_postpre = array2table(table_allch_postpre,'VariableNames',{'Ch_name',sprintf('post_pre_%s',contact_config),sprintf('SEM_%s',contact_config),sprintf('STD_%s',contact_config)}); 
tbl_allch_postpre.Ch_name = ch_labels; 