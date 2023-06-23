%% get_FOIidx_avgpow.m - Handy fxn that can give you indices of frequency bands
% so it doesn't have to be incorporated into your script everytime. 
% Also averages power in each freq band 

% Additional info: N/A
% Inputs: Freq band of interest, power in different stim states. Power data
% can be structured as time x ch x frequency OR ch x tr x freq. Only req is
% THIRD DIMENSION IS FREQUENCY 
% Outputs: Avg power in freq of interest 
% Dependencies: UH3 github repo 
% Sub-functions: N/A 

% Anusha Allawala, 10/2021

%------------ START OF CODE ---------------% 

function [FOI_prestim_avg,FOI_stim1_avg,FOI_stim2_avg,...
    FOI_stim3_avg, FOI_poststim_avg] = get_FOI_idx_avgpow(FOI,...
    pow_prestim,pow_stim1,pow_stim2,pow_stim3,pow_poststim,...
    f_prestim,f_stim,f_poststim)
switch FOI 
    case 'delta'
        lower_freq = 1; 
        upper_freq = 4; 
    case 'alpha'
        lower_freq = 8; 
        upper_freq = 13; 
    case 'theta'
        lower_freq = 4; 
        upper_freq = 7; 
    case 'beta'       
        lower_freq = 13; 
        upper_freq = 30; 
    case 'low gamma' 
        lower_freq = 40; 
        upper_freq = 65; 
    case 'high gamma' 
        lower_freq = 70; 
        upper_freq = 90; 
end 


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
FOI_stim1_avg = avgFOI(pow_stim1,freq_stim_idx); 
FOI_stim2_avg = avgFOI(pow_stim2,freq_stim_idx); 
FOI_stim3_avg = avgFOI(pow_stim3,freq_stim_idx); 
FOI_poststim_avg = avgFOI(pow_poststim,freq_poststim_idx); 
disp('......Averaging data in FOI')


end 