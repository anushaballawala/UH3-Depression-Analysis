function [output] = compute_FOI(FOI,data,f)

%% Define FOI freq band and input 
switch FOI 
    case 'delta'
        lower_freq = 1; 
        upper_freq = 4; 
    case 'theta'
        lower_freq = 4; 
        upper_freq = 7; 
    case 'alpha'
        lower_freq = 8; 
        upper_freq = 13; 
    case 'beta'
        lower_freq = 13; 
        upper_freq = 30; 
    case 'gamma' 
        lower_freq = 50;
        upper_freq = 80;
end

%% Get FOI indices 

freq_idx = find(f >lower_freq & f < upper_freq); 
freq = f(freq_idx);

%% Avg PSD across FOI 

avgFOI = @(pow,FOI_idx)(mean(pow(:,:,FOI_idx),3)); 

avg_FOI = avgFOI(data,freq_idx); 
output = avg_FOI; 
end 