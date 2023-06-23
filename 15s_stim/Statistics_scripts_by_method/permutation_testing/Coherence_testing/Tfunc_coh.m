function[t_stat] = Tfunc_coh(data,shuffled_labels) 

% Input = data matrix that is trial x ch x ch x freq 

% Output = data matrix containing T-statistic values that is ch x ch x freq OR
% ROI x freq 

%% -----------------------START OF CODE------------------------

%perform t-test for data1 vs data2 
num_ch_roi = size(data,2); 
num_freq = size(data,4);

t_stat = compute_Tstatistic_coh(num_ch_roi,num_freq,data,shuffled_labels); 
end 