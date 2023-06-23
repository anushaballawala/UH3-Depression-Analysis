            
function[t_stat] = Tfunc(data,shuffled_labels) 

% Input = data matrix that is trial x ch x freq 

% Output = data matrix containing T-statistic values that is ch x freq OR
% ROI x freq 

%% -----------------------START OF CODE------------------------

%perform t-test for data1 vs data2 
num_ch_roi = size(data,2); 
num_freq = size(data,3);

t_stat = compute_Tstatistic(num_ch_roi,num_freq,data,shuffled_labels); 



%     switch data_dim
% 
%         case 'Ch'
% 
%             t_stat = compute_Tstatistic(num_ch_roi,num_freq,data,shuffled_labels); 
% 
%         case 'ROI'
%             %compute t-statistic for ROIs 
%             
%             t_stat = compute_Tstatistic(num_ch_roi,num_freq,matrix_ROI,shuffled_labels); 
%             
%             %** add case for left and right matrix ROI 
% 
%     end   
     
end 