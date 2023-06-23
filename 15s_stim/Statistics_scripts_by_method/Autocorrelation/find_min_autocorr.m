function [first_lag,first_lag_allROI] = find_min_autocorr(matrix_ROI, ROI_labels,PatientID,freqband)
%% Fxn for autocorr 


num_ROI = numel(ROI_labels);
acf = {}; 
lags = {}; 
for i = 1:num_ROI 
    figure() 
    [acf{i},lags{i},bounds,h] = autocorr(double(matrix_ROI(i,:)),'NumLags',304080,'NumSTD',3); 
    title(ROI_labels{i})
    cd /Users/anushaallawala/Data/Baseline_Data/003/autocorr_figs
    saveas(gcf, sprintf('%s_%s_%s.fig',PatientID,freqband,ROI_labels{i}))
    saveas(gcf, sprintf('%s_%s_%s.png',PatientID,freqband,ROI_labels{i}))

end 



%% Find minimum # of lags 
%**** THIS GIVES A ROUGH ESTIMATE FOR WHEN AUTOCORRELATION DROPS CLOSE TO
%ZERO, I MANUALLY CHECKED THE DISTRIBUTION OF VALUES ACROSS ROI'S AND
%FREQUENCIES AND THE GENERATED AUTOCORRELATION PLOTS TO GET A SENSE FOR
%WHETHER THIS MINIMUM VALUE IS ACTUALLY REPRESENTATIVE OF THE # OF LAGS
%WHEN AUTOCORRELATION DROPS TO ZERO *******

for i = 1:num_ROI 
    acf_idx{i} = find(acf{i}<0.005);
    first_lag_idx{i} = lags{1,i}(acf_idx{i});                                                                                                 
    first_lag_allROI{i} = first_lag_idx{i}(1); 
   
end 
first_lag = cell2mat(first_lag_allROI); 

