
function [] = plot_time_series_channel(t, data , ...
    num_ch, PatientID,contact_config, fs, ch_labels, experiment_name,outputdir) 

% num_samples = size(data,3); 
% num_ch = size(data,2); 
% num_trials = size(data,1); 
% t = (0:(num_samples-1))./fs; 


% If less than 5 trials 
if size(data,1) == 4  
    num_trials = 4 ; 
elseif size(data,1) == 5  
    num_trials = 5; 
end 

parfor i = 1:num_ch 
    figure()
    for j = 1:num_trials         
        subplot(num_trials,1,j)
        plot(t,squeeze(data(j,i,:)))
        ylim([-120 120])
        title(sprintf('%s',ch_labels(i)))
    end    
    %filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/Time-Series Plots/%s/timeseries_%s_%s.fig',PatientID,contact_config,ch_labels(i),experiment_name); 
    name = sprintf('%s/%s.fig',outputdir,ch_labels(i)); 
    saveas (gcf,name); 
    %filename_png = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/Time-Series Plots/%s/timeseries_%s_%s.png',PatientID, contact_config,ch_labels(i),experiment_name); 
    name_png = sprintf('%s/%s.png',outputdir,ch_labels(i)); 
    saveas (gcf,name_png); 
    
end 