%%%%% The goal of this script is to generate PSD graphs showing differences
%%%%% between ON vs OFF and generate 1/f trend lines for each stim
%%%%% condition b/w ON vs OFF states
% *** might have to epoch time domain data 
% A) do with no PARRM (ON vs OFF) 
% B) do with PARRM (ON vs OFF) 
% C) do with PARRM (ON vs ON with PARRM, OFF vs OFF with PARRM) 

clear all 
PatientID = 'DBSTRD001'; 
hemi = 'r' ; 
DBS_target = 'SCC'; 
freq = 130; 
stim_freq = freq; 
% Load single trial time series data 
timeseriesfile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
    PatientID,hemi,DBS_target,freq,hemi,DBS_target,freq); 
timeseries_data = load(timeseriesfile); 
fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 
%% Assign data for each condition 

switch DBS_target 
    case 'VCVS'
        [elec1, elec25, elec36, elec47, elec8] = assign_VCVS_conditions(timeseries_data,hemi); 
    case 'SCC'
        [elec1, elec234, elec25, elec36, elec47,...
            elec567,elec8] = assign_SCC_conditions(timeseries_data,hemi); 
end 

%% Describe the data as a sanity check.
disp("The data consists of the power in *physical units* for each bin");
disp("Size of data (trials x channels x freq x time):");
data = elec25; 
current_config = 'elec25'; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 

disp(size(data));
num_trials = size(data, 1);
num_channels = size(data, 2);
num_freq = size(data, 3);
num_samples = size(data, 4);

 %% Plot time-series traces for windows for inspection 
 
 %confirm that there is no artifact during the time windows defined 

num_samples = size(data,3); 
num_ch = size(data,2); 
num_trials = size(data,1); 
t = (0:(num_samples-1))./fs; 


parfor i = 1:num_ch 
    figure()
    for j = 1:num_trials 
        subplot(5,1,j)
    	plot(t,squeeze(data(j,i,:))) 
        title(sprintf('Trial %d %s',j,ch_labels(i))) 
    end   
    
    filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/Time-series_trace/Timeseries_%s.png',PatientID,experiment_name,ch_labels(i)); 
    saveas(gcf, filename)
end 

clear filename 

%% Define the different time windows of interest, in samples.
pre_stim_win = 1:5000; 
stim_win = 5001:20000; 
post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
post_stim_win2 = 25001:30000; 
post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
%% Get different stim windows with data averaged across trials 

avg_data = squeeze(mean(data,1)); 
size(avg_data)

 pre_stim_data = avg_data(:,pre_stim_win); 
 stim_data = avg_data(:,stim_win); 
 poststim1_data = avg_data(:,post_stim_win1); 
 poststim2_data = avg_data(:,post_stim_win2); 
 total_poststim_data = avg_data(:,post_stim_total_win);  
 
 
 %% Get different stim windows with indiv trial data 
 
 pre_stim_data_indivtr = data(:,:,pre_stim_win); 
 stim_data_indivtr = data(:,:,stim_win); 
 poststim_data_total_indivtr = data(:,:,post_stim_total_win); 
 
 %% Generate periodogram for ON vs OFF for all channels 
parfor i = 1:num_ch
    
    figure('Position',[500 200 1500 1000])
    subplot(2,1,1)
    [Pxx,f,mdl,slope] = plot_1overf(pre_stim_data(i,:),fs,'180'); 
    hold on 
    [Pxx_stim,f_stim,mdl_stim,slope_stim] = plot_1overf(stim_data(i,:),fs,'180') 
    [Pxx_poststim,f_poststim,mdl_poststim,slope_poststim] = plot_1overf(total_poststim_data(i,:),fs,'180'); 
    title(ch_labels(i))
    legend('prestim','prestim','prestim','stim','stim','stim','poststim','poststim','poststim','Location','southwest')
    dim_1 = [0.122231270358306 0.926500000491738 0.332247547344199 0.0264999995082617];
    rsqr = mdl.Rsquared.Ordinary; 
    slope = mdl.Coefficients.Estimate(2); 
    str = sprintf('slope pre = %s; slope stim = %s; slope post = %s',...
        num2str(slope), num2str(slope_stim), num2str(slope_poststim)); 
    annotation('textbox',dim_1,'String',str,'FitBoxToText','on'); 
    

    subplot(2,1,2)
    [Pxx,f,mdl,slope] = plot_1overf(pre_stim_data(i,:),fs,'50'); 
    hold on 
    [Pxx_stim,f_stim,mdl_stim,slope_stim] = plot_1overf(stim_data(i,:),fs,'50');
    [Pxx_poststim,f_poststim,mdl_poststim,slope_poststim] = plot_1overf(total_poststim_data(i,:),fs,'50'); 
    title(ch_labels(i))
    rsqr = mdl.Rsquared.Ordinary; 
    slope = mdl.Coefficients.Estimate(2); 
    dim_2 = [0.0749106078665077 0.472500000491738 0.495232405143216 0.0264999995082617];

    str = sprintf('slope pre = %s; slope stim = %s; slope post = %s',...
        num2str(slope), num2str(slope_stim), num2str(slope_poststim)); 
    annotation('textbox',dim_2,'String',str,'FitBoxToText','on'); 
    
    filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/1overf_channels/%s_PSD.png',PatientID,experiment_name,ch_labels(i)); 
    saveas(gcf, filename)

end

%% Generate ROIs and make periodograms
% load ROI 

load('E:/DBSTRD/DBSTRD001/Experiments/ROI_labels_DBSTRD001.mat'); 
%extract ROIs for individual trial data 

[ROI_pre_stim_indivtr] = generate_ROI_indiv_tr(pre_stim_data_indivtr); % gives ROI x current dir matrix 
[ROI_stim_indivtr] = generate_ROI_indiv_tr(stim_data_indivtr); 
[ROI_poststim_total_indivtr] = generate_ROI_indiv_tr(poststim_data_total_indivtr); 

%extract ROIs for averaged data 
[ROI_pre_stim] = generate_ROI(pre_stim_data); 
[ROI_stim] = generate_ROI(stim_data);
[ROI_poststim_total] = generate_ROI(total_poststim_data); 

%
ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG','lVPF',...
    'rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rVPF'}; 
%% Generate ROI periodograms 
    figure('Position',[400,100,2000,800])
    
    slope = zeros(14,1); 
    slope_stim = zeros(14,1);  
    slope_poststim = zeros(14,1); 
    
for i = 1:14%num_ROI 
    
    subplot(2,7,i)
    [Pxx,f,mdl,slope(i)] = plot_1overf(ROI_pre_stim(i,:),fs,'50'); 
    hold on 
    [Pxx_stim,f_stim,mdl_stim,slope_stim(i)] = plot_1overf(ROI_stim(i,:),fs,'50');
    [Pxx_poststim,f_poststim,mdl_poststim,slope_poststim(i)] = plot_1overf(ROI_poststim_total(i,:),fs,'50'); 
    title(ROI_labels(i))

    rsqr = mdl.Rsquared.Ordinary; 
    slope = mdl.Coefficients.Estimate(2); 
    dim_2 = [0.0749106078665077 0.472500000491738 0.495232405143216 0.0264999995082617];

%     str = sprintf('slope pre = %s; slope stim = %s; slope post = %s',...
%         num2str(slope), num2str(slope_stim), num2str(slope_poststim)); 
%     annotation('textbox',dim_2,'String',str,'FitBoxToText','on'); 
    title(ROI_labels(i))
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
end 
    legend('prestim','prestim','prestim','stim','stim','stim','poststim','poststim','poststim','Position',[0.928110457482198 0.468125004526228 0.0512318022357257 0.161874995473772]);
    name = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/1overf_ROI/ROI_1overf50hz_noPARRM.fig',PatientID,experiment_name); 
    saveas(gcf,name); 
    

    outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s',PatientID, experiment_name); 

if ~exist(outputdir,'dir'), mkdir(outputdir), end

thisfile = sprintf('%s_ROI_1overf_slopevalues_50Hz.mat', experiment_name); 
fulldestination = fullfile(outputdir, thisfile); 

save(fulldestination, 'ROI_labels','slope','slope_poststim','slope_stim');


%% Generate ROI periodograms 
    figure('Position',[400,100,2000,800])
    
    slope = zeros(14,1); 
    slope_stim = zeros(14,1);  
    slope_poststim = zeros(14,1); 
    
for i = 1:14%num_ROI 
    
    subplot(2,7,i)
    [Pxx,f,mdl,slope(i)] = plot_1overf(ROI_pre_stim(i,:),fs,'180'); 
    hold on 
    [Pxx_stim,f_stim,mdl_stim,slope_stim(i)] = plot_1overf(ROI_stim(i,:),fs,'180');
    [Pxx_poststim,f_poststim,mdl_poststim,slope_poststim(i)] = plot_1overf(ROI_poststim_total(i,:),fs,'180'); 
    title(ROI_labels(i))

    rsqr = mdl.Rsquared.Ordinary; 
    slope = mdl.Coefficients.Estimate(2); 
    dim_2 = [0.0749106078665077 0.472500000491738 0.495232405143216 0.0264999995082617];

%     str = sprintf('slope pre = %s; slope stim = %s; slope post = %s',...
%         num2str(slope), num2str(slope_stim), num2str(slope_poststim)); 
%     annotation('textbox',dim_2,'String',str,'FitBoxToText','on'); 
    title(ROI_labels(i))
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
end 
    legend('prestim','prestim','prestim','stim','stim','stim','poststim','poststim','poststim','Position',[0.928110457482198 0.468125004526228 0.0512318022357257 0.161874995473772]);
    name = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/1overf_ROI/ROI_1overf180hz_noPARRM.fig',PatientID,experiment_name); 
    saveas(gcf,name); 
    

    outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s',PatientID, experiment_name); 

if ~exist(outputdir,'dir'), mkdir(outputdir), end

thisfile = sprintf('%s_ROI_1overf_slopevalues_180Hz.mat', experiment_name); 
fulldestination = fullfile(outputdir, thisfile); 

save(fulldestination, 'ROI_labels','slope','slope_poststim','slope_stim');






%% look at variability in ROI across trials 
    figure('Position',[400,100,2000,800])

for i = 1:14
    subplot(2,7,i)
    for j = 1:num_trials 
    
    plot_periodogram(ROI_poststim_total_indivtr(j,i,:),fs,'semilog') 
    hold on 
    plot_periodogram(ROI_poststim_total_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_poststim_total_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_poststim_total_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_poststim_total_indivtr(j,i,:),fs,'semilog') 
    xlabel('Frequency (Hz)')
        ylabel('PSD (dB/Hz)')

    legend('Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Location','southwest'); 
    title(ROI_labels{i})
    end 
    
    name = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/poststim_semilog_indivtr.fig',PatientID,experiment_name); 
    saveas(gcf,name); 
end 

%% repeat for stim 
    figure('Position',[400,100,2000,800])

for i = 1:14
         subplot(2,7,i)

    for j = 1:num_trials 
        

    plot_periodogram(ROI_stim_indivtr(j,i,:),fs,'semilog') 
    hold on 
    plot_periodogram(ROI_stim_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_stim_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_stim_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_stim_indivtr(j,i,:),fs,'semilog') 
    legend('Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Location','southwest'); 
    title(ROI_labels{i})
      xlabel('Frequency (Hz)')
        ylabel('PSD (dB/Hz)')
    end 
    
    name = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/stim_semilog_indivtr.fig',PatientID,experiment_name); 
    saveas(gcf,name); 
end 



%% repeat for prestim 
    figure('Position',[400,100,2000,800])

for i = 1:14
         subplot(2,7,i)

    for j = 1:num_trials 

    plot_periodogram(ROI_pre_stim_indivtr(j,i,:),fs,'semilog') 
    
    hold on 
    plot_periodogram(ROI_pre_stim_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_pre_stim_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_pre_stim_indivtr(j,i,:),fs,'semilog') 
    plot_periodogram(ROI_pre_stim_indivtr(j,i,:),fs,'semilog') 
    legend('Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Location','southwest');
      xlabel('Frequency (Hz)')
        ylabel('PSD (dB/Hz)')
    title(ROI_labels{i})
    end 
    name = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/PSD_channels/%s/prestim_semilog_indivtr.fig',PatientID,experiment_name); 
    saveas(gcf,name); 
end

%% Save variables 

metadata = timeseries_data.metadata; 
outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s',PatientID, experiment_name); 

if ~exist(outputdir,'dir'), mkdir(outputdir), end

thisfile = sprintf('%s_time_series_ROI.mat', experiment_name); 
fulldestination = fullfile(outputdir, thisfile); 


filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s/%s_time_series_ROI.mat',PatientID,experiment_name, experiment_name); 
save(fulldestination, 'ROI_labels', 'ROI_poststim_total','ROI_pre_stim',...
    'ROI_stim','ROI_poststim_total_indivtr','ROI_stim_indivtr',...
    'ROI_pre_stim_indivtr', 'metadata', 'timeseriesfile') 




