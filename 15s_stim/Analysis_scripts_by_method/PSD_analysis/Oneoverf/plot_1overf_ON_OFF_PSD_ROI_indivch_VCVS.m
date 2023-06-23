%%%% The goal of this script is to generate PSD graphs showing differences
%%%% between ON vs OFF and generate 1/f trend lines for each stim
% %%%% condition b/w ON vs OFF states
% *** might have to epoch time domain data 
% A) do with no PARRM (ON vs OFF) 
% B) do with PARRM (ON vs OFF) 
% C) do with PARRM (ON vs ON with PARRM, OFF vs OFF with PARRM) 

% clear all 
% PatientID = 'DBSTRD001'; 
% DBS_target = 'SCC'; 
% hemi = 'l'; 
% stim_freq = 130; 
% experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq);  
% 
% % Load single trial time series data 
% timeseriesfile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
%     PatientID,hemi,DBS_target,stim_freq,hemi,DBS_target,stim_freq); 
% timeseries_data = load(timeseriesfile); 
% fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 
% ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 
% experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 

%% Run again 

clearvars -except PatientID DBS_target hemi stim_freq experiment_name timeseriesfile timeseries_data fs ch_labels experiment_name
contact_config = 'elec234'; 

%% Assign data for each condition 
disp(experiment_name)
switch DBS_target 
    case 'VCVS'
         [elec1, elec234, elec25, elec36, elec47,...
            elec567, elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi); 
    case 'SCC'
        [elec1, elec234, elec25, elec36, elec47,...
            elec567, elec8] = assign_SCC_conditions(timeseries_data,hemi); 
end 

%% Describe the data as a sanity check.
disp("The data consists of the power in *physical units* for each bin");
disp("Size of data (trials x channels x freq x time):");


switch contact_config 
    case 'elec25'
        data = elec25; 
    case 'elec1'
        data = elec1; 
    case 'elec8'
        data = elec8; 
    case 'elec36'
        data = elec36; 
    case 'elec47'
        data = elec47; 
    case 'elec234'
        data = elec234; 
    case 'elec567'
        data = elec567; 
end 

%% load ROI info 

load(sprintf('E:/DBSTRD/%s/Experiments/ROI_labels_%s.mat',PatientID,PatientID)); 
 %% Plot time-series traces for windows for inspection 
 
%confirm that there is no artifact during the time windows defined 
type = 'figure'; 
dir_name = 'Time series Plots'; 
%generate output directory 
outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config); 

disp(size(data));
%num_trials = size(data, 1);
num_trials = 5; 
num_ch = size(data, 2);
num_samples = size(data, 3);
t = (0:(num_samples-1))./fs; 

plot_time_series_channel(t, data , num_trials, ...
    num_ch, PatientID,contact_config, fs, ch_labels, experiment_name, outputdir) 

%% Define the different time windows of interest, in samples.

switch PatientID 
    case 'DBSTRD001' 
        pre_stim_win = 1:5000; 
        stim_win = 5001:20000; 
        post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
        post_stim_win2 = 25001:30000; 
        post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
    case 'DBSTRD002'     
        pre_stim_win = 1:4500; 
        stim_win = 5001:20000; 
        post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
        post_stim_win2 = 25001:30000; 
        post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
        %post_stim_total_win = 21000:26000; 
end 


%% Get different stim windows with data averaged across trials 

avg_data = squeeze(mean(data,1)); 
size(avg_data)

avg_pre_stim = avg_data(:,pre_stim_win); 
avg_poststim_total = avg_data(:,post_stim_total_win); 
avg_stim = avg_data(:,stim_win); 

poststim1_data = avg_data(:,post_stim_win1); 
poststim2_data = avg_data(:,post_stim_win2); 
 
 
 %% Get different stim windows with indiv trial data 
 
 indivtr_pre_stim_data = data(:,:,pre_stim_win); 
 indivtr_stim_data = data(:,:,stim_win); 
 indivtr_poststim_total_data = data(:,:,post_stim_total_win); 
 
 
 %% Generate periodogram using mtspectrumc for all freqs (multitaper chronix fxn) 
 
params.Fs = 1000;
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 1;
params.tapers = [4 7]; 
prestim_color = [0.6706    0.1882    0.1882]; 
stim_color = [0.6549    0.3765    0.7098]; 
poststim_color = [0.5020    0.5020    0.5020]; 
type = 'figure'; 
dir_name = 'Multitaper PSD with CI for channels'; 
%generate output directory 
outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config) ; 

parfor i=1:num_ch
    figure()   
    %prestim 
        [S,f,Serr] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(:,i,:),[3,1,2])),params);
        f = f(2:end);
        S = S(2:end);
        Serr = Serr(:,2:end);       
        fz = f; 
        mean_pow = 10*log(S)'; %convert to dB 
        Serr = 10*log(Serr); % convert to dB 
        h = plot(fz,mean_pow,'color',prestim_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz,Serr,'color',prestim_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz, fliplr(fz)], ...
        [Serr(2, :), fliplr(Serr(1, :))], ...
        prestim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
     % stim 
        [S_stim,f_stim,Serr_stim] = mtspectrumc(squeeze(permute(indivtr_stim_data(:,i,:),[3,1,2])),params);
        f_stim = f_stim(2:end); 
        S_stim = S_stim(2:end); 
        Serr_stim = Serr_stim(:,2:end);
        Serr_stim = 10*log(Serr_stim); % convert to dB 
        fz_stim = f_stim; 
        mean_pow_stim = 10*log(S_stim)'; 
        h3 = plot(fz_stim,mean_pow_stim,'color',stim_color); 
        h3.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h4 = plot(fz_stim,Serr_stim,'color',stim_color) ; 
        h4(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        
        h4(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        patch(...
        [fz_stim, fliplr(fz_stim)], ...
        [Serr_stim(2, :), fliplr(Serr_stim(1, :))], ...
        stim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);     
        % shade in confidence intervals  
        
     % post stim 
        [S_poststim,f_poststim,Serr_poststim] = mtspectrumc(squeeze(permute(indivtr_poststim_total_data(:,i,:),[3,1,2])),params);
        f_poststim = f_poststim(2:end); 
        S_poststim = S_poststim(2:end); 
        Serr_poststim = Serr_poststim(:,2:end); 
        fz_poststim = f_poststim; 
        mean_pow_poststim = 10*log(S_poststim)';
        Serr_poststim = 10*log(Serr_poststim); % convert to dB 
        h5 = plot(fz_poststim,mean_pow_poststim,'color',poststim_color); 
        h5.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
    disp('processing')
        h6 = plot(fz_poststim,Serr_poststim,'color',poststim_color); 
        h6(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
       
        h6(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        
        %shade in confidence intervals 
        patch(...
        [fz_poststim, fliplr(fz_poststim)], ...
        [Serr_poststim(2, :), fliplr(Serr_poststim(1, :))], ...
        poststim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);      
        
        title(sprintf('%s',ch_labels(i))); 
        xlabel('Frequency(Hz)'); 
        ylabel('PSD(db/Hz'); 
        
        legend('prestim','stim','poststim','Location','northeast'); 
        
       filename = sprintf('%s/%s.fig',outputdir,ch_labels(i));
       saveas(gcf, filename); 
       filename = sprintf('%s/%s.png',outputdir,ch_labels(i));
       saveas(gcf, filename); 
       
        
end 

% add slope to data 
%make confidence intervals the same color as the main averaged PSD 

%% Generate periodograms for smaller frequency band 

params.Fs = 1000;
params.fpass = [1 40];
params.err = [1 0.05];
params.trialave = 1;
params.tapers = [4 8]; 


prestim_color = [0.6706    0.1882    0.1882]; 
stim_color = [0.6549    0.3765    0.7098]; 
poststim_color = [0.5020    0.5020    0.5020]; 

dir_name = 'Multitaper PSD with CI for channels low freqs'; 
%generate output directory 
outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config); 

parfor i=1:num_ch
    figure()   
    %prestim 
        [S,f,Serr] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(:,i,:),[3,1,2])),params);
        f = f(2:end);
        S = S(2:end);
        Serr = Serr(:,2:end);       
        fz = f; 
        mean_pow = 10*log(S)'; %convert to dB 
        Serr = 10*log(Serr); % convert to dB 
        h = plot(fz,mean_pow,'color',prestim_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz,Serr,'color',prestim_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz, fliplr(fz)], ...
        [Serr(2, :), fliplr(Serr(1, :))], ...
        prestim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
     % stim 
        [S_stim,f_stim,Serr_stim] = mtspectrumc(squeeze(permute(indivtr_stim_data(:,i,:),[3,1,2])),params);
        f_stim = f_stim(2:end); 
        S_stim = S_stim(2:end); 
        Serr_stim = Serr_stim(:,2:end);
        Serr_stim = 10*log(Serr_stim); % convert to dB 
        fz_stim = f_stim; 
        mean_pow_stim = 10*log(S_stim)'; 
        h3 = plot(fz_stim,mean_pow_stim,'color',stim_color); 
        h3.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h4 = plot(fz_stim,Serr_stim,'color',stim_color) ; 
        h4(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        
        h4(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        patch(...
        [fz_stim, fliplr(fz_stim)], ...
        [Serr_stim(2, :), fliplr(Serr_stim(1, :))], ...
        stim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);     
        % shade in confidence intervals  
        
     % post stim 
        [S_poststim,f_poststim,Serr_poststim] = mtspectrumc(squeeze(permute(indivtr_poststim_total_data(:,i,:),[3,1,2])),params);
        f_poststim = f_poststim(2:end); 
        S_poststim = S_poststim(2:end); 
        Serr_poststim = Serr_poststim(:,2:end); 
        fz_poststim = f_poststim; 
        mean_pow_poststim = 10*log(S_poststim)';
        Serr_poststim = 10*log(Serr_poststim); % convert to dB 
        h5 = plot(fz_poststim,mean_pow_poststim,'color',poststim_color); 
        h5.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        h6 = plot(fz_poststim,Serr_poststim,'color',poststim_color); 
        h6(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
       
        h6(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        
        %shade in confidence intervals 
        patch(...
        [fz_poststim, fliplr(fz_poststim)], ...
        [Serr_poststim(2, :), fliplr(Serr_poststim(1, :))], ...
        poststim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);      
        
        title(sprintf('%s',ch_labels(i))); 
        xlabel('Frequency(Hz)'); 
        ylabel('PSD(db/Hz'); 
        
        legend('prestim','stim','poststim','Location','northeast'); 
        
       filename = sprintf('%s/%s.fig',outputdir,ch_labels(i));
       saveas(gcf, filename); 
       filename = sprintf('%s/%s.png',outputdir,ch_labels(i));
       saveas(gcf, filename); 
       
        
end 
% add slope to data 
% add params for mtspectrumc used for metadata 
%% plot PSD without confidence intervals for low freq 


dir_name = 'Multitaper PSD without CI with oneoverf for channels low freqs'; 
%generate output directory 
outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config); 

% initialize variables 
       f_all_ch = {};  
       S_all_ch=  {}; 
       S_upper_all_ch = {}; 
       S_lower_all_ch = {};  
       f_stim_all_ch= {}; 
       S_stim_all_ch = {}; 
       Serr_stim_upper_all_ch = {};  
       Serr_stim_lower_all_ch = {}; 
      
       f_poststim_all_ch = {}; 
       S_poststim_all_ch = {}; 
       Serr_poststim_upper_all_ch = {}; 
       Serr_poststim_lower_all_ch= {};


parfor i=1:num_ch
    figure()   
    %prestim 
        [S,f,Serr] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(:,i,:),[3,1,2])),params);
        f = f(2:end);
        S = S(2:end);
        fz = f; 
        mean_pow = 10*log(S)'; %convert to dB 
        h = plot(fz,mean_pow,'color',prestim_color,'linewidth',1.25); 
        hold on
        % 1/f line 
        [mdl,slope,r_squared] = get_1overf_slope(fz,S);
        j = plot(mdl); 
        delete(j(1));
     % stim 
        [S_stim,f_stim,Serr_stim] = mtspectrumc(squeeze(permute(indivtr_stim_data(:,i,:),[3,1,2])),params);
        f_stim = f_stim(2:end); 
        S_stim = S_stim(2:end); 

        fz_stim = f_stim; 
        mean_pow_stim = 10*log(S_stim)'; 
        h3 = plot(fz_stim,mean_pow_stim,'color',stim_color,'linewidth',1.25);
        % get 1/f line 
        [mdl_stim,slope_stim,r_squared_stim] = get_1overf_slope(fz_stim,S_stim);
        j1 = plot(mdl_stim);        
        delete(j1(1));
        
     % post stim 
        [S_poststim,f_poststim,Serr_poststim] = mtspectrumc(squeeze(permute(indivtr_poststim_total_data(:,i,:),[3,1,2])),params);
        f_poststim = f_poststim(2:end); 
        S_poststim = S_poststim(2:end); 
        fz_poststim = f_poststim; 
        mean_pow_poststim = 10*log(S_poststim)';
        h5 = plot(fz_poststim,mean_pow_poststim,'color',poststim_color,'linewidth',1.25); 
        [mdl_poststim,slope_poststim,r_squared_poststim] = get_1overf_slope(fz_poststim,S_poststim);
        j2 = plot(mdl_poststim); 
        delete(j2(1));
       
        
        title(sprintf('%s',ch_labels(i))); 
        xlabel('Frequency(Hz)'); 
        ylabel('PSD(db/Hz'); 
        
        legend('prestim','stim','poststim','Location','northeast'); 
        
        dim_1 = [0.129374127501163 0.130380953551758 0.702678551019303 0.0630952369244326];
        str = sprintf('slope pre = %s; slope stim = %s; slope post = %s',...
            num2str(slope), num2str(slope_stim), num2str(slope_poststim));
        annotation('textbox',dim_1,'String',str,'FitBoxToText','on');
        
       filename = sprintf('%s/%s.fig',outputdir,ch_labels(i));
       saveas(gcf, filename); 
       filename = sprintf('%s/%s.png',outputdir,ch_labels(i));
       saveas(gcf, filename); 
        
       f_all_ch{i} = f; 
       S_all_ch{i} =  mean_pow; 
       Serr_upper_all_ch{i} = Serr(1,:); 
       Serr_lower_all_ch{i} = Serr(2,:); 
       f_stim_all_ch{i} = f_stim; 
       S_stim_all_ch{i} = mean_pow_stim; 
       Serr_stim_upper_all_ch{i} = Serr_stim(1,:); 
       Serr_stim_lower_all_ch{i} = Serr_stim(2,:); 
      
       f_poststim_all_ch{i} = f_poststim; 
       S_poststim_all_ch{i} = mean_pow_poststim;
       Serr_poststim_upper_all_ch{i} = Serr_poststim(1,:); 
       Serr_poststim_lower_all_ch{i} = Serr_poststim(2,:); 
       
      
end 
%% 
type = 'data'; 
dir_name = 'multitaper_spectrum_all_ch'; 
outputdir = make_directory(PatientID, type, dir_name, experiment_name,...
    contact_config) ; 
metadata = timeseries_data.metadata; 
metadata.good_ch = good_ch_labels; 
metadata.ROI = ROI_labels; 

thisfile = sprintf('%s_time_series_ROI.mat', experiment_name); 
fulldestination = fullfile(outputdir, thisfile); 

filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s/%s_time_series_indivch.mat',PatientID,experiment_name, experiment_name); 
save(fulldestination, 'good_ch_labels', 'metadata', 'timeseriesfile',...
    'f_all_ch','S_all_ch','Serr_upper_all_ch','Serr_lower_all_ch',...
    'f_stim_all_ch','S_stim_all_ch','Serr_stim_upper_all_ch',...
    'Serr_stim_lower_all_ch','f_poststim_all_ch','S_poststim_all_ch',...
    'Serr_poststim_upper_all_ch','Serr_poststim_lower_all_ch'); 

%% Generate ROIs and make periodograms
% load ROI 
type = 'figure'; 

%extract ROIs for individual trial data 

[ROI_pre_stim_indivtr] = generate_ROI_indiv_tr(indivtr_pre_stim_data,PatientID); % gives ROI x current dir matrix 
[ROI_stim_indivtr] = generate_ROI_indiv_tr(indivtr_stim_data,PatientID); 
[ROI_poststim_total_indivtr] = generate_ROI_indiv_tr(indivtr_poststim_total_data,PatientID); 

%extract ROIs for averaged data 
[ROI_pre_stim] = generate_ROI(avg_pre_stim,PatientID); 
[ROI_stim] = generate_ROI(avg_stim,PatientID);
[ROI_poststim_total] = generate_ROI(avg_poststim_total,PatientID); 

%% Generate ROI periodograms 
num_ROI = size(ROI_pre_stim,1); 
disp(num_ROI)
dir_name = 'Multitaper PSD without CI for ROI'; 
%generate output directory 
outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config); 
figure('Position',[400,100,2000,800])
for i = 1:num_ROI 
    
    subplot(2,9,i)
        [S,f,Serr] = mtspectrumc(squeeze(permute(ROI_pre_stim_indivtr(:,i,:),[3,1,2])),params);
        f = f(2:end);
        S = S(2:end);
        fz = f; 
        mean_pow = 10*log(S)'; %convert to dB 
        h = plot(fz,mean_pow,'color',prestim_color,'linewidth',1.25); 
        hold on
        % 1/f line 
%         [mdl,slope,r_squared] = get_1overf_slope(fz,S);
%         j = plot(mdl); 
%         delete(j(1));
     % stim 
        [S_stim,f_stim,Serr_stim] = mtspectrumc(squeeze(permute(ROI_stim_indivtr(:,i,:),[3,1,2])),params);
        f_stim = f_stim(2:end); 
        S_stim = S_stim(2:end); 

        fz_stim = f_stim; 
        mean_pow_stim = 10*log(S_stim)'; 
        h3 = plot(fz_stim,mean_pow_stim,'color',stim_color,'linewidth',1.25);
        % get 1/f line 
%         [mdl_stim,slope_stim,r_squared_stim] = get_1overf_slope(fz_stim,S_stim);
%         j1 = plot(mdl_stim);        
%         delete(j1(1));
        
     % post stim 
        [S_poststim,f_poststim,Serr_poststim] = mtspectrumc(squeeze(permute(ROI_poststim_total_indivtr(:,i,:),[3,1,2])),params);
        f_poststim = f_poststim(2:end); 
        S_poststim = S_poststim(2:end); 
        fz_poststim = f_poststim; 
        mean_pow_poststim = 10*log(S_poststim)';
        h5 = plot(fz_poststim,mean_pow_poststim,'color',poststim_color,'linewidth',1.25); 
%         [mdl_poststim,slope_poststim,r_squared_poststim] = get_1overf_slope(fz_poststim,S_poststim);
%         j2 = plot(mdl_poststim); 
%         delete(j2(1));
       
        
        title(sprintf('%s',ROI_labels{i})); 
        xlabel('Frequency(Hz)'); 
        ylabel('PSD(db/Hz'); 
        
        legend('prestim','stim','poststim','Location','southwest'); 
        
       
       filename = sprintf('%s/%s.fig',outputdir,ROI_labels{i});
       saveas(gcf, filename); 
       filename = sprintf('%s/%s.png',outputdir,ROI_labels{i});
       saveas(gcf, filename); 

end 



%% Generate ROI periodograms with CI 
dir_name = 'Multitaper PSD with CI for ROI'; 

%generate output directory 
outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config); 
figure('Position',[400,100,2200,800])


for i = 1:num_ROI
      
        subplot(2,9,i)
    %prestim 
        [S,f,Serr] = mtspectrumc(squeeze(permute(ROI_pre_stim_indivtr(:,i,:),[3,1,2])),params);
        f = f(2:end);
        S = S(2:end);
        Serr = Serr(:,2:end);       
        fz = f; 
        mean_pow = 10*log(S)'; %convert to dB 
        Serr = 10*log(Serr); % convert to dB 
        h = plot(fz,mean_pow,'color',prestim_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz,Serr,'color',prestim_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz, fliplr(fz)], ...
        [Serr(2, :), fliplr(Serr(1, :))], ...
        prestim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
     % stim 
        [S_stim,f_stim,Serr_stim] = mtspectrumc(squeeze(permute(ROI_stim_indivtr(:,i,:),[3,1,2])),params);
        f_stim = f_stim(2:end); 
        S_stim = S_stim(2:end); 
        Serr_stim = Serr_stim(:,2:end);
        Serr_stim = 10*log(Serr_stim); % convert to dB 
        fz_stim = f_stim; 
        mean_pow_stim = 10*log(S_stim)'; 
        h3 = plot(fz_stim,mean_pow_stim,'color',stim_color); 
        h3.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h4 = plot(fz_stim,Serr_stim,'color',stim_color) ; 
        h4(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        
        h4(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        patch(...
        [fz_stim, fliplr(fz_stim)], ...
        [Serr_stim(2, :), fliplr(Serr_stim(1, :))], ...
        stim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);     
        % shade in confidence intervals  
        
     % post stim 
        [S_poststim,f_poststim,Serr_poststim] = mtspectrumc(squeeze(permute(ROI_poststim_total_indivtr(:,i,:),[3,1,2])),params);
        f_poststim = f_poststim(2:end); 
        S_poststim = S_poststim(2:end); 
        Serr_poststim = Serr_poststim(:,2:end); 
        fz_poststim = f_poststim; 
        mean_pow_poststim = 10*log(S_poststim)';
        Serr_poststim = 10*log(Serr_poststim); % convert to dB 
        h5 = plot(fz_poststim,mean_pow_poststim,'color',poststim_color); 
        h5.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        h6 = plot(fz_poststim,Serr_poststim,'color',poststim_color); 
        h6(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
       
        h6(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        
        %shade in confidence intervals 
        patch(...
        [fz_poststim, fliplr(fz_poststim)], ...
        [Serr_poststim(2, :), fliplr(Serr_poststim(1, :))], ...
        poststim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);      
        
        title(sprintf('%s',ROI_labels{i})); 
        xlabel('Frequency(Hz)'); 
        ylabel('PSD(db/Hz'); 
        
        legend('prestim','stim','poststim','Location','southwest'); 
        
       filename = sprintf('%s/%s.fig',outputdir,ROI_labels{i});
       saveas(gcf, filename); 
       filename = sprintf('%s/%s.png',outputdir,ROI_labels{i});
       saveas(gcf, filename); 
       
        
end 

%% Save spectral output over pre-stim, stim and post-stim periods from multitaper






%% Get spectral output for 1-s bins 








%% save data 
% type = 'data'; 
% outputdir = make_directory(PatientID, type, dir_name, experiment_name,...
%     contact_config) ; 
% metadata = timeseries_data.metadata; 
% metadata.good_ch = good_ch_labels; 
% metadata.ROI = ROI_labels; 
% 
% thisfile = sprintf('%s_time_series_ROI.mat', experiment_name); 
% fulldestination = fullfile(outputdir, thisfile); 
% 
% filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/%s/%s_time_series_ROI.mat',PatientID,experiment_name, experiment_name); 
% save(fulldestination, 'good_ch_labels', 'metadata', 'timeseriesfile', ...
%     'f', 'S', 'f_stim','S_stim','f_poststim','S_poststim',...
%     'Serr','Serr_stim','Serr_poststim') 






