%% oneoverf_periodogram.m - Gives 1/f slope for each stim state 

% Additional info: *** pwelch params need to be changed and code has to be 
% cleaned up
% Inputs: N/A
% Outputs: plots of 1/f overlaid with PSD 
% Dependencies: UH3 github repo 
% Sub-functions: none 

% Anusha Allawala, 10/2021

%------------ START OF CODE ---------------% 

%% Generate periodogram for ON vs OFF for all channels 
parfor i = 1:num_ch
    
    figure('Position',[500 200 1500 1000])
    subplot(2,1,1)
    [Pxx,f,mdl,slope] = plot_1overf(avg_pre_stim(i,:),fs,'180'); 
    hold on 
    [Pxx_stim,f_stim,mdl_stim,slope_stim] = plot_1overf(avg_stim(i,:),fs,'180') 
    [Pxx_poststim,f_poststim,mdl_poststim,slope_poststim] = plot_1overf(avg_poststim_total(i,:),fs,'180'); 
    title(ch_labels(i))
    legend('prestim','prestim','prestim','stim','stim','stim','poststim','poststim','poststim','Location','southwest')
    dim_1 = [0.122231270358306 0.926500000491738 0.332247547344199 0.0264999995082617];
    rsqr = mdl.Rsquared.Ordinary; 
    slope = mdl.Coefficients.Estimate(2); 
    str = sprintf('slope pre = %s; slope stim = %s; slope post = %s',...
        num2str(slope), num2str(slope_stim), num2str(slope_poststim)); 
    annotation('textbox',dim_1,'String',str,'FitBoxToText','on'); 
    

    subplot(2,1,2)
    [Pxx,f,mdl,slope] = plot_1overf(avg_pre_stim(i,:),fs,'50'); 
    hold on 
    [Pxx_stim,f_stim,mdl_stim,slope_stim] = plot_1overf(avg_stim(i,:),fs,'50');
    [Pxx_poststim,f_poststim,mdl_poststim,slope_poststim] = plot_1overf(avg_poststim_total(i,:),fs,'50'); 
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