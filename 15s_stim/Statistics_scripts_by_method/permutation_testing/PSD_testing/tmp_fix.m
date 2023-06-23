%%%%%%%%%Find channels that are significantly different. 
clear baseline_plot_data stim_plot_data

%format p value 
alpha = 0.05; 
metadata.alpha = alpha; 

% Get p-values less than 0.05 
flat_p_idx = find(p_value < alpha);
[Ch,freq,conds,dbstarg] = ind2sub(size(p_value), flat_p_idx) ;
p_idx = [Ch,freq,conds,dbstarg];

%------------INFO------------% 
% size
%     10     6    24     2
%     Ch X freqband X conds X DBS targ

% dbstargvalue --> SCC = 1, VCVS = 2
% condvalue -----> 1:SCC_idx (24), VCVS_idx(25):end. 
% freqvalue -----> delta = 1, theta = 2, alpha = 3, beta = 4, low gamma =
% 5, high gamma = 6. 

%---------END OF INFO--------% 

dbstargvalue = 2; 
condvalue = 1; 
freqvalue = 6; 

% Get the indices for significant ROIs/Channels 
logical_condition = (freq == freqvalue & conds == condvalue & dbstarg == dbstargvalue);
corresponding_channels = Ch(logical_condition);
corresponding_freqs = freq(logical_condition); % But this is all 1
flat_p_idx_for_any_channel_but_fcd_is_111 = flat_p_idx(logical_condition);
significant_p_vals_for_any_channel_but_111 = p_value(flat_p_idx_for_any_channel_but_fcd_is_111);

%pseudocode *** 
if dbstargvalue == 1 
    plot_data = mean_data_SCC;  
    plot_labels = comp_labels_SCC; 
elseif dbstargvalue == 2
    plot_data = mean_data_VCVS; 
    plot_labels = comp_labels_VCVS;
end 

baseline_plot_data = mean_data_baseline(corresponding_channels,freqvalue); 

stim_plot_data = plot_data{condvalue}(corresponding_channels,freqvalue); 

% can change logical condition based on plot we want to generate for freq,
% condition and DBS target - we have the legend for what each entry is
% because we made that earlier 

%Metadata ** 

cond1 = plot_labels{1, condvalue}.DBS  ; 
cond2 = 'baseline'; 
elec = plot_labels{1, condvalue}.elec  ; 

metadata.cond1 = plot_labels{1, condvalue}.DBS  ; 
metadata.cond2 = 'baseline'; 

metadata.stim_cond = plot_labels{1, condvalue}.Stim_state ; 
metadata.elec_cond_type = plot_labels{1, condvalue}.elec  ;

                
%%%%%%% select which figures to plot 

%list_freqs = {'delta','theta','alpha','beta','lowgamma','highgamma'}; 

if freqvalue == 1
    freq = 'delta'
elseif freqvalue == 2
    freq = 'theta'
elseif freqvalue == 3
    freq = 'alpha'
elseif freqvalue == 4
    freq = 'beta'
elseif freqvalue == 5 
    freq = 'lowgamma'
elseif freqvalue == 6
    freq = 'highgamma'
end 

metadata.freqband = freq; 

 figure()
        colormap = [0.6353    0.0784    0.1843; 0    0.4471    0.7412];  %red, blue -> above, below
        
        % scatter(allstimdata, allbaseline)
        
        % get idx for channels significant 
        
        % circle channels that are significant 
        
        %plot all channels/ROIs. 
                        
                A = stim_plot_data; % A and B are examples of one of the n plot data
                B = baseline_plot_data;

                idx = 1 + (A < B);
                color = colormap(idx,:);
                pointsize = 110;
                s = scatter (A, B, pointsize, color, 'o','filled', 'MarkerEdgeColor',[0.5 .5 .5]);
                s.MarkerFaceAlpha = 0.6;
                s.LineWidth = 2;
                xlabel(sprintf('%s poststim',cond1))
                ylabel(cond2)
                
                hold on 
                
                A_all = plot_data{condvalue}(:,freqvalue); 
                B_all = mean_data_baseline(:,freqvalue); 
                
                idx_all = 1 + (A_all < B_all);
                color_all = colormap(idx_all,:);
                
                s1  = scatter(A_all,B_all,85,color_all,'o','filled');  
                s1.MarkerFaceAlpha = 0.6; 
                

                dx = 0.005; dy = 0.005; % displacement so the text does not overlay the data points
                text(A_all+dx, B_all+dy, ROI_labels,'FontSize',12,'FontName','Avenir');
                min(A)
                
                Diag_X = [min(A_all),max(A_all)]; % Diag_X and Diag_Y form a reference line that is common in all the plot data
                Diag_Y = [min(A_all),max(A_all)];
                
                plot (Diag_X, Diag_Y, 'k');
                
                title(sprintf('%s power in stim in sig. ROI %s',freq,metadata.elec_cond_type));
                addpath('/Users/anushaallawala/Data/stats_results/')
                cd(sprintf('/Users/anushaallawala/Data/stats_results/two-tail_pos_neg_Tstat_zscore/figs/%s/bl_poststim', data_dim))
                
                figname = sprintf('new_%s_%s_%s_%s_elec_%s_elec-%s.png',PatientID,metadata.cond2,...
                    metadata.cond1,freq,metadata.stim_cond{1},metadata.elec_cond_type);
                
                %saveas(gcf,figname);
