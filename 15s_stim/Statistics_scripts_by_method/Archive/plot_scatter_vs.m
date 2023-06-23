function [] =  plot_scatter_vs(PatientID,metadata,p_idx,freqband_idx,cond1,cond2,freqband,data_cond1_mean,data_cond2_mean,data_dim)

        figure()
        colormap = [0.6353    0.0784    0.1843; 0    0.4471    0.7412];  %red, blue -> above, below
        
        A = data_cond1_mean; % A and B are examples of one of the n plot data
        B = data_cond2_mean;
                idx = 1 + (A < B);
                color = colormap(idx,:);
                pointsize = 60;
                s = scatter (A, B, pointsize, color, 'o','filled', 'MarkerEdgeColor',[0.5 .5 .5]);
                s.MarkerFaceAlpha = 0.6;
                xlabel(cond1)
                ylabel(cond2)
                
                min(A)
                
                hold on
%         switch data_dim
%             case 'Ch'
                Diag_X = [min(A),max(A)]; % Diag_X and Diag_Y form a reference line that is common in all the plot data
                Diag_Y = [min(A),max(A)];
                
                plot (Diag_X, Diag_Y, 'k');
                
                title(sprintf('%s power in stim in sig. channels',freqband));
                
                %save
                cd(sprintf('/Users/anushaallawala/Data/stats_results/two-tail_pos_neg_Tstat/figs/%s', data_dim))
                
                figname = sprintf('%s_%s_%s_%s_%s_%s-elec_%s.png',PatientID,metadata.comparison,...
                    metadata.cond1,freqband,metadata.stim_cond,metadata.elec_cond_type,metadata.elec_config);
                
                saveas(gcf,figname);
                
%             case 'ROI'
                
                title(sprintf('%s power in stim in sig. ROIs',freqband));
                
                cd(sprintf('/Users/anushaallawala/Data/stats_results/two-tail_pos_neg_Tstat/figs/%s', data_dim))
                
                figname = sprintf('%s_%s_%s_%s_%s_%s-elec_%s.png',PatientID,metadata.comparison,...
                    metadata.cond1,freqband,metadata.stim_cond,metadata.elec_cond_type,metadata.elec_config);
                
                saveas(gcf,figname);
                
%         end
end
