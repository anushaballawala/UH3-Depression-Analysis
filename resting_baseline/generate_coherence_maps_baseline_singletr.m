


function t = generate_coherence_maps_baseline_singletr(data_for_fig,FOI,PatientID,...
    experiment_name,caxis_limits,trial_type,num_trials)

f = figure; 
f.Position = [100 100 1400 800]; 
t = tiledlayout(3,3); %change to num trials 
    for i = 1:num_trials
        nexttile 
        imagesc(data_for_fig{i}) 
        title(sprintf('Trial Num%d %s',i,FOI))
        colormap redbluecmap ; 
        colorbar 
        caxis(caxis_limits); 
        t.TileSpacing = 'compact'; 
        t.Padding = 'compact'; 

    end 

        outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/Coherence_Maps/%s_%s/%s',PatientID,experiment_name,trial_type,date); 
            if ~exist(outputdir,'dir')
                mkdir(outputdir)
            end 
                
        filename = sprintf('%s/%s_%s_%s',outputdir,experiment_name,FOI,trial_type); 
        filename_png = sprintf('%s/%s_%s_%s.png',outputdir,experiment_name,FOI,trial_type); 
        saveas(gcf,filename); 
        saveas(gcf,filename_png); 
     
end 