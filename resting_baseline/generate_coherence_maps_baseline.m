
function t = generate_coherence_maps_baseline(data_for_fig,FOI,PatientID,...
    experiment_name,caxis_limits,trial_type)
f = figure; 
f.Position = [100 100 800 600]; 

        imagesc(data_for_fig) 
        title(sprintf('%s',FOI))
        colormap redbluecmap ; 
        colorbar 
        caxis(caxis_limits); 
        t.TileSpacing = 'compact'; 
        t.Padding = 'compact'; 

        outputdir = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/Coherence_Maps/%s/%s',PatientID,experiment_name,date); 
            if ~exist(outputdir,'dir')
                mkdir(outputdir)
            end 
                
        filename = sprintf('%s/%s_%s_%s',outputdir,experiment_name,FOI,trial_type); 
        filename_png = sprintf('%s/%s_%s_%s.png',outputdir,experiment_name,FOI,trial_type); 
        saveas(gcf,filename); 
        saveas(gcf,filename_png); 
     
end 