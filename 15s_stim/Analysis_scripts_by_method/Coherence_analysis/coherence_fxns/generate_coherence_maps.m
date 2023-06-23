
function t = generate_coherence_maps(data_for_fig,FOI,PatientID,...
    num_contact_configs,experiment_name,stim_state,caxis_limits)
f = figure; 
f.Position = [100 100 1400 800]; 
t = tiledlayout(2,4); 
    for i = 1:num_contact_configs 
        nexttile 
        imagesc(data_for_fig{i}) 
        title(sprintf('%s %s',FOI))
        colormap redbluecmap ; 
        colorbar 
        caxis(caxis_limits); 
        t.TileSpacing = 'compact'; 
        t.Padding = 'compact'; 

        outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Figures/AA/Coherence_Maps/%s/%s',PatientID,experiment_name,date); 
            if ~exist(outputdir,'dir')
                mkdir(outputdir)
            end 
                
        filename = sprintf('%s/%s_%s_%s',outputdir,stim_state,experiment_name,FOI); 
%         filename_png = sprintf('%s/%s_%s_%s.png',outputdir,stim_state,experiment_name,FOI); 
%         saveas(gcf,filename); 
%         saveas(gcf,filename_png); 
    end 
end 