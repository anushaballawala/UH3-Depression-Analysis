%% 
clear all 
PatientID = 'DBSTRD001'; 
experiment = 'BaselineFix'; 
experiment_name = 'BaselineFix_run-01_blk-01'; 

name_of_file = '12-Oct-2021_BaselineFix_run-07_blk-01_mtspectrum_pow_indivtr.mat'; 
ROI_labels = load(sprintf('ROI_labels_%s.mat',PatientID)); 
num_ROI = length(ROI_labels.ROI_labels); 
load(name_of_file)
%% Generate ROI periodograms with CI 

%generate output directory 
% outputdir = make_directory(PatientID, type, dir_name, experiment_name, contact_config); 

%prestim_color = uisetcolor(); 
stim1_color = uisetcolor();
stim2_color = uisetcolor();
stim3_color = uisetcolor();
poststim_color = uisetcolor();

%% Trial averaged 

S_poststim = squeeze(mean(S_poststim,1)); 
S_stim1 = squeeze(mean(S_stim1,1));  
S_stim2= squeeze(mean(S_stim2,1)); 
S_stim3= squeeze(mean(S_stim3,1)); 
S_prestim = squeeze(mean(S_prestim,1)); 

Serr_poststim = squeeze(mean(Serr_poststim,1)); 
Serr_stim1 = squeeze(mean(Serr_stim1,1));
Serr_stim2 = squeeze(mean(Serr_stim2,1));
Serr_stim3 = squeeze(mean(Serr_stim3,1));
Serr_prestim = squeeze(mean(Serr_prestim,1));

%% Get ROI 

prestim_ROI = generate_ROI_tmp(S_prestim,PatientID); 
stim1_ROI = generate_ROI_tmp(S_stim1,PatientID); 
stim2_ROI = generate_ROI_tmp(S_stim2,PatientID); 
stim3_ROI = generate_ROI_tmp(S_stim3,PatientID); 
poststim_ROI = generate_ROI_tmp(S_poststim,PatientID);
%% 
Serr_prestim_ROI = generate_ROI_Serr(Serr_prestim,PatientID);
Serr_stim1_ROI = generate_ROI_Serr(Serr_stim1,PatientID);
Serr_stim2_ROI = generate_ROI_Serr(Serr_stim2,PatientID);
Serr_stim3_ROI = generate_ROI_Serr(Serr_stim3,PatientID);
Serr_poststim_ROI = generate_ROI_Serr(Serr_poststim,PatientID);


%% 
figure('Position',[400,100,2200,800])

for i = 1:num_ROI
        subplot(2,9,i)
        %%%%%%%%%prestim 
%         fz = f_prestim(2:408);
%         S_prestim_ROI = 10*log(prestim_ROI(i,2:408)); 
%         Serr_prestim_ROI_plot = squeeze(10*log(Serr_prestim_ROI(i,:,2:408)));
%         
%         
%         h = plot(fz,S_prestim_ROI,'color',prestim_color); 
%         h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
%         hold on
%         h1 = plot(fz,Serr_prestim_ROI_plot,'color',prestim_color); 
%         h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
%         h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
% 
%         % shade confidence intervals 
%         patch(...
%         [fz, fliplr(fz)], ...
%         [Serr_prestim_ROI_plot(2, :), fliplr(Serr_prestim_ROI_plot(1, :))], ...
%         prestim_color, ...
%         'FaceAlpha', 0.3, 'EdgeAlpha',0);
    
    %%%%%%%%%%%%%stim 
        fz_stim1 = f_stim1(2:408);
        S_stim1_ROI = 10*log(stim1_ROI(i,2:408)); 
        Serr_stim1_ROI_plot = squeeze(10*log(Serr_stim1_ROI(i,:,2:408)));
        
        
        h = plot(fz_stim1,S_stim1_ROI,'color',stim1_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz_stim1,Serr_stim1_ROI_plot,'color',stim1_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz_stim1, fliplr(fz_stim1)], ...
        [Serr_stim1_ROI_plot(2, :), fliplr(Serr_stim1_ROI_plot(1, :))], ...
        stim1_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
    
    
    %%%%%%%%%% stim 2
    
        fz_stim1 = f_stim1(2:408);
        S_stim2_ROI = 10*log(stim2_ROI(i,2:408)); 
        Serr_stim2_ROI_plot = squeeze(10*log(Serr_stim2_ROI(i,:,2:408)));
        
        
        h = plot(fz_stim1,S_stim2_ROI,'color',stim2_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz_stim1,Serr_stim2_ROI_plot,'color',stim2_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz_stim1, fliplr(fz_stim1)], ...
        [Serr_stim2_ROI_plot(2, :), fliplr(Serr_stim2_ROI_plot(1, :))], ...
        stim2_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
    
    %%%%%%%%% stim3
    
        fz_stim1 = f_stim1(2:408);
        S_stim3_ROI = 10*log(stim3_ROI(i,2:408)); 
        Serr_stim3_ROI_plot = squeeze(10*log(Serr_stim3_ROI(i,:,2:408)));
        
        
        h = plot(fz_stim1,S_stim3_ROI,'color',stim3_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz_stim1,Serr_stim3_ROI_plot,'color',stim3_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz_stim1, fliplr(fz_stim1)], ...
        [Serr_stim3_ROI_plot(2, :), fliplr(Serr_stim3_ROI_plot(1, :))], ...
        stim3_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
    
    %%%%%%%%% post-stim 
    
        fz_poststim = f_poststim(2:408);
        S_poststim_ROI = 10*log(poststim_ROI(i,2:408)); 
        Serr_poststim_ROI_plot = squeeze(10*log(Serr_poststim_ROI(i,:,2:408)));
        
        
        h = plot(fz_poststim,S_poststim_ROI,'color',poststim_color); 
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        hold on
        h1 = plot(fz_poststim,Serr_poststim_ROI_plot,'color',poststim_color); 
        h1(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        h1(2,1).Annotation.LegendInformation.IconDisplayStyle = 'off'; 

        % shade confidence intervals 
        patch(...
        [fz_poststim, fliplr(fz_poststim)], ...
        [Serr_poststim_ROI_plot(2, :), fliplr(Serr_poststim_ROI_plot(1, :))], ...
        poststim_color, ...
        'FaceAlpha', 0.3, 'EdgeAlpha',0);
  
    title(ROI_labels.ROI_labels{i}); 
    legend ('stim1','stim2','stim3','poststim')

    
end 