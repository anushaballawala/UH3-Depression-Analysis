%% extract ROIs
load('E:/DBSTRD/DBSTRD001/Experiments/ROI_labels_DBSTRD001.mat'); 
ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG','lVPF',...
    'rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rVPF'}; 
% This gives trials x ROI x time 
[ROI_rSCC_stim] = generate_ROI_indiv_tr(rSCC_stim_e8);
[ROI_rSCC_poststim1] = generate_ROI_indiv_tr(rSCC_poststim1_e8); 
[ROI_rSCC_poststim2] = generate_ROI_indiv_tr(rSCC_poststim2_e8); 

% Average across time. 

ROI_rSCC_stim_avgtime = mean(ROI_rSCC_stim,3); 
ROI_rSCC_poststim1_avgtime = mean(ROI_rSCC_poststim1,3); 
ROI_rSCC_poststim2_avgtime = mean(ROI_rSCC_poststim2,3); 

%% 

x1 = 1*ones(1,length(ROI_rSCC_stim_avgtime(:,1))); 
x2 = 2*ones(1,length(ROI_rSCC_stim_avgtime(:,1))); 
x3 = 3*ones(1,length(ROI_rSCC_stim_avgtime(:,1))); 

%% Plot individual trials for each ROI 
g = [0.5020 0.5020 0.5020]; 
num_ROI = length(ROI_labels); 
f = figure; 
f.Position = [1,250,1500,1000]; 
prompt = 'enter electroe config'; 
econfig = inputdlg(prompt,'str'); 
econfig = econfig{1}; 
for i = 1:num_ROI
    subplot(2,7,i)
    CT=cbrewer('seq', 'YlGnBu', 3);
    scatter(x1,ROI_rSCC_stim_avgtime(:,i),40,CT(1,:),'filled','MarkerEdgeColor',g,'LineWidth',0.75)
    hold on
    scatter(x2,ROI_rSCC_poststim1_avgtime(:,i),40,CT(2,:),'filled','MarkerEdgeColor',g,'LineWidth',0.75)
    scatter(x3,ROI_rSCC_poststim2_avgtime(:,i),40,CT(3,:),'filled','MarkerEdgeColor',g,'LineWidth',0.75)
    xticks([1 2 3])
    xticklabels({'stim ON', 'OFF1', 'OFF2'})
    ylabel('% change from BL')
    title(sprintf('%s-%s',ROI_labels{i},freqband))
    filename = (sprintf('%s_%s_%s.png',DBStarget,freqband,econfig)); 
    
end
saveas(gcf, filename)
close all 
