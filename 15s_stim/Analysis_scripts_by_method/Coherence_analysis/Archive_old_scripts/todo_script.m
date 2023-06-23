% 1. save pre-stim data from first round of stimulation (check log) 
% 2. load saved matrices with coherence data for prestim, stim, poststim 
% 3. subtract stim from original pre, and post from original pre
% 4. generate coherence matrices 
% 5. get ROI for coherence matrices 
% 6. perform preprocessing on baselien data before any stim was performed
% (baseline recording). 
%. perform preprocessing on baseline data before any stim was performed 

%% load 
clear all 
PatientID = 'DBSTRD001'; 
experiment_name = 'lVCVS_f130'; 
filename_prestim = dir(sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/%s/prestim*',...
    PatientID, experiment_name)); 
filename_stim = dir(sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/%s/stim*',...
    PatientID, experiment_name)); 
filename_poststim = dir(sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/Coherence/%s/poststim*',...
    PatientID, experiment_name)); 


prestim = load(filename_prestim.name); 
stim = load(filename_stim.name); 
poststim = load(filename_poststim.name); 
names = poststim.metadata.coherence.contact_config_order ; 
%% Create difference coherence matrix  

for i = 1:7                                                                                                                                      
alpha_stimminuspre{1,i} = stim.matrix_coh_alpha{1,i} - prestim.matrix_coh_alpha{1,i}; 
alpha_postminuspre{1,i} = poststim.matrix_coh_alpha{1,i} - prestim.matrix_coh_alpha{1,i}; 

beta_stimminuspre{1,i} = stim.matrix_coh_beta{1,i} - prestim.matrix_coh_beta{1,i}; 
beta_postminuspre{1,i} = poststim.matrix_coh_beta{1,i} - prestim.matrix_coh_beta{1,i};

theta_stimminuspre{1,i} = stim.matrix_coh_theta{1,i} - prestim.matrix_coh_theta{1,i}; 
theta_postminuspre{1,i} = poststim.matrix_coh_theta{1,i} - prestim.matrix_coh_theta{1,i};

delta_stimminuspre{1,i} = stim.matrix_coh_delta{1,i} - prestim.matrix_coh_delta{1,i}; 
delta_postminuspre{1,i} = poststim.matrix_coh_delta{1,i} - prestim.matrix_coh_delta{1,i};

end 


%% generate difference plots 

stim_state = 'stimminuspre'; 
t_delta = generate_coh_maps(delta_stimminuspre,names,'delta',PatientID,experiment_name,stim_state); 
t_theta = generate_coh_maps(theta_stimminuspre,names,'theta',PatientID,experiment_name,stim_state); 
t_alpha = generate_coh_maps(alpha_stimminuspre,names,'alpha',PatientID,experiment_name,stim_state); 
t_beta  = generate_coh_maps(beta_stimminuspre,names,'beta',PatientID,experiment_name,stim_state); 

stim_state_post = 'postminuspre'; 
t_delta_post = generate_coh_maps_post(delta_postminuspre,names,'delta',PatientID,experiment_name,stim_state_post); 
t_theta_post = generate_coh_maps_post(theta_postminuspre,names,'theta',PatientID,experiment_name,stim_state_post); 
t_alpha_post = generate_coh_maps_post(alpha_postminuspre,names,'alpha',PatientID,experiment_name,stim_state_post); 
t_beta_post  = generate_coh_maps_post(beta_postminuspre,names,'beta',PatientID,experiment_name,stim_state_post); 


% function 
function t = generate_coh_maps(data_for_fig,names,FOI,PatientID,experiment_name,stim_state)
f = figure; 
f.Position = [100 100 1400 800]; 
t = tiledlayout(2,4); 
for i = 1:7 
nexttile 
imagesc(data_for_fig{i}) 
title(sprintf('%s %s',names{i},FOI))
colormap redbluecmap ; 
colorbar 
caxis([0 0.3]); 
t.TileSpacing = 'compact'; 
t.Padding = 'compact'; 

% filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081621_Coherence/%s_%s_%s',...
%     PatientID,stim_state,DBS_target,FOI); 
filename_png = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081721_CoherenceDiffPlots/%s_%s_%s.png',...
    PatientID,stim_state,experiment_name,FOI); 
% saveas(gcf,filename); 
saveas(gcf,filename_png); 
end 
end 
 

% function 
function t = generate_coh_maps_post(data_for_fig,names,FOI,PatientID,experiment_name,stim_state_post)
f = figure; 
f.Position = [100 100 1400 800]; 
t = tiledlayout(2,4); 
for i = 1:7 
nexttile 
imagesc(data_for_fig{i}) 
title(sprintf('%s %s',names{i},FOI))
colormap redbluecmap ; 
colorbar 
caxis([0 0.3]); 
t.TileSpacing = 'compact'; 
t.Padding = 'compact'; 

% filename = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081621_Coherence/%s_%s_%s',...
%     PatientID,stim_state,DBS_target,FOI); 
filename_png = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Figures/AA/081721_CoherenceDiffPlots/%s_%s_%s.png',...
    PatientID,stim_state_post,experiment_name,FOI); 
% saveas(gcf,filename); 
saveas(gcf,filename_png); 
end 
end 



%% create ROI plot 

%load ROI labels 


%average coherence numbers across ROI for pre-stim and stim 
%use generate ROI_average trial function 
%% 

names = metadata.coherence.contact_config_order; 
PatientID = 'DBSTRD001'; 
num_contact_configs = 7; 
experiment_name = 'lSCC_f130'; 
stim_state = 'poststim'; 
caxis_limits = [0 0.5]; 

t_theta = generate_coherence_maps(matrix_coh_theta,names,...
    'theta',PatientID,num_contact_configs,experiment_name,stim_state,caxis_limits);






