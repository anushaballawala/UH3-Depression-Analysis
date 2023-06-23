%% Load baseline data (prestim windows) 
PatientID = 'DBSTRD001'; 
experiment_name = 'rVCVS_f130'; 
stim_state = 'poststim'; 
%% Load data 
baseline_data = load(sprintf('BaselineFix_run-07_blk-01_avgtrials_fake_%s_matlabcohfxn.mat',...
    stim_state)); 
stim_data = load(sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed data/Coherence/%s/%s_%s_matlabcohfxn.mat',...
    PatientID,experiment_name,stim_state,experiment_name)); 

num_contact_configs = length(stim_data.matrix_coh_alpha); 
bl_alpha_coh = baseline_data.avg_alpha; 
bl_delta_coh = baseline_data.avg_delta; 
bl_theta_coh = baseline_data.avg_theta;

%% subtract average of coherence matrices to generate difference plot 

coh_diff_alpha = cell(1,num_contact_configs); 
coh_diff_theta = cell(1,num_contact_configs); 
coh_diff_delta = cell(1,num_contact_configs); 
for i = 1:num_contact_configs 
    
    coh_diff_alpha{i} = stim_data.matrix_coh_alpha{1, i} -  bl_alpha_coh; 
    coh_diff_theta{i} = stim_data.matrix_coh_theta{1, i} - bl_theta_coh; 
    coh_diff_delta{i} = stim_data.matrix_coh_delta{1, i} - bl_delta_coh;
    
end 

%% make coherence plots and save 

caxis_limits = [-0.2 0.1]; 
names = stim_data.metadata.coherence.contact_config_order ; 
fig_name = sprintf('%s_diffmaps',experiment_name); 

t_delta = generate_coherence_maps(coh_diff_delta,names,'delta',PatientID,...
    num_contact_configs,fig_name,stim_state,caxis_limits); 

t_theta = generate_coherence_maps(coh_diff_theta,names,'theta',PatientID,...
    num_contact_configs,fig_name,stim_state,caxis_limits); 

t_alpha = generate_coherence_maps(coh_diff_alpha,names,'alpha',PatientID,...
    num_contact_configs,fig_name,stim_state,caxis_limits); 

