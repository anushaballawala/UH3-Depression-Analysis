function [S_prestim,f_prestim,...
    S_stim,f_stim, S_poststim, f_poststim] = compute_multitaper_spectrum_power(num_ch,params,...
    num_trials,indivtr_pre_stim_data,...
    indivtr_stim_data, indivtr_poststim_data)
S_prestim = []; 
S_stim = []; 
S_poststim = []; 
for i=1:num_ch
    fprintf('....Processing Channel %d \n',i)
    for j = 1:num_trials 
    %prestim 
    fprintf('....Processing Trial %d \n',j)
        [S_prestim(j,i,:),f_prestim,~] = mtspectrumc(squeeze(permute(indivtr_pre_stim_data(j,i,:),[3,1,2])),params);
        [S_stim(j,i,:),f_stim,~] = mtspectrumc(squeeze(permute(indivtr_stim_data(j,i,:),[3,1,2])),params);
        [S_poststim(j,i,:),f_poststim,~] = mtspectrumc(squeeze(permute(indivtr_poststim_data(j,i,:),[3,1,2])),params);
    end    
end
