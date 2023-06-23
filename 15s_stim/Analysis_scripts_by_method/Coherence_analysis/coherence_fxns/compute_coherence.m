function [matrix_coh_delta,matrix_err_delta,...
    matrix_coh_theta,matrix_err_theta,...
    matrix_coh_alpha,matrix_err_alpha,...
    matrix_coh_beta,matrix_err_beta,...
    matrix_coh_lowgamma, matrix_err_lowgamma,matrix_coh_highgamma,matrix_err_highgamma,...
    delta_freqs,theta_freqs,alpha_freqs,...
    beta_freqs,lowgamma_freqs,highgamma_freqs,...
    delta_indiv_tr,theta_indiv_tr,alpha_indiv_tr,...
    beta_indiv_tr,lowgamma_indiv_tr,highgamma_indiv_tr] = compute_coherence(num_channels,num_trials,data,win_len,noverlap,f,fs)
  %% 
    %initialize coherence matrices 
    matrix_coh_delta = zeros(num_channels, num_channels);
    matrix_coh_theta = zeros(num_channels, num_channels);
    matrix_coh_alpha = zeros(num_channels, num_channels);
    matrix_coh_beta = zeros(num_channels, num_channels);
    matrix_coh_lowgamma = zeros(num_channels, num_channels);
    matrix_coh_highgamma = zeros(num_channels, num_channels);

    %initialize std error matrices 
    matrix_err_delta  = zeros(num_channels, num_channels);
    matrix_err_theta = zeros(num_channels, num_channels);
    matrix_err_alpha  = zeros(num_channels, num_channels);
    matrix_err_beta  = zeros(num_channels, num_channels);
    matrix_err_lowgamma  = zeros(num_channels, num_channels);
    matrix_err_highgamma = zeros(num_channels, num_channels);
    
    %Compute coherence across all channels and trials. 
    for i = 1:num_channels
        for j = i:num_channels
            for itrial = 1:num_trials
                
                first_ch = squeeze(data(:,itrial,i)); 
                second_ch = squeeze(data(:,itrial,j)); 
                [cxy, freqs] = mscohere(first_ch,second_ch,...
                    hamming(win_len),noverlap,[],fs); 
                
                % get indices for frequencies 
                    %get from other script 
                    theta_idx = find(freqs >3 & freqs < 8); 
                    theta_freqs = freqs(theta_idx); 
                    theta_cxy = cxy(theta_idx); 
                    
                    delta_idx = find(freqs > 1 & freqs < 4); 
                    delta_freqs = freqs(delta_idx); 
                    delta_cxy = cxy(delta_idx); 
                    
                    alpha_idx = find(freqs > 8 & freqs < 13); 
                    alpha_freqs = freqs(alpha_idx); 
                    alpha_cxy = cxy(alpha_idx); 
                    
                    beta_idx = find(freqs > 13 & freqs < 30); 
                    beta_freqs = freqs(beta_idx); 
                    beta_cxy = cxy(beta_idx); 
                    
                    lowgamma_idx = find(freqs >30 & freqs <70); 
                    lowgamma_freqs = freqs(lowgamma_idx); 
                    lowgamma_cxy = cxy(lowgamma_idx); 
                    
                    highgamma_idx = find(freqs >70 & freqs <150); 
                    highgamma_freqs = freqs(highgamma_idx); 
                    highgamma_cxy = cxy(highgamma_idx);                 
                    
                % take avg of coh values across each freq band 
                
                theta_val = mean(theta_cxy); 
                theta_val_alltrials(itrial) = theta_val; 
                
                theta_indiv_tr(itrial,i,j) = theta_val_alltrials(itrial); 
                theta_indiv_tr(itrial,j,i) = theta_val_alltrials(itrial); 

                delta_val = mean(delta_cxy); 
                delta_val_alltrials(itrial) = delta_val; 
                
                delta_indiv_tr(itrial,i,j) = delta_val_alltrials(itrial); 
                delta_indiv_tr(itrial,j,i) = delta_val_alltrials(itrial); 
                
                alpha_val = mean(alpha_cxy); 
                alpha_val_alltrials(itrial) = alpha_val; 
                
                alpha_indiv_tr(itrial,i,j) = alpha_val_alltrials(itrial); 
                alpha_indiv_tr(itrial,j,i) = alpha_val_alltrials(itrial); 
                
                beta_val = mean(beta_cxy); 
                beta_val_alltrials(itrial) = beta_val; 
                
                beta_indiv_tr(itrial,i,j) = beta_val_alltrials(itrial); 
                beta_indiv_tr(itrial,j,i) = beta_val_alltrials(itrial); 
                
                lowgamma_val = mean(lowgamma_cxy); 
                lowgamma_val_alltrials(itrial) = lowgamma_val; 
                
                lowgamma_indiv_tr(itrial,i,j) = lowgamma_val_alltrials(itrial); 
                lowgamma_indiv_tr(itrial,j,i) = lowgamma_val_alltrials(itrial); 
                
                highgamma_val = mean(highgamma_cxy); 
                highgamma_val_alltrials(itrial) = highgamma_val;  
                
                highgamma_indiv_tr(itrial,i,j) = highgamma_val_alltrials(itrial); 
                highgamma_indiv_tr(itrial,j,i) = highgamma_val_alltrials(itrial); 
                                
            end 
           
                % take average of coherence values across trials 
               
               
               theta_mean = mean(theta_val_alltrials); 
               theta_std = std(theta_val_alltrials); 
               
               delta_mean = mean(delta_val_alltrials); 
               delta_std = std(delta_val_alltrials); 
               
               alpha_mean = mean(alpha_val_alltrials); 
               alpha_std = std(alpha_val_alltrials); 
               
               beta_mean = mean(beta_val_alltrials); 
               beta_std = std(beta_val_alltrials); 
               
               lowgamma_mean = mean(lowgamma_val_alltrials); 
               lowgamma_std = std(lowgamma_val_alltrials); 
               
               highgamma_mean = mean(highgamma_val_alltrials); 
               highgamma_std = std(highgamma_val_alltrials);              
                               
                 
                % generate matrix of coherence values and std dev for each foi   
                matrix_coh_theta(i,j) = theta_mean;
                matrix_coh_theta(j,i) = theta_mean;
                
                matrix_err_theta(i,j) = theta_std; 
                matrix_err_theta(j,i) = theta_std; 
                
                matrix_coh_delta(i,j) = delta_mean; 
                matrix_coh_delta(j,i) = delta_mean; 
                matrix_err_delta(i,j) = delta_std; 
                matrix_err_delta(j,i) = delta_std; 
                
                matrix_coh_alpha(i,j) = alpha_mean; 
                matrix_coh_alpha(j,i) = alpha_mean; 
                matrix_err_alpha(i,j) = alpha_std; 
                matrix_err_alpha(j,i) = alpha_std; 
                
                matrix_coh_beta(i,j) = beta_mean; 
                matrix_coh_beta(j,i) = beta_mean; 
                matrix_err_beta(i,j) = beta_std ; 
                matrix_err_beta(j,i) = beta_std; 
                
                matrix_coh_lowgamma(i,j) = lowgamma_mean; 
                matrix_coh_lowgamma(j,i) = lowgamma_mean; 
                matrix_err_lowgamma(i,j) = lowgamma_std; 
                matrix_err_lowgamma(j,i) = lowgamma_std; 
                
                matrix_coh_highgamma(i,j) = highgamma_mean; 
                matrix_coh_highgamma(j,i) = highgamma_mean; 
                matrix_err_highgamma(i,j) = highgamma_std; 
                matrix_err_highgamma(j,i) = highgamma_std;          
                
        end
    end
end 
