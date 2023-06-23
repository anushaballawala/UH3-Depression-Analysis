function [matrix_coh_delta,...
    matrix_coh_theta,...
    matrix_coh_alpha,...
    matrix_coh_beta,...
    matrix_coh_gamma,...
    delta_freqs,theta_freqs,alpha_freqs,...
    beta_freqs,gamma_freqs] = compute_coherence_single_trial(num_channels,data,win_len,noverlap,f,fs)
  
    %Compute coherence across all channels and trials. 
    for i = 1:num_channels
        for j = i:num_channels
           
                
                first_ch = squeeze(data(:,i));
                second_ch = squeeze(data(:,j));
                [cxy, freqs] = mscohere(first_ch,second_ch,...
                    hamming(win_len),noverlap,f,fs); 
                
                % get indices for frequencies 
                    %get from other script 
                    theta_idx = find(freqs >3 & freqs < 7.1); 
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
                    
                    gamma_idx = find(freqs >30 & freqs <80); 
                    gamma_freqs = freqs(gamma_idx); 
                    gamma_cxy = cxy(gamma_idx); 
                    
                % take avg of coh values across each freq band 
                
                theta_val = mean(theta_cxy); 
                theta_val_alltrials = theta_val; 
                
                delta_val = mean(delta_cxy); 
                delta_val_alltrials = delta_val; 
                
                alpha_val = mean(alpha_cxy); 
                alpha_val_alltrials = alpha_val; 
                
                beta_val = mean(beta_cxy); 
                beta_val_alltrials = beta_val; 
                
                gamma_val = mean(gamma_cxy); 
                gamma_val_alltrials = gamma_val; 
                                
        end 
          
                                              
                % generate matrix of coherence values for each foi   
                matrix_coh_theta(i,j) = theta_val;
                matrix_coh_theta(j,i) = theta_val;

                
                matrix_coh_delta(i,j) = delta_val; 
                matrix_coh_delta(j,i) = delta_val; 

                
                matrix_coh_alpha(i,j) = alpha_val; 
                matrix_coh_alpha(j,i) = alpha_val; 
  
                
                matrix_coh_beta(i,j) = beta_val; 
                matrix_coh_beta(j,i) = beta_val; 
   
                
                matrix_coh_gamma(i,j) = gamma_val; 
                matrix_coh_gamma(j,i) = gamma_val; 

                
    end
   
end 
