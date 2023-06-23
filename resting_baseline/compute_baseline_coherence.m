% turn into function once we know this works 
check_nan = @(x) assert(~any(isnan(x(:))));
check_nan(output);
    for i = 1:num_ch
        for j = i:5
            for itrial = 1:num_trials
               
                
                first_ch = squeeze(output(:,itrial,i));
                check_nan(first_ch);
                second_ch = squeeze(output(:,itrial,j));
                check_nan(second_ch)
                [cxy, freqs] = mscohere(first_ch,second_ch,...
                    hamming(win_len),noverlap,[],fs); 
                 check_nan(cxy);
                 check_nan(freqs); 
                                % get indices for frequencies 
                    %get from other script 
                    theta_idx = find(freqs >3 & freqs < 7.1); 
                    theta_freqs = freqs(theta_idx); 
                    theta_cxy = cxy(theta_idx);                   
                 
                % take avg of coh values across each freq band 
                
                theta_val = mean(theta_cxy); 
                theta_val_alltrials(itrial) = theta_val; 
                                
            end 
           
                % take average of coherence values across trials 
                                
               theta_mean = mean(theta_val_alltrials); 
               theta_std = std(theta_val_alltrials); 

                              
                % generate matrix of coherence values and std dev for each foi   
                matrix_coh_theta(i,j) = theta_mean;
                matrix_coh_theta(j,i) = theta_mean;
                matrix_err_theta(i,j) = theta_std; 
                matrix_err_theta(j,i) = theta_std; 

        end 
    end 
 