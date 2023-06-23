function [t_stat] = compute_Tstatistic_coh(num_ch_roi,num_freq,data,shuffled_labels)

stats = []; 

    for i = 1:num_ch_roi
        for k = i:num_ch_roi 
            for j = 1:num_freq
                [~,~,~,stats] = ttest2(data(shuffled_labels==1,i,k,j),...
                    data(shuffled_labels==0,i,k,j),'Tail','both');
                t_stat(i,k,j) = stats.tstat;
                if isnan(t_stat(i,k,j)), t_stat(i,k,j)=0; end
            end
        end 
    end
    
    
end


