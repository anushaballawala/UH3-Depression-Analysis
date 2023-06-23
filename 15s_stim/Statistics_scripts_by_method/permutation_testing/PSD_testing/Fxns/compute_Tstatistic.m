function [t_stat] = compute_Tstatistic(num_ch_roi,num_freq,data,shuffled_labels)

stat = []; 

    for i = 1:num_ch_roi
        for j = 1:num_freq
            [~,~,~,stats(i,j)] = ttest2(data(shuffled_labels==1,i,j),...
                data(shuffled_labels==0,i,j),'Tail','both');
            t_stat(i,j) = stats(i,j).tstat;
        end
    end
end

        
       