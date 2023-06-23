function [t_stat,shuffled_labels] = permutations_multiple_cond(data,labels,num_perm,fxn_handle) 

%output t_stat is array of M X num_ch X freqs

M = num_perm+1; % first permutation is just actual data 

shuffled_labels = zeros(numel(labels),M);
%t_stat = cell(num_perm+1,1); 

shuffled_labels(:,1) = labels; 

%do permutations 
    for k = 2:M
            shuffled_labels(:,k) = labels(randperm(length(labels)));
            disp(k)
    end 
    
t = fxn_handle(data,shuffled_labels(:,1)); % ch x num freqbands 


t_stat = zeros([M,size(t)]); 
t_stat(1,:) = reshape(t,1,numel(t)); 


    %get t-statistic for shuffled data 
    for k = 2:num_perm+1
            %t_stat(k) = fxnhandle(data(labels==1),data(labels==0),shuffled_labels(:,k),'Ch') ; 
            t = fxn_handle(data,shuffled_labels(:,k));
            t_stat(k,:) = reshape(t,1,numel(t)); % M x numch x num freqbands x num_SCC conds 
            disp(k)
    end
    

end


