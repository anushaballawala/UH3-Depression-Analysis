%combine lSCC stim on and off data 



%% create artificial trials 


for i = 1:7 
    delta{i} = repmat(matrix_coh_delta{i},1,1,5);
    theta{i} = repmat(matrix_coh_theta{i},1,1,5);
    alpha{i} = repmat(matrix_coh_alpha{i},1,1,5); 
    beta{i} = repmat(matrix_coh_beta{i},1,1,5);
    gamma{i} = repmat(matrix_coh_gamma{i},1,1,5);
    
end 


%% combine all current directions 
all_alpha = cat(3,alpha{:}); 
all_theta = cat(3,theta{:}); 
all_beta = cat(3,beta{:}); 
all_delta = cat(3,delta{:}); 
all_gamma = cat(3,gamma{:}); 

%% concatenate freqbands 

all_data = cat(4,all_alpha,all_theta,all_beta,all_delta); 

%% labels 

labels1 = repmat(1,17,1); 
labels2 = repmat(0,18,1); 

labels = [labels1;labels2]; 







