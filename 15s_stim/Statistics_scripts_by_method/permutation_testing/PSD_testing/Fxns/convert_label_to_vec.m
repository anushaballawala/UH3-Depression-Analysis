%% convert stim/presstim state labels to vector 
function[stim_vec_labels] = convert_label_to_vec(stim_labels)

% Input = labels containing info about stim state 

% Output = vector corresponding to stim state labels 

stim_vec_labels = zeros(size(stim_labels)); 


for i = 1:length(stim_vec_labels) 
    switch string(stim_labels(i))

        case 'prestim'
            stim_vec_labels(i) = -1 ; 
            
        case 'stim1' 
            stim_vec_labels(i) = 1.1; 
            
        case 'stim2'
            stim_vec_labels(i) = 1.2; 
            
        case 'stim3'
            stim_vec_labels(i) = 1.3; 
            
        case 'poststim'
            stim_vec_labels(i) = 2; 
            
        case 'Baseline'
            stim_vec_labels(i) = 0; 
            
    end 
end 

end 






