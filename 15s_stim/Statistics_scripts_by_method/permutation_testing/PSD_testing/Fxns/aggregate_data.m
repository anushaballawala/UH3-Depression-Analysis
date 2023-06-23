function [new_data,new_data_labels,data1,data2,...
    data1_idx, data2_idx] = aggregate_data(data,labels,stimtype)

%Aggregate data - anything with 1 or 0, etc. will be permuted together 
pre_stim_idx = find(labels == -1); 
stim1_idx = find(labels == 1.1); 
stim2_idx = find(labels == 1.2); 
stim3_idx = find(labels == 1.3); 
poststim_idx = find(labels == 2); 
allstim_idx = find(round(labels) == 1); 
baseline_idx = find(labels == 0); 

switch stimtype 
    case 'prestimvsstimall'
        data1 = data(pre_stim_idx,:,:); data2 = data(allstim_idx,:,:); 
        new_labels1 = labels(pre_stim_idx,:); new_labels2 = labels(allstim_idx,:); 
        disp('comparing prestim and stim all')
    case 'prestimvsstim1'
        data1 = data(pre_stim_idx,:,:); data2 = data(stim1_idx,:,:); 
        new_labels1 = labels(pre_stim_idx,:); new_labels2 = labels(stim1_idx,:); 
        disp('comparing pre-stim and stim1')
    case 'prestimvspoststim'
        data1 = data(pre_stim_idx,:,:); data2 = data(poststim_idx,:,:); 
        new_labels1 = labels(pre_stim_idx,:); new_labels2 = labels(poststim_idx); 
        disp('comparing prestim and poststim')
    case 'baselinevspoststim'
        data1 = data(baseline_idx,:,:); data2 = data(poststim_idx,:,:); 
        new_labels1 = labels(baseline_idx,:); new_labels2 = labels(poststim_idx); 
        disp('comparing poststim and baseline') 
end 
        
new_data = vertcat(data1,data2); 
new_data_labels = vertcat(new_labels1, new_labels2); 
data1_idx = 1:size(data1,1); 
data2_idx = size(data1,1)+1:size(data1,1)+size(data2,1) ; 

end 