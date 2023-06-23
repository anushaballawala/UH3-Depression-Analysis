function [output_and_metadata,current_condition,conditions,tbl_idx,num_conditions] = bin_time_series_data(measured_data,tbl,fs,length_time,trial_length,stim_win,num_channels) 
    
    conditions = unique(tbl.Condition); % An array of the unique condition summaries in the column.
    condition_bytrial = tbl.Condition;
    num_conditions = length(conditions); %num of stim conditions tested in this file 

    output = cell(num_conditions,1); % each cell in cell array is data for one condition
    output_and_metadata = {}; 
    tbl_idx = {}; 
    output_and_metadata = {}; 

       % Process one condition bin at a time.
    for condition_idx = 1:num_conditions
        current_condition = conditions(condition_idx); % Get the condition numbered `condition_idx` which is an element of the array conditions;
        tbl_for_condition = tbl((tbl.Condition == current_condition),:); 
        tbl_idx{condition_idx} = tbl_for_condition; % may throw error ****
        num_trials_per_condition = height(tbl_for_condition); 
        
        % Pull out only the tbl rows for the current condition.
        output_for_condition = zeros(num_trials_per_condition, num_channels,trial_length);

        % Pull out one trial at a time.
        for trial_idx = 1:num_trials_per_condition
            % Time that trial started.
            start_time = tbl_for_condition.Time(trial_idx);
            start_idx = (start_time-3)*fs;
            start_idx = round(start_idx); % Need to make sure we get an integer
              %%%%%%% This code snippet needs to fixed to address when %%%
              % pre-stim is too short at start of recording 
%             %if start idx negative, remove first timestamp 
%             if start_idx <0 
%                 tbl_for_condition(1,:) = []; 
%             elseif 
%                 continue % continue as is 
%             end 
            %%%%%%%% %%%%%%%%%
            end_idx = start_idx + trial_length - 1;
            if end_idx > length_time
                error('Trial is asking for data after the end of the session (check sample rate)');
            end
            all_channels_for_this_trial = measured_data(:,start_idx:end_idx);
            output_for_condition(trial_idx, :, :) = all_channels_for_this_trial;
        end
        output{condition_idx} = output_for_condition; 
        output_and_metadata{condition_idx} = [{output_for_condition},{current_condition}];
        fprintf('finished generating epoch for %s', current_condition)

    end
    disp('output cell created with diff conditions'); 
    
 
end 