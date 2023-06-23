function [labels] = get_labels(comparison, num_trials, stim_cond, tbl, cond_type,...
cond1, cond2, varargin)

elec_type = varargin{1} ; 
elec_config = varargin{2} ;  


labels = zeros(num_trials,1); 
switch comparison   
    case 'targetvstarget'     
        switch elec_type             
            case 'all' 
                switch stim_cond 
                    case 'stim on'
                        for i = 1:num_trials
                            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim1') == 1
                                labels(i,1) = 1;
                            elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim2') == 1
                                labels(i,1) = 1;
                            elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim3') == 1
                                labels(i,1) = 1;
                            elseif contains(tbl.DBS(i),cond2)==1 && strcmp(tbl.Stim_state(i),'stim1') == 1
                                labels(i,1) = 0;
                            elseif contains(tbl.DBS(i),cond2)==1 && strcmp(tbl.Stim_state(i),'stim2') == 1
                                labels(i,1) = 0;
                            elseif contains(tbl.DBS(i),cond2)==1 && strcmp(tbl.Stim_state(i),'stim3') == 1
                                labels(i,1) = 0;
                            else
                                labels(i,1) = NaN;
                            end
                        end
                    case 'poststim'                                              
                        for i = 1:num_trials
                            if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1
                                labels(i,1) = 1;
                            elseif contains(tbl.DBS(i),cond2)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1
                                labels(i,1) = 0;
                            else
                                labels(i,1) = NaN;
                            end
                        end
                end 
                        
            case 'indiv'
                switch stim_state 
                    case 'stim on' 
                        for i = 1:num_trials
                            if contains(DBS_labels(i),cond1)==1 && strcmp(stim_labels(i),'stim1') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 1;
                            elseif contains(DBS_labels(i),cond1)==1 && strcmp(stim_labels(i),'stim2') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 1;
                            elseif contains(DBS_labels(i),cond1)==1 && strcmp(stim_labels(i),'stim3') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 1;
                            elseif contains(DBS_labels(i),cond2)==1 && strcmp(stim_labels(i),'stim1') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 0;
                            elseif contains(DBS_labels(i),cond2)==1 && strcmp(stim_labels(i),'stim2') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 0;
                            elseif contains(DBS_labels(i),cond2)==1 && strcmp(stim_labels(i),'stim3') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 0;
                            else
                                labels(i,1) = NaN;
                            end
                        end
                    case 'poststim'                        
                        for i = 1:num_trials
                            if contains(DBS_labels(i),cond1)==1 && strcmp(stim_labels(i),'poststim') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 1;
                            elseif contains(DBS_labels(i),cond2)==1 && strcmp(stim_labels(i),'poststim') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                                labels(i,1) = 0;
                            else
                                labels(i,1) = NaN;
                            end
                        end                        
                                             
                end 
        end 
    case 'baselinevsDBStarget'
        switch cond_type 
            case 'all'
            for i = 1:num_trials
                if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim1') == 1
                    labels(i,1) = 1;
                elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim2') == 1
                    labels(i,1) = 1;
                elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim3') == 1
                    labels(i,1) = 1;
                elseif contains(tbl.DBS(i),cond2)== 1 
                    labels(i,1) = 0;
                else
                    labels(i,1) = NaN;
                end
            end
            case 'indiv' 
                for i = 1:num_trials
                    if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim1') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim2') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim3') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond2)== 1
                        labels(i,1) = 0;
                    else
                        labels(i,1) = NaN;
                    end
                end
        end 

    case 'stimstate' 

        switch cond1
            case 'prestim'
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),cond1)==1 
                        labels(i,1) = 1;
                    end 
                end 
            case 'Baseline' 
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),cond1)==1 
                        labels(i,1) = 1;
                    end
                end
        end 
                
        switch cond2
            case 'stim1' 
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),cond2)==1
                        labels(i,1) = 0;
                    end
                end
            case 'stim2'
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),cond2)==1
                        labels(i,1) = 0;
                    end
                end
            case 'stim3'
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),cond2)==1
                        labels(i,1) = 0;
                    end
                end
            case 'poststim'
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),cond2)==1
                        labels(i,1) = 0;
                    end
                end
            case 'stimall'
                for i = 1:num_trials
                    if contains(tbl.Stim_state(i),'stim1')==1 || contains(tbl.Stim_state(i),'stim2')==1 || contains(tbl.Stim_state(i),'stim3')==1
                        labels(i,1) = 0;
                    end
                end
        end 
                
        
end 

end 