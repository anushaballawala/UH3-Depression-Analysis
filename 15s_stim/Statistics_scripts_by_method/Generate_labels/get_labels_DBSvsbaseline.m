%% 
function [labels] = get_labels_DBSvsbaseline(comparison, num_trials, stim_cond, tbl, elec_cond_type,...
    cond1, cond2, varargin)

%%%%%%% 
%Function get_labels_DBSvsbaseline.m 

% The goal of this fxn is to generate binary labels when statistically 
% testing b/w trials where stim was on, or poststim, vs baseline recording.
% You can test baseline vs all DBS targets, or a specific DBS target, as
% well as against all electrode configurations (ring, single, etc.) or 
% against all of them together. 
% Support is not included for combinations of electrode configurations or 
% stimulation frequency. 

%%%%%%%
%% 

%check that correct function has been chosen for analysis: 
if strcmp(comparison,'DBSvsbaseline') == 1 
    disp('generating labels for DBS vs baseline')
elseif strcmp(comparison,'DBSvsbaseline') == 0 
    error('check the conditions chosen for comparison') 
end 

labels = zeros(num_trials,1); %initialize. 
if nargin == 1 
    elec_config = varargin{1};
else 
    disp('no elec config input')
end 

%test between all electrode configs or a specific one. 
switch elec_cond_type
    case 'all'  % all elec stim configurations
        
        % test between stim states 
        switch stim_cond            
            case 'stimall'  % all stim states where stim was on
                for i = 1:num_trials
                    if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim1') == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim2') == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'stim3') == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond2) == 1
                        labels(i,1) = 0;
                    else
                        labels(i,1) = NaN;
                    end
                end
                
            case 'poststim' % all poststim trials for all elec stim configurations
                for i = 1:num_trials
                    if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond2)== 1
                        labels(i,1) = 0;
                    else
                        labels(i,1) = NaN;
                    end
                end
        end
        
    case 'indiv' %trials for a specific electrode configuration
        switch stim_cond
            case 'stimall' %all trials where stim was on
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
                
            case 'poststim' % all poststim trials for specific electrode configuration
                for i = 1:num_trials
                    if contains(tbl.DBS(i),cond1)==1 && strcmp(tbl.Stim_state(i),'poststim') == 1 && strcmp(tbl.elec(i),elec_config) == 1
                        labels(i,1) = 1;
                    elseif contains(tbl.DBS(i),cond2)== 1
                        labels(i,1) = 0;
                    else
                        labels(i,1) = NaN;
                    end
                end
        end
        
        
end
end 


