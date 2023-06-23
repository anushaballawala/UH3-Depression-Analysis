% Define input parameters 

clear all 
PatientID = 'DBSTRD001'; 
DBS_target = 'VCVS'; 
hemi = 'r'; 
stim_freq = 130; 
experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq);  

% Load single trial time series data 
timeseriesfile = sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
    PatientID,hemi,DBS_target,stim_freq,hemi,DBS_target,stim_freq); 
timeseries_data = load(timeseriesfile); 
fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 
ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 

addpath(genpath('C:\CODE\Toolboxes\chronux_2_12')); 
%% Run function for all contact configurations 
disp(experiment_name) 
switch PatientID
    case 'DBSTRD001'
         ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lMTG',...
            'lVPF','rACC','rAMY','rDPF','rLOF','rMOF','rMTG','rVPF'}; 
        switch DBS_target
            case 'VCVS'
                % get different contact configuration data. 
                [elec1, elec25, elec36, elec47, ....
                    elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi);
                % assign names(string) for contact configs 
                all_config_names = {'elec1','elec25', 'elec36', 'elec47',...
                     'elec8'};
            case 'SCC'
              [elec1,elec234,elec25,elec36,elec47,...
                   elec567,elec8]  = assign_SCC_conditions(timeseries_data,hemi); 
        end
    case 'DBSTRD002'
          ROI_labels = {'lACC','lAMY','lDPF','lLOF','lMOF','lSFG',...
      'lSTG','lVPF','rACC','rAMY','rDPF','rLOF','rMOF',...
      'rPHG','rSFG','rSTG','rVPF'}; 
        switch DBS_target
            case 'VCVS'
                % get different contact configuration data 
                [elec1, elec234, elec25, elec36, elec47,...
                    elec567, elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi);
                % assign names (String) for contact configs 
                all_config_names = {'elec1','elec234', 'elec25', 'elec36', 'elec47',...
                    'elec567', 'elec8'};
            case 'SCC' 
               [elec1,elec234,elec25,elec36,elec47,...
                   elec567,elec8]  = assign_SCC_conditions(timeseries_data,hemi); 
                
        end
    case 'DBSTRD003'
        switch DBS_target
            case 'VCVS'
                 % get different contact configuration data 
                [elec1, elec234, elec25, elec36, elec47,...
                    elec567, elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi);
                 % assign names (String) for contact configs 
                all_config_names = {'elec1','elec234', 'elec25', 'elec36', 'elec47',...
                    'elec567', 'elec8'};
            case 'SCC'
                 % get different contact configuration data 
                [elec1, elec234, elec25, elec36, elec47,...
                    elec567, elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi);
                 % assign names (String) for contact configs 
                all_config_names = {'elec1','elec234', 'elec25', 'elec36', 'elec47',...
                    'elec567', 'elec8'};
        end
end

%% 

num_contact_configs = numel(timeseries_data.output_and_metadata); 
disp(num_contact_configs);
%% 

for i = 1:num_contact_configs
    contact_config = all_config_names{i}; 
    switch contact_config
        case 'elec25'
            data = elec25;
        case 'elec1'
            data = elec1;
        case 'elec8'
            data = elec8;
        case 'elec36'
            data = elec36;
        case 'elec47'
            data = elec47;
        case 'elec234'
            data = elec234;
        case 'elec567'
            data = elec567;
    end 
    disp(all_config_names{i}); 
    names{i} = all_config_names{i}; 
    
    plot_1overf_ON_OFF_PSD_ROI_indivch(PatientID, contact_config, DBS_target, hemi, stim_freq,...
        experiment_name, timeseriesfile, timeseries_data, fs, ch_labels,...
        ROI_labels, data); 
 
    
end 
