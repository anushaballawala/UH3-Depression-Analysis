% Define input parameters 
function[timeseriesfile, timeseries_data, fs, ch_labels,...
    all_config_names,metadata,...
    elec1,elec25,elec36,elec47,...
    elec8,varargout] = run_setup_15s(PatientID, hemi, DBS_target, stim_freq, experiment_name,filename) %* commented out rings 


num_inputs = 11; %if number of inputs to fxn is larger change this num 
    %**CHANGE 
%timeseriesfile = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
%     PatientID,hemi,DBS_target,stim_freq,hemi,DBS_target,stim_freq); 
 

timeseriesfile = filename; 
timeseries_data = load(filename); 

if strcmp(DBS_target,'VCVS') ==1 && strcmp(hemi,'r') ==1 && strcmp(PatientID, 'DBSTRD002') ==1 
    timeseries_data.metadata = timeseries_data.metadata.metadata_file1; 
else 
    disp('do nothing')
end 

fs = timeseries_data.metadata.preprocessing.New_SamplingRate; 

ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels); 

%get metadata & initialize for processing 
metadata = timeseries_data.metadata; 
metadata.processing = []; 
%% Run function for all contact configurations 


disp(experiment_name) 
switch PatientID
    case 'DBSTRD001'

        switch DBS_target            
            case 'VCVS'
                num_outputs_VCVScond = 5;
                % get different contact configuration data. 
                [elec1, elec25, elec36, elec47, ....
                    elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi,num_outputs_VCVScond);
                % assign names(string) for contact configs 
                all_config_names = {'elec1','elec25', 'elec36', 'elec47',...
                     'elec8'};                 
                 
            case 'SCC'

              [elec1,elec234,elec25,elec36,elec47,...
                   elec567,elec8]  = assign_SCC_conditions(timeseries_data,hemi); 
               all_config_names = {'elec1','elec234','elec25','elec36','elec47',...
                   'elec567','elec8'}; 

               if nargout > num_inputs
                     varargout{1} = elec234;
                     varargout{2} = elec567;
              end
        end
    case 'DBSTRD002'

        switch DBS_target
            case 'VCVS' 
                num_outputs_VCVScond = 7;
                % get different contact configuration data 
                [elec1, elec25, elec36, elec47, elec8,elec234,elec567] = assign_VCVS_conditions(PatientID,timeseries_data,hemi,num_outputs_VCVScond);
                % assign names (String) for contact configs 
                all_config_names = {'elec1','elec25', 'elec36', 'elec47', 'elec8',...
                    'elec234', 'elec567'};

               
                     varargout{1} = elec234;
                     varargout{2} = elec567;
                 
                 
                 
            case 'SCC' 
               [elec1,elec234,elec25,elec36,elec47,...
                   elec567,elec8]  = assign_SCC_conditions(timeseries_data,hemi);
               all_config_names = {'elec1','elec234','elec25','elec36','elec47',...
                   'elec567','elec8'}; 

                if nargout > num_inputs
                     varargout{1} = elec234;
                     varargout{2} = elec567;
                end
                
        end
    case 'DBSTRD003'
        switch DBS_target
            case 'VCVS'
                 % get different contact configuration data 
                [elec1, elec234, elec25, elec36, elec47,...
                    elec567, elec8] = assign_VCVS_conditions(PatientID,timeseries_data,hemi,num_outputs_VCVScond);
                 % assign names (String) for contact configs 
                all_config_names = {'elec1','elec234', 'elec25', 'elec36', 'elec47',...
                    'elec567', 'elec8'};
            case 'SCC'
                 % get different contact configuration data 
                [elec1, elec234, elec25, elec36, elec47,...
                    elec567, elec8] = assign_SCCconditions(PatientID,timeseries_data,hemi);
                 % assign names (String) for contact configs 
                all_config_names = {'elec1','elec234', 'elec25', 'elec36', 'elec47',...
                    'elec567', 'elec8'};
        end
end

%% 
num_contact_configs = numel(timeseries_data.output_and_metadata); 
fprintf('Number of contact configurations = %d',num_contact_configs);


end