% Define input parameters 

function[] = run_timeseries_mtspectrumc_power(timeseriesfile)

% timeseriesfile = dir('/gpfs/data/dborton/TRD_Project/DBSTRD/DBSTRD003/Experiments/15s_stim/Epoched Data/*/15s_stim_all_currdir_timeseries_singletrial_VCVS_f130_567.mat'); 
% load(sprintf('%s/%s',timeseriesfile.folder,timeseriesfile.name)); 
% timeseriesfile = sprintf('%s/%s',timeseriesfile.folder,timeseriesfile.name)
%% 
% get info from filename 
subjectID_pattern= '/DBSTRD00*\w*'; 
[start_idx_subid end_idx_subid] = regexp(timeseriesfile,subjectID_pattern);
PatientID = extractBetween(timeseriesfile,start_idx_subid+1,end_idx_subid);
PatientID = PatientID{1}; 

braintarget = '/(l|r)(VCVS|SCC)_';
[start_idx_br, end_idx_br] = regexp(timeseriesfile,braintarget);
DBS_target = extractBetween(timeseriesfile,start_idx_br+1,end_idx_br-1);
DBS_target = DBS_target{1}; 

hemi = DBS_target(1); 

DBS_target = DBS_target(2:end); 
stim_freq = 130; 

experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq);  

% % Load single trial time series data 
% timeseriesfile = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Epoched Data/%s%s_f%d/15s_stim_all_currdir_timeseries_singletrial_%s%s_f%d.mat',...
%     PatientID,hemi,DBS_target,stim_freq,hemi,DBS_target,stim_freq); 


%% Load data 
timeseries_data = load(timeseriesfile); 

if strcmp(experiment_name,'rVCVS_f130')==1 && strcmp(PatientID,'DBSTRD002') ==1 
    fs = timeseries_data.metadata.metadata_file1.preprocessing.New_SamplingRate;
    ch_labels = deblank(timeseries_data.metadata.metadata_file1.preprocessing.GoodChannelLabels);
else
    fs = timeseries_data.metadata.preprocessing.New_SamplingRate;
    ch_labels = deblank(timeseries_data.metadata.preprocessing.GoodChannelLabels);
end

% Add chronux to filepath.
addpath(genpath('/gpfs/home/aallawa1/Documents/MATLAB/Toolboxes/chronux_2_12.v03')); 

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
               all_config_names = {'elec1','elec234','elec25','elec36','elec47',...
                   'elec567','elec8'}; 
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
               all_config_names = {'elec1','elec234','elec25','elec36','elec47',...
                   'elec567','elec8'}; 
                
        end
    case 'DBSTRD003'
        switch DBS_target
            case 'VCVS'
                 % get different contact configuration data 
                [elec1, elec2] = assign_VCVS_conditions(PatientID,timeseries_data,hemi);
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

%% Metadata 

metadata = timeseries_data.metadata; 
metadata.processing.date = date; 

%% Run multitaper spectrum fxn 
num_contact_configs = numel(all_config_names); 
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
    
[pow_prestim,pow_stim1,pow_stim2,pow_stim3,...
        pow_poststim,f_prestim, S_prestim, f_stim1, S_stim1,...
        S_stim2,S_stim3,f_poststim,S_poststim,Serr_prestim, Serr_stim1, Serr_stim2,...
        Serr_stim3, Serr_poststim,metadata] = mtspectrum_pow_from_timeseries(PatientID,...
        contact_config, DBS_target, hemi, stim_freq,experiment_name,...
                timeseriesfile, timeseries_data, fs, ch_labels,...
                data,metadata); 

outputdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/Experiments/15s_stim/Processed Data/PSD',...
    PatientID);
 if ~exist(outputdir,'dir')
     mkdir(outputdir); end 
 
    filename = sprintf('%s_%s_mtspectrum_pow_indivtr_elec%d.mat',experiment_name,contact_config); 
    filedir = sprintf('%s/%s',outputdir,filename); 
    save(filedir,'pow_prestim','pow_stim1','pow_stim2','pow_stim3',...
        'pow_poststim','f_prestim', 'S_prestim', 'f_stim1', 'S_stim1',...
        'S_stim2','S_stim3','f_poststim','S_poststim','Serr_prestim', 'Serr_stim1', 'Serr_stim2',...
        'Serr_stim3', 'Serr_poststim','metadata') 


             

end 
end 
