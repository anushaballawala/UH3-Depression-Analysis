function [contact_configs,num_contact_configs] = output_contact_config_conditions(PatientID, DBS_target)

switch PatientID 
    case 'DBSTRD001' 
        switch DBS_target 
            case 'VCVS' 
                disp('VCVS')
                contact_configs = {'elec1','elec25','elec36','elec47','elec8'}; 
                num_contact_configs = numel(contact_configs); 
                
            case 'SCC' 
                disp('SCC')
                contact_configs = {'elec1','elec25','elec234','elec36','elec47','elec234','elec8'}; 
                num_contact_configs = numel(contact_configs); 

        end 
         
    case 'DBSTRD002'
        switch DBS_target 
            case 'VCVS' 
                disp('VCVS')
                contact_configs = {'elec1','elec25','elec234','elec36','elec47','elec234','elec8'}; 
                num_contact_configs = numel(contact_configs);
            case 'SCC'
                disp('SCC') 
                contact_configs = {'elec1','elec25','elec234','elec36','elec47','elec234','elec8'}; 
                num_contact_configs = numel(contact_configs);
        end 
        
    case 'DBSTRD003' 

        switch DBS_target 
            case 'VCVS' 
                disp('VCVS')
            case 'SCC'
                disp('SCC') 
                contact_configs = {'elec1','elec25','elec234','elec36','elec47','elec234','elec8'}; 
                num_contact_configs = numel(contact_configs);
        end 
         
end         
      

    end 