function [elec1,elec25,elec36,elec47,elec8] = assign_VCVS_conditions(data,hemi)
% 
switch hemi
    case 'l' 
        for i = 1:numel(data.output_and_metadata)
         if find(contains(data.output_and_metadata{1, i}{1, 2},'elec1')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                disp('elec1')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec2_5')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                disp('elec25')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec3_6')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                disp('elec36')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec4_7')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                disp('elec47')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec8')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                disp('elec8')
         end
        end 

              
    case 'r' 
        for i = 1: numel(data.output_and_metadata) 
        if find(contains(data.output_and_metadata{1, i}{1, 2},'elec9')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                disp('elec1')

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec10_13')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                disp('elec25')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec11_14')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                disp('elec36')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec12_15')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                disp('elec47')

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec16')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                disp('elec8')
        end
        end 
end 
end 
