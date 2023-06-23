function [elec1,elec234,elec25,elec36,elec47,elec567,elec8] = assign_SCC_conditions(data,hemi)

switch hemi
    case 'l'         
        for i = 1:numel(data.output_and_metadata)
            
            if find(contains(data.output_and_metadata{1, i}{1, 2},'elec17')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                disp('elec1')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec18_19_20')) == 1
                elec234 = data.output_and_metadata{1, i}{1, 1};
                disp('elec234')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec18_21')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                disp('elec25')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec19_22')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                disp('elec36')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec20_23')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                disp('elec47')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec21_22_23')) == 1
                elec567 = data.output_and_metadata{1, i}{1, 1};
                disp('elec567')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec24')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                disp('elec8')
            end
        end 
            
   case 'r'         
    
                
                for i = 1:numel(data.output_and_metadata)
                    
                    if find(contains(data.output_and_metadata{1, i}{1, 2},'elec25')) == 1
                        elec1 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec1')
                    elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec26_27_28')) == 1
                        elec234 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec234')
                    elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec26_29')) == 1
                        elec25 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec25')
                    elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec27_30')) == 1
                        elec36 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec36')
                    elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec28_31')) == 1
                        elec47 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec47')
                    elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec29_30_31')) == 1
                        elec567 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec567')
                    elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec32')) == 1
                        elec8 = data.output_and_metadata{1, i}{1, 1};
                        disp('elec8')
                        
                    end
                end
                
        end
end 
 