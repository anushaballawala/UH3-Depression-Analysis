
function [stim_e1,poststim1_e1,poststim2_e1,...
    stim_e234,poststim1_e234,poststim2_e234,...
    stim_e25,poststim1_e25,poststim2_e25,...
    stim_e36,poststim1_e36,poststim2_e36,...
    stim_e47,poststim1_e47,poststim2_e47,...
    stim_e567,poststim1_e567,poststim2_e567,...
    stim_e8,poststim1_e8,poststim2_e8,...
    stim_e, poststim1_e, poststim2_e] = compute_for_each_contact_config(DBStarget,data,freqs) 
switch DBStarget 
    case 'lSCC' 
    for i = 1:numel(data.output_and_metadata)

        if find(contains(data.output_and_metadata{1, i}{1, 2},'elec17')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                [stim_e1,poststim1_e1,poststim2_e1] = perform_single_tr_norm(elec1,freqs);
                disp('finished elec1')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec18_19_20')) == 1
                elec234 = data.output_and_metadata{1, i}{1, 1};
                [stim_e234,poststim1_e234,poststim2_e234] = perform_single_tr_norm(elec234,freqs);
                disp('elec234')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec18_21')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                [stim_e25,poststim1_e25,poststim2_e25] = perform_single_tr_norm(elec25,freqs);
                disp('elec25')
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec19_22')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                [stim_e36,poststim1_e36,poststim2_e36] = perform_single_tr_norm(elec36,freqs);
                disp('elec36')

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec20_23')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                [stim_e47,poststim1_e47,poststim2_e47] = perform_single_tr_norm(elec47,freqs);
                disp('elec47')

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec21_22_23')) == 1
                elec567 = data.output_and_metadata{1, i}{1, 1};
                [stim_e567,poststim1_e567,poststim2_e567] = perform_single_tr_norm(elec567,freqs);
                disp('elec567')

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec24')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                [stim_e8,poststim1_e8,poststim2_e8] = perform_single_tr_norm(elec8,freqs);
                disp('elec8')
                
          
        end
    end 
           stim_e = horzcat(stim_e1, stim_e25, stim_e36, stim_e47, stim_e234, stim_e567, stim_e8); 
       poststim1_e = horzcat(poststim1_e1, poststim1_e25,...
           poststim1_e36, poststim1_e47, poststim1_e234, poststim1_e567, stim_e8); 
       poststim2_e = horzcat(poststim2_e1, poststim2_e25,...
           poststim2_e36, poststim2_e47, poststim2_e234, poststim2_e567, stim_e8);
    
    %%%%%%%%%%%%%%% rSCC 
    
        case 'rSCC' 
       for i = 1:numel(data.output_and_metadata)

        if find(contains(data.output_and_metadata{1, i}{1, 2},'elec25')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                [stim_e1,poststim1_e1,poststim2_e1] = perform_single_tr_norm(elec1,freqs);
                
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec26_27_28')) == 1
                elec234 = data.output_and_metadata{1, i}{1, 1};
                [stim_e234,poststim1_e234,poststim2_e234] = perform_single_tr_norm(elec234,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec26_29')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                [stim_e25,poststim1_e25,poststim2_e25] = perform_single_tr_norm(elec25,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec27_30')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                [stim_e36,poststim1_e36,poststim2_e36] = perform_single_tr_norm(elec36,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec28_31')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                [stim_e47,poststim1_e47,poststim2_e47] = perform_single_tr_norm(elec47,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec29_30_31')) == 1
                elec567 = data.output_and_metadata{1, i}{1, 1};
                [stim_e567,poststim1_e567,poststim2_e567] = perform_single_tr_norm(elec567,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec32')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                [stim_e8,poststim1_e8,poststim2_e8] = perform_single_tr_norm(elec8,freqs);
        end
        
       end % for loop 
       % concatenate 
       stim_e = horzcat(stim_e1, stim_e25, stim_e36, stim_e47, stim_e234, stim_e567, stim_e8); 
       poststim1_e = horzcat(poststim1_e1, poststim1_e25,...
           poststim1_e36, poststim1_e47, poststim1_e234, poststim1_e567, stim_e8); 
       poststim2_e = horzcat(poststim2_e1, poststim2_e25,...
           poststim2_e36, poststim2_e47, poststim2_e234, poststim2_e567, stim_e8);
       
       
        %%%%%%%%%%% lVCVS 
        case 'lVCVS'
        for i = 1:numel(data.output_and_metadata)
         if find(contains(data.output_and_metadata{1, i}{1, 2},'elec1')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                [stim_e1,poststim1_e1,poststim2_e1] = perform_single_tr_norm(elec1,freqs);
                
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec2_3_4')) == 1
                elec234 = data.output_and_metadata{1, i}{1, 1};
                [stim_e234,poststim1_e234,poststim2_e234] = perform_single_tr_norm(elec234,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec2_5')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                [stim_e25,poststim1_e25,poststim2_e25] = perform_single_tr_norm(elec25,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec3_6')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                [stim_e36,poststim1_e36,poststim2_e36] = perform_single_tr_norm(elec36,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec4_7')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                [stim_e47,poststim1_e47,poststim2_e47] = perform_single_tr_norm(elec47,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec5_6_7')) == 1
                elec567 = data.output_and_metadata{1, i}{1, 1};
                [stim_e567,poststim1_e567,poststim2_e567] = perform_single_tr_norm(elec567,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec8')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                [stim_e8,poststim1_e8,poststim2_e8] = perform_single_tr_norm(elec8,freqs);
         end
        end 
           stim_e = horzcat(stim_e1, stim_e25, stim_e36, stim_e47, stim_e8); 
       poststim1_e = horzcat(poststim1_e1, poststim1_e25,...
           poststim1_e36, poststim1_e47, poststim1_e8); 
       poststim2_e = horzcat(poststim2_e1, poststim2_e25,...
           poststim2_e36, poststim2_e47, poststim2_e8);
         %%%%%%%%%%%% rVCVS
    
             
    case 'rVCVS' 
        for i = 1: numel(data.output_and_metadata) 
        if find(contains(data.output_and_metadata{1, i}{1, 2},'elec9')) == 1
                elec1 = data.output_and_metadata{1, i}{1, 1};
                [stim_e1,poststim1_e1,poststim2_e1] = perform_single_tr_norm(elec1,freqs);
                
            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec10_11_12')) == 1
                elec234 = data.output_and_metadata{1, i}{1, 1};
                [stim_e234,poststim1_e234,poststim2_e234] = perform_single_tr_norm(elec234,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec10_13')) == 1
                elec25 = data.output_and_metadata{1, i}{1, 1};
                [stim_e25,poststim1_e25,poststim2_e25] = perform_single_tr_norm(elec25,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec11_14')) == 1
                elec36 = data.output_and_metadata{1, i}{1, 1};
                [stim_e36,poststim1_e36,poststim2_e36] = perform_single_tr_norm(elec36,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec12_15')) == 1
                elec47 = data.output_and_metadata{1, i}{1, 1};
                [stim_e47,poststim1_e47,poststim2_e47] = perform_single_tr_norm(elec47,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec13_14_15')) == 1
                elec567 = data.output_and_metadata{1, i}{1, 1};
                [stim_e567,poststim1_e567,poststim2_e567] = perform_single_tr_norm(elec567,freqs);

            elseif find(contains(data.output_and_metadata{1, i}{1, 2},'elec16')) == 1
                elec8 = data.output_and_metadata{1, i}{1, 1};
                [stim_e8,poststim1_e8,poststim2_e8] = perform_single_tr_norm(elec8,freqs);
         end
        end  
          
        % concatenate 
                   stim_e = horzcat(stim_e1, stim_e25, stim_e36, stim_e47, stim_e8); 
       poststim1_e = horzcat(poststim1_e1, poststim1_e25,...
           poststim1_e36, poststim1_e47, poststim1_e8); 
       poststim2_e = horzcat(poststim2_e1, poststim2_e25,...
           poststim2_e36, poststim2_e47, poststim2_e8);
            
         
end % end for switch case statement 

end 