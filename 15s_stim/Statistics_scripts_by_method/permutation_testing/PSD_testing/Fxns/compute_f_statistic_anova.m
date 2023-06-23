%% compute_f_statistic_anova.m 

%%% PURPOSE: This fxn runs anova for each ROI/ch and FOI and obtains the
%%% p-value, fstat and degrees of freedom from anova for each ROI, FOI and
%%% is used to obtain these statistic values for permutation testing. 

function [p_value] = compute_f_statistic_anova(data,shuffled_labels)

    for i = 1:size(data,2) % number of channels or ROI
        for j = 1:size(data,3) % number of frequency bands 
            tbl_anova = {}; 
            stats_anova = []; 
            [p(i,j),tbl_anova,stats_anova] = anova1(data(:,i,j),shuffled_labels,'off');
            close all hidden 
            f_stat(i,j) = tbl_anova{2, 5}; 
            p_value(i,j) = tbl_anova{2, 6}  ; 
%             df(i,j) = tbl_anova{2, 3} ;  
        end
    end
    
end
