function [mdl,slope,r_squared]= get_1overf_slope(f,Pxx)

    mdl = fitlm(f,Pxx);
    slope = mdl.Coefficients.Estimate(2);
    r_squared = mdl.Rsquared.Ordinary; 
    
    
end 
