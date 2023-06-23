
%tbl = py_table((py_table.ROI == 'ACC'|py_table.ROI == 'aCC') & py_table.Target == 'lSCC',:); 

function [tbl] = Jan16_make_table_py_viz(ROI_name,DBS_target,py_table,PatientID)
        
    crit_idx = (py_table.ROI == ROI_name & py_table.Target == DBS_target) | (py_table.ROI == ROI_name & py_table.Target == 'BaselineFix'); 
    tbl = py_table(crit_idx,:); 
    addpath(genpath('/Users/anushaallawala/Python_projects/Stim_paper/01-16-23-paper/'))
    writetable(tbl,sprintf('/Users/anushaallawala/Python_projects/Stim_paper/01-16-23-paper/%s_%s_%s_offvson.csv',PatientID,ROI_name,DBS_target));
end 

