function[pre_stim_win, stim_win, post_stim_win, stim_win1, stim_win2, stim_win3] = define_window_of_interest(PatientID, prestim, stim, poststim)

% Define the different time windows of interest, in samples.
if strcmp(prestim,'full window') == 1 && strcmp(stim, 'full window') == 1 && strcmp(poststim,'full window') ==1 
    disp('Using 5seconds post stim') 
     switch PatientID 
    case 'DBSTRD001' 
        pre_stim_win = 1:5000; 
        stim_win = 5010:20000; 
        post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
        post_stim_win2 = 25001:30000; 
        %post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
        post_stim_win = 21000:26000; 
    case 'DBSTRD002'     
        pre_stim_win = 1:4500; 
        stim_win = 5001:20000; 
        post_stim_win1 = 21000:25000; % stimulation ends a bit after the 20 second mark 
        post_stim_win2 = 25001:30000; 
        %post_stim_total_win = horzcat(post_stim_win1, post_stim_win2);
        post_stim_win = 21000:26000; 
    end 
else
    disp('splitting stim window into 5s windows') 
    
    
    switch PatientID 
        case 'DBSTRD001'
            pre_stim_win = 1:4900; 
            stim_win = 5000:2000;
            stim_win1 = 5010:10000;
            stim_win2 = 10000:15000;
            stim_win3 = 15000:19990;
            post_stim_win = 21000:26000; 
            
        case 'DBSTRD002'
            pre_stim_win = 1:4900; 
            stim_win = 5000:2000;
            stim_win1 = 5010:10000;
            stim_win2 = 10000:15000;
            stim_win3 = 15000:19990;
            post_stim_win = 21000:26000; 
    end 
    
    
end 

end 