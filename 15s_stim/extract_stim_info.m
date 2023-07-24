function [PatientID,PatientType,block_file,block_num,...
    freq,hemi, DBStarget] = extract_stim_info(fullfile) 

PatientID = getPatientID(fullfile); 
PatientType = getPatientType(fullfile); 

%% get subblock 
if contains(fullfile,'1-1')
    subblk = '1-1';
    
elseif contains(fullfile,'e1_repeat')
    subblk = 'e1_repeat';

else
    subblk = 'subblk-\d*';
    [start_idx, end_idx] = regexp(fullfile,  subblk);
    subblk = extractBetween(fullfile, start_idx, end_idx);
    subblk = subblk{1};
    subblk = subblk(8);
end

if contains(fullfile,'e1_repeat')
    block_file = subblk; 
    block_num = subblk; 
else 
    block_file = subblk; 
    block_num = sprintf('0%s',block_file); 
end 

    %% get freq 
    
    freq_txt = 'f\d*';
    [start_idx, end_idx] = regexp(fullfile,  freq_txt);
    start_idx = start_idx(1); 
    end_idx = end_idx(1);
    freq = extractBetween(fullfile, start_idx, end_idx);
    freq = freq{1}; 
    freq = freq(2:end); 
    
%% get DBStarget 

braintarget = '(VCVS|SCC)';
[start_idx_br, end_idx_br] = regexp(fullfile,braintarget);
start_idx_br = start_idx_br(1); 
end_idx_br = end_idx_br(1); 
DBStarget_full = extractBetween(fullfile,start_idx_br-1,end_idx_br);
DBStarget_full = DBStarget_full{1}; 
hemi = DBStarget_full(1); 

DBStarget = DBStarget_full(2:end); 

end 




