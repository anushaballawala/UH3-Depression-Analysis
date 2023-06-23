%% 
clear 
stim_info_input = 'VCVSf130_234'; 
neuraldatafile = dir(sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD003/Raw Data/15s stim/Stim15s_%s/*.mat',stim_info_input)); 
neuraldata = load(neuraldatafile.name); 

stim_info_output = 'VCVS_f130_234'; 
outputdir = sprintf('/users/aallawa1/data/TRD_Project/DBSTRD/DBSTRD003/Experiments/15s_stim/Preprocessed Data/%s/%s_goodch.mat',stim_info_output,stim_info_output); 

%% load montage info 

montageInfo = getMontageInfo(neuraldata); 

%% 

good_ch_idx = [2:5,8,11:14,18:19,22,24:26,30:39,41:43,46:47,50:52,...
    102:117,120,122:125,128:130,133:135,136:141]; 

good_ch_labels = montageInfo.ElectrodeLabels(good_ch_idx);

%% Save data 

save(outputdir, 'good_ch_idx','good_ch_labels'); 



