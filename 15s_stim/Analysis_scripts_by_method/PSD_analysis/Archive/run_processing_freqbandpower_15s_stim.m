
file = load('15s_stim_elec16_lVCVS_f130.mat');
filename = 'lVCVS_f130_ele8'; 
freqs = 3:7; 
[tbl] = processing_freqbandpower_15s_stim(file,filename,freqs); 


filename_csv = 'E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/050521/lVCVS_f130_elec8_theta.csv';

writetable(tbl,filename_csv);

%% 
clear all 
clc 
file = load('15s_stim_elec16_rVCVS_f130.mat');
filename = 'rVCVS_f130_ele8'; 
[tbl] = processing_freqbandpower_15s_stim(file,filename); 

filename_csv = 'E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/050521/rVCVS_f130_elec8_theta.csv';
writetable(tbl,filename_csv);

%% 
file = load('15s_stim_cond7_lSCC_f130.mat');
filename = 'lSCC_f130_ele8'; 
[tbl] = processing_freqbandpower_15s_stim(file,filename); 

filename_csv = 'E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/050521/lSCC_f130_elec8_theta.csv';
writetable(tbl,filename_csv);

%% 
file = load('15s_stim_cond7_rSCC_f130.mat');
filename = 'rSCC_f130_ele8'; 
[tbl] = processing_freqbandpower_15s_stim(file,filename); 

filename_csv = 'E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Processed Data/050521/rSCC_f130_elec8_theta.csv';
writetable(tbl,filename_csv);
