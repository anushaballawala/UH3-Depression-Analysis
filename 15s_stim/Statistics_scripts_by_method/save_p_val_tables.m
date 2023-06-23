clear all 
clc

PatientID = 'DBSTRD002'; 

dbs = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/SCC_vs_VCVS_multcompare/%s_mult_compare_data.mat',...
    PatientID)); 

offon_lSCC = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/off_vs_on_multcompare/%s_lSCC_offvson_pvals.mat',...
    PatientID));

offon_rSCC = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/off_vs_on_multcompare/%s_rSCC_offvson_pvals.mat',...
    PatientID));

offon_lVCVS = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/off_vs_on_multcompare/%s_lVCVS_offvson_pvals.mat',...
    PatientID));

offon_rVCVS = load(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/stats/PSD/stats_data/off_vs_on_multcompare/%s_rVCVS_offvson_pvals.mat',...
    PatientID));

dbs.ROI_labels1
