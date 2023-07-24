function [timestamps,cerestim,stimsync] = extract_timestamps_analog(PatientID,stim_info)

%% get timestamps from analog channel 


%load .ns3 file 
%fprintf('Select NS5 files %s',NEVfilename)
[filename,filepath]=uigetfile('*.ns5','Select NS5 files');
NS5 = openNSx(fullfile(filepath,filename));
%% plot comments and analog channel together 
fs_orig = 30000;
fs_new = 1000;
sync_ch = 185; 
cerestim_ch = 186; 
%downsample data
stimsync = decimate(double(NS5.Data(sync_ch,:)),fs_orig / fs_new);
cerestim = decimate(double(NS5.Data(cerestim_ch,:)),fs_orig / fs_new);

%%
y = stimsync;
x = 1/fs_new * ((1:length(stimsync)) - 1);
plot(x,y)
N = 10;
stimsync_filt = medfilt1(stimsync,N); 

maxi = max(stimsync); 
minpeakh = maxi - 5000; 
minpeakd = 15; %start of stim should be at least 15 seconds apart 

[pks,locs] = findpeaks(stimsync_filt,1000,'MinPeakDistance',minpeakd,'MinPeakHeight',minpeakh);

findpeaks(stimsync_filt,1000,'MinPeakDistance',minpeakd,'MinPeakHeight',minpeakh);

plot(x,stimsync_filt)
hold on 
scatter(locs,pks)
figdir = sprintf('/gpfs/data/dborton/TRD_Project/DBSTRD/%s/EXP/15s_stim/Figures/FindPeaks_Timestamps_%s.fig',PatientID,stim_info); 
saveas(gcf, figdir)
timestamps = locs; 

%[~,locs,~,~] = findpeaks(cerestim,fs_new,'MinPeakDistance',0.005,'MinPeakHeight',3000);



end 