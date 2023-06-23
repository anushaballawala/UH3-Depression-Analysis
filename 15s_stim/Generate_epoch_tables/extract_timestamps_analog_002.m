function [timestamps,cerestim,stimsync] = extract_timestamps_analog(subjectname,stim_info)

%% get timestamps from analog channel 


%load .ns3 file 
disp('Select NS5 files for 15s STIM')
[filename,filepath]=uigetfile('*.ns5','Select NS3 File');
NS5 = openNSx(fullfile(filepath,filename));
%% plot comments and analog channel together 
fs_orig = 30000;
fs_new = 1000;
%downsample data 
stimsync = decimate(double(NS5.Data(3,:)),fs_orig / fs_new); 
cerestim = decimate(double(NS5.Data(4,:)),fs_orig / fs_new); 

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
figdir = sprintf('E/DBSTRD/%s/Experiments/15s_stim/Figures/AA/FindPeaks_Timestamps/%s.png',subjectname,stim_info); 
saveas(gcf, figdir)
timestamps = locs; 

%[~,locs,~,~] = findpeaks(cerestim,fs_new,'MinPeakDistance',0.005,'MinPeakHeight',3000);



end 