%% Plot all normalized values across time 

% set up directory 

files = dir('E:/DBSTRD/DBSTRD001/Experiments/15s_stim/Epoched Data/*/*trialaveraged*.mat');
[datafile, data] = load_files_from_dir(files); 

% define baseline period 
pre_stim_win = 1:4800; %1:5000
theta = 4:7; 

% create baseline variables 
lSCC = data{1, 1}.trial_averaged_pw{7, 1} ; 
rSCC = data{3, 1}.trial_averaged_pw{7, 1} ; 
lVCVS = data{2, 1}.trial_averaged_pw{2, 1} ; 
rVCVS = data{4, 1}.trial_averaged_pw{2, 1} ; 

% compute baseline average for each target in theta band 

lSCC_theta = lSCC(:,theta,pre_stim_win); 
lSCC_avg = mean(lSCC_theta, 2); 
lSCC_avg_avg = mean(lSCC_avg, 3); 

rSCC_theta = rSCC(:,theta,pre_stim_win);
rSCC_avg = mean(rSCC_theta,2); 
rSCC_avg_avg = mean(rSCC_avg, 3); 

lVCVS_theta = lVCVS(:,theta,pre_stim_win);
lVCVS_avg = mean(lVCVS_theta,2);
lVCVS_avg_avg = mean(lVCVS_avg,3); 

rVCVS_theta = rVCVS(:,theta,pre_stim_win);
rVCVS_avg = mean(rVCVS_theta,2); 
rVCVS_avg_avg = mean(rVCVS_avg,3); 
%% 
% plot log power 
figure()
scatter(1:136,10*log10(lSCC_avg_avg),'filled'); 
hold on 
scatter(1:136,10*log10(rSCC_avg_avg),'filled'); 
scatter(1:136,10*log10(lVCVS_avg_avg),'filled'); 
scatter(1:136,10*log10(rVCVS_avg_avg),'filled'); 

ylabel('Power in dB')
xlabel('Channel number') 
legend('lSCC','rSCC','lVCVS','rVCVS')

% plot raw power 
figure() 
scatter(1:136,lSCC_avg_avg,'filled'); 
hold on 
scatter(1:136,rSCC_avg_avg,'filled'); 
scatter(1:136,lVCVS_avg_avg,'filled'); 
scatter(1:136,rVCVS_avg_avg,'filled'); 
xlabel('Channel number') 
ylabel('Raw power in uV^2/Hz') 
legend('lSCC','rSCC','lVCVS','rVCVS')
ylim([0 10^6])

