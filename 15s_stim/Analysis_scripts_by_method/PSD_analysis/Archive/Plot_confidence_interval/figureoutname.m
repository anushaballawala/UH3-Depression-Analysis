% load indiv files 

%load one contact configuration for each target (loading contact 8)
lVCVS_file = load('15s_stim_elec16_lVCVS_f130.mat'); 
rVCVS_file = load('15s_stim_elec16_rVCVS_f130.mat'); 
lSCC_file = load('15s_stim_cond7_lSCC_f130.mat'); 
rSCC_file = load('15s_stim_cond7_rSCC_f130.mat'); 

%% 
lVCVS = lVCVS_file.indiv2; 
rVCVS = rVCVS_file.indiv2; 
lSCC = lSCC_file.indiv7; 
rSCC = rSCC_file.indiv7; 

pre_stim = 1:5000; 
stim_win = 5001:20000; 
post_stim_win1 = 20001:25000; 
post_stim_win2 = 25001:30000; 
post_stim_total_win = 20001:30000; 

%% normalize data to dB 

theta = 4:7; 
%get average baseline of 5 trials 
bl_win_pw_lVCVS = mean(mean(lVCVS(:,:,:,pre_stim),4),1);
bl_win_pw_rVCVS = mean(mean(rVCVS(:,:,:,pre_stim),4),1);
bl_win_pw_lSCC = mean(mean(lSCC(:,:,:,pre_stim),4),1);
bl_win_pw_rSCC = mean(mean(rSCC(:,:,:,pre_stim),4),1);

%normalize to dB for each trial in each channel 

parfor i = 1:5
    norm_pw_lVCVS(i,:,:,:) = 10*log10( bsxfun(@rdivide,lVCVS(i,:,:,:),bl_win_pw_lVCVS)); %divide power by baseline 
    norm_pw_rVCVS(i,:,:,:) = 10*log10( bsxfun(@rdivide,rVCVS(i,:,:,:),bl_win_pw_rVCVS));
    norm_pw_lSCC(i,:,:,:) = 10*log10( bsxfun(@rdivide,lSCC(i,:,:,:),bl_win_pw_lSCC));
    norm_pw_rSCC(i,:,:,:) = 10*log10( bsxfun(@rdivide,rSCC(i,:,:,:),bl_win_pw_rSCC));
end 


%% index theta power for each trial 

% theta_tbl = array2table(squeeze(norm_pw_lVCVS(:,:,theta,pre_stim)),...
%     norm_pw_lVCVS(:,:,theta,stim_win),norm_pw_lVCVS(:,:,theta,post_stim_win1),...
%     norm_pw_lVCVS(:,:,theta,post_stim_win2),norm_pw_rVCVS(:,:,theta,pre_stim),...
%     norm_pw_rVCVS(:,:,theta,stim_win),norm_pw_rVCVS(:,:,theta,post_stim_win1),...
%     norm_pw_rVCVS(:,:,theta,post_stim_win2),norm_pw_lSCC(:,:,theta,pre_stim),...
%     norm_pw_lSCC(:,:,theta,stim_win),norm_pw_lSCC(:,:,theta,post_stim_win1),...
%     norm_pw_lSCC(:,:,theta,post_stim_win2),norm_pw_rSCC(:,:,theta,pre_stim),...
%     norm_pw_rSCC(:,:,theta,stim_win),norm_pw_rSCC(:,:,theta,post_stim_win1),...
%     norm_pw_rSCC(:,:,theta,post_stim_win2),'VariableNames',...
%     {'lVCVS_theta_pre','lVCVS_theta_stim','lVCVS_theta_post1','lVCVS_theta_post2',...
%     'rVCVS_theta_pre','rVCVS_theta_stim','rVCVS_theta_post1','rVCVS_theta_post2',...
%     'lSCC_theta_pre','lSCC_theta_stim','lSCC_theta_post1','lSCC_theta_post2',...
%     'rSCC_theta_pre','rSCC_theta_stim','rSCC_theta_post1','lSCC_theta_post2'}); 

LVCVS_pre_theta = mean(mean(norm_pw_lVCVS(:,:,theta,pre_stim),4),3); 
RVCVS_pre_theta = mean(mean(norm_pw_rVCVS(:,:,theta,pre_stim),4),3); 
LSCC_pre_theta = mean(mean(norm_pw_lSCC(:,:,theta,pre_stim),4),3); 
RSCC_pre_theta = mean(mean(norm_pw_rSCC(:,:,theta,pre_stim),4),3); 

LVCVS_theta = mean(mean(norm_pw_lVCVS(:,:,theta,stim_win),4),3); 
RVCVS_theta = mean(mean(norm_pw_rVCVS(:,:,theta,stim_win),4),3); 
LSCC_theta = mean(mean(norm_pw_lSCC(:,:,theta,stim_win),4),3); 
RSCC_theta = mean(mean(norm_pw_rSCC(:,:,theta,stim_win),4),3); 

LVCVS_theta_alltrials = mean(LVCVS_theta,1); 
RVCVS_theta_alltrials = mean(RVCVS_theta,1); 
LSCC_theta_alltrials = mean(LSCC_theta,1); 
RSCC_theta_alltrials = mean(RSCC_theta,1); 

channellabel = deblank(rSCC_file.metadata.preprocessing.GoodChannelLabels);  
%% 
num_ch = 136; %**
condition = {'LSCC_pre_theta','LSCC_thetastim','RSCC_pre_theta','RSCC_thetastim',...
    'LVCVS_pre_theta','LVCVS_thetastim','RVCVS_pre_theta','RVCVS_thetastim'}; 


%% 

x1 = repmat(1,1,5); 
x2 = repmat(2,1,5); 
x3 = repmat(3,1,5); 
x4= repmat(4,1,5); 
x5= repmat(5,1,5); 
x6= repmat(6,1,5); 
x7= repmat(7,1,5); 
x8= repmat(8,1,5); 


parfor i = 1:136
    figure() 
    scatter(x1,LSCC_pre_theta(:,i),'filled'); 
    hold on 
    scatter(x2,LSCC_theta(:,i),'filled'); 
    scatter(x3,RSCC_pre_theta(:,i),'filled'); 
    scatter(x4,RSCC_theta(:,i),'filled'); 
    scatter(x5,LVCVS_pre_theta(:,i),'filled');
    scatter(x6,LVCVS_theta(:,i),'filled'); 
    scatter(x7,RVCVS_pre_theta(:,i),'filled'); 
    scatter(x8,RVCVS_theta(:,i),'filled'); 
    hold off 
    ax = gca; 
    ax.XTick = 1:1:8; 
    ax.XTickLabel = condition; 
    xtickangle(45)
    set(gca,'TickLabelInterpreter','none'); 
    ylabel('dB')
    title(sprintf('theta power - Ch%s',channellabel(i))); 
    filename = sprintf('Ch%s.png',channellabel(i)); 
    saveas(gcf,filename); 

end 
%% 

% scatter(x,total); 
% 
% parfor i = 1:2 %num_ch
%     total = [LSCC_pre_theta(:,i),LSCC_theta(:,i),RSCC_pre_theta(:,i),RSCC_theta(:,i),...
%         LVCVS_pre_theta(:,i),LVCVS_theta(:,i),RVCVS_pre_theta(:,i),RVCVS_theta(:,i)]; 
%     figure()
%     scatter(x,LSCC_pre_theta(:,1));
%     hold on 
%     scatter(x1,LSCC_theta(:,1)); 
%     ylabel('dB')
%     title(sprintf('theta power - Ch%s',channellabel(i))); 
%     filename = sprintf('Ch%s.png',channellabel(i)); 
%     %saveas(gcf,filename); 
%     ylabel('dB')
% end 
%% plot theta power over time 

%get theta over time for each stim condition 

norm_LSCC_theta = norm_pw_lSCC(:,:,theta,:); 
norm_RSCC_theta = norm_pw_rSCC(:,:,theta,:); 
norm_LVCVS_theta = norm_pw_lVCVS(:,:,theta,:);
norm_RVCVS_theta = norm_pw_rVCVS(:,:,theta,:);

LSCC_theta_avg = squeeze(mean(norm_LSCC_theta,3)); 
RSCC_theta_avg = squeeze(mean(norm_RSCC_theta,3)); 
LVCVS_theta_avg = squeeze(mean(norm_LVCVS_theta,3)); 
RVCVS_theta_avg = squeeze(mean(norm_RVCVS_theta,3)); 

mu_LSCC_theta = squeeze(mean(mean(norm_LSCC_theta,3),1)); 
mu_RSCC_theta = squeeze(mean(mean(norm_LSCC_theta,3),1));
mu_LVCVS_theta = squeeze(mean(mean(norm_LSCC_theta,3),1));
mu_RVCVS_theta  = squeeze(mean(mean(norm_LSCC_theta,3),1));

sd_LSCC_theta = squeeze(std(LSCC_theta_avg,0,1)); 
sdmn_LSCC_theta = sd_LSCC_theta./sqrt(5);
CI_upper = mn+2*sdmn_LSCC_theta; 
CI_lower = mn-2*sdmn_LSCC_theta; 

%% compute for every channel 

norm_LSCC_theta = norm_pw_lSCC(:,:,theta,:); 

LCC_theta_avg = []; 
mu_LSCC = []; 
sd_LSCC = []; 
sdmn_LSCC = []; 
CI_upper = []; 
CI_lower = []; 

for i = 1:136 
    LSCC_theta_avg(:,i,:,:) = mean(norm_LSCC_theta(:,i,:,:),3); 
    mu_LSCC(1,i,:) = mean(LSCC_theta_avg(:,i,:),1); 
    sd_LSCC(1,i,:) = std(LSCC_theta_avg(:,i,:,:),0,1); 
    sdmn_LSCC(1,i,:) = sd_LSCC(1,i,:)./sqrt(5); 
    CI_upper(1,i,:) = mu_LSCC(1,i,:)+2*sdmn_LSCC(1,i,:); 
    CI_lower(1,i,:) = mu_LSCC(1,i,:)-2*sdmn_LSCC(1,i,:); 
end 
            
mu = squeeze(mu_LSCC); 
CI_U = squeeze(CI_upper); 
CI_L = squeeze(CI_lower); 
%% RSCC
RSCC_theta_avg = []; 
mu_RSCC = []; 
sd_RSCC = []; 
sdmn_RSCC = []; 
CI_upper = []; 
CI_lower = []; 

for i = 1:136 
    RSCC_theta_avg(:,i,:,:) = mean(norm_RSCC_theta(:,i,:,:),3); 
    mu_RSCC(1,i,:) = mean(RSCC_theta_avg(:,i,:),1); 
    sd_RSCC(1,i,:) = std(RSCC_theta_avg(:,i,:,:),0,1); 
    sdmn_RSCC(1,i,:) = sd_RSCC(1,i,:)./sqrt(5); 
    CI_upper(1,i,:) = mu_RSCC(1,i,:)+2*sdmn_RSCC(1,i,:); 
    CI_lower(1,i,:) = mu_RSCC(1,i,:)-2*sdmn_RSCC(1,i,:); 
end 
            
mu = squeeze(mu_RSCC); 
CI_U = squeeze(CI_upper); 
CI_L = squeeze(CI_lower);


%% 
% parfor i = 1:136 
% figure('Position',[200 200 2200 600])
% plot(timeAx,mu(i,x),'color',[0.6353    0.0784    0.1843])
% hold on 
% xline(5,'linewidth',2.5); 
% xline(20,'linewidth',2.5); 
% ylabel('Normalized theta power')
% xlabel('Time[s]'); 
% plot(timeAx,CI_U(i,x))
% plot(timeAx,CI_L(i,x))
% patch([timeAx,fliplr(timeAx)], [CI_U(i,:), fliplr(CI_L(i,:))],[ 1     0     0], 'FaceAlpha',0.3)
% filename = (sprintf('RSCC_theta_%s.png',ch_labels(i))); 
% saveas(gcf,filename); 



%end 


%% 
x = 1:30000; 
    srate_new = 1000; fs = 1000;
    dataDuration = 30000;
timeAx = [ 0: (1/fs):dataDuration/srate_new ]; 
timeAx = timeAx(1:end-1);

%% 
ch_labels = deblank(rSCC_file.metadata.preprocessing.GoodChannelLabels);   

% parfor i = 1:136 
% figure('Position',[200 200 2200 600])
% plot(timeAx,mu(i,x),'color',[0.6353    0.0784    0.1843])
% hold on 
% xline(5,'linewidth',2.5); 
% xline(20,'linewidth',2.5); 
% ylabel('Normalized theta power')
% xlabel('Time[s]'); 
% plot(timeAx,CI_U(i,x))
% plot(timeAx,CI_L(i,x))
% patch([timeAx,fliplr(timeAx)], [CI_U(i,:), fliplr(CI_L(i,:))],[ 1     0     0], 'FaceAlpha',0.3)
% filename = (sprintf('LSCC_theta_%s.png',ch_labels(i))); 
% saveas(gcf,filename); 



%end 


%% 

















