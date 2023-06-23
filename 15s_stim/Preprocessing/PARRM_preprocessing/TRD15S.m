trialDuration=15;
fs=2000;
stimFreq=130;
guessPer=fs/stimFreq;
buffer=round(guessPer/2);


Filtered=highpass(NS3.Data(1,:),1,2000);

for i=1:length(trialStartTimes)
    temp=Filtered(1,round(trialStartTimes(i)*fs)-fs:round(trialStartTimes(i)*fs)+trialDuration*fs-1+fs);
    [pks,locs]=findpeaks(-temp,'MinPeakDistance',0.9*guessPer);
    idx=kmeans(pks',2);
    m1=abs(mean(pks(idx==1)));
    m2=abs(mean(pks(idx==2)));
    if m1>m2
        locs=locs(idx==1);
    else
        locs=locs(idx==2);
    end
    temp=temp(1,locs(1)-buffer:locs(end)+buffer);
    if mod(i,5)==1
        % add in scrollbar to visualize which window is best 
        Period=FindPeriodLFP(temp,[1,length(temp)-1],2000/130);
        PARRM=PeriodicFilter(Period,4000,0.01,0);
    end
    baseline=movmean(temp,4000);
    temp=temp-baseline;
    start=round(trialStartTimes(i)*fs)-fs+locs(1)-buffer-1;
    Filtered(start:start+length(temp)-1)=((filter2(PARRM.',temp','same')-temp')./(1-filter2(PARRM.',ones(size(temp')),'same'))+temp')'+baseline;
end

Filtered=lowpass(Filtered,100,2000);