function [data, idx] = removeStochStimArtifact(dataIn, sync, fs,  win)
% datafile - row vector, data
% sync     - row vector, Cerestim sync channel
% fs       - sampling frequency
% win      - blank window in samples


if size(dataIn, 1)>1
    val = mean(dataIn, 2);
    [~, ch] = max(val);
    datafile = dataIn(ch, :);
else
    datafile = dataIn;
end
datafile(datafile>max(datafile)*0.8) = max(datafile);


sync(sync>(max(sync)*0.8)) = 5e3;
idx =  find(edge(sync));
onset =  idx(1:2:end);
offset = idx(2:2:end);

% Find stim onset indices for trials with regular IPI (F>16Hz)
idx = [];
for n = 1:numel(onset)
    if (offset(n) - onset(n))/fs>0.1
        tmpTrace = datafile(onset(n):offset(n));
        tmpTrace(tmpTrace>max(tmpTrace)*0.7) = max(tmpTrace);
        [~, loc] = findpeaks(tmpTrace, 'MinPeakHeight', max(tmpTrace)*0.7, 'MinPeakDistance', fs*1/100);
        loc = loc - loc(1) + onset(n);
        idx = [idx loc];
    else
        idx = [idx onset(n)];
    end
end

data = dataIn;
if size(data, 1)>1
    for n = 1:numel(idx)
        data(:, idx(n)-3:idx(n)+win) = nan;
    end
    data = fillmissing(data, 'pchip', 2);
else
    for n = 1:numel(idx)
        data(idx(n)-3:idx(n)+win) = nan;
    end
    data = fillmissing(data, 'pchip');
end

        


