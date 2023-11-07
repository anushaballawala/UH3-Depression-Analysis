
fname = 'subj-DBSTRD006_task-Stim15s_run-lSCC_blk-f130_subblk-1.ns3';
NS3 = openNSx(fname, 'uv');

elecLabels = cat(1, NS3.ElectrodesInfo.Label);
elecLabels = mat2cell(elecLabels, ones(size(NS3.Data, 1), 1), 16);
elecLabels = strrep(elecLabels, elecLabels{1}(end), '');

chList = find(contains(elecLabels, 'Amy'));
chSel = chList(2);
data = NS3.Data;
trace = data(chSel, :);

fs = NS3.MetaTags.SamplingFreq;
stimFreq = 130;
[pks, loc] = findpeaks(data(65, 1:6.669e5), 'MinPeakHeight', 1e3, 'MinPeakDistance', fs*1/(2*stimFreq));


lpFilt = designfilt('lowpassfir','PassbandFrequency', 80, ...
'StopbandFrequency', 100, 'PassbandRipple', 0.5, ...
'StopbandAttenuation', 40, 'SampleRate', fs);

filtdata_lp80Hz = filtfilt(lpFilt, trace);

snipSamples = ceil(fs/130);
preSamples = 6;
postSamples = snipSamples - preSamples;
snips = [];
for n = 1:numel(loc)
    snips(n, :) = data(65, loc(n)-6:loc(n+11);
end

idx = loc;
[~, ix] = min(diff(snips', 2));
idx = idx - (ix(1) - ix');
snips = [];
for n = 1:numel(loc)
    snips(n, :) = data(65, idx(n)-preSamples:idx(n)+postSamples);
end
g = kmeans(snips, 4);
idx = idx + 1;


for k = 1:4
    vals = idx(g==k);
    for n = 1:numel(vals)
        trace(vals(n)-preSamples:vals(n)+postSamples) = data(65, vals(n)-preSamples:vals(n)+postSamples) - mean(snips(g==k, :));
    end
end


lpFilt = designfilt('lowpassfir','PassbandFrequency', 150, ...
    'StopbandFrequency', 200, 'PassbandRipple', 0.5, ...
    'StopbandAttenuation', 40, 'SampleRate', fs);

filtdata_lp150Hz = filtfilt(trace);

