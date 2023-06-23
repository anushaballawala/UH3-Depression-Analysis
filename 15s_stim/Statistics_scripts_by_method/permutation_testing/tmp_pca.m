
% trials X brain region x frequency. 
% first 35 trials are SCC 
% last 24 trials are VCVS 


matrix_ROI_left;

lVCVS_mean = squeeze(mean(lVCVS_data, 1));
lSCC_mean = squeeze(mean(lSCC_data, 1));



close all;
figure(42);

freq_label = 1:6;
freq_label = reshape(freq_label, 1, 1, length(freq_label));
freq_label = repmat(freq_label, 8, 1);

scatter(lSCC_mean(:), lVCVS_mean(:), [], freq_label(:), 'filled');

%%
freq_idx = 1;

Y = squeeze(mean(lVCVS_data(:,:,freq_idx), 1));
X = squeeze(mean(lSCC_data(:,:,freq_idx), 1));

close all;
figure(42);

roi_label = 1:8;
scatter(X(:), Y(:), [], roi_label(:), 'filled');


%%

 % trials X brain region x frequency. 
 freq_idx = 3;

v = squeeze(lVCVS_data(:,:,freq_idx));
s = squeeze(lSCC_data(:,:,freq_idx));

% trials x brain region

data = vertcat(v, s);
 [coeff, score, latent] = pca(data);
 
 score_v = score(1:size(v, 1), :);
 score_s = score(size(v, 1) + 1:end, :);
 
 close all;
 figure (45);
 scatter(score_v(:, 1), score_v(:,2), [], 'r', 'filled');
 hold on;
 scatter(score_s(:, 1), score_s(:,2), [], 'b', 'filled');