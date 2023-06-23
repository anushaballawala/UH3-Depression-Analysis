% Take a 128x128x6 matrix and produce an 8x8x6 matrix where rows and cols
% have been averaged over their respective brain regions.
n_regions = 8;
n_channels = 128;
n_freqs = 6;
region_projector = zeros(n_regions, n_channels);

roi_vecs = cell(n_regions, 1);
roi_vecs{1} = ROI_labels.MTG_idx_ch;
roi_vecs{2} = ROI_labels.ACC_idx_ch;
roi_vecs{3} = ROI_labels.AMY_idx_ch;
roi_vecs{4} = ROI_labels.DPF_idx_ch;
roi_vecs{5} = ROI_labels.LOF_idx_ch;
roi_vecs{6} = ROI_labels.VPF_idx_ch;
roi_vecs{7} = ROI_labels.SFG_idx_ch;
roi_vecs{8} = ROI_labels.MOF_idx_ch;

for i = 1:n_regions
    v = roi_vecs{i};
    region_projector(i, v) = 1 / length(v);
end


data_rVCVS_diff_region = zeros(n_regions, n_regions, n_freqs);
for i = 1:n_freqs
    data_rVCVS_diff_region(:, :, i) = region_projector * squeeze(data_rVCVS_diff(:,:,i)) * region_projector';
end

for i = 1:n_freqs
    figure()
    title(sprintf('Freq band %d', i));
    heatmap(squeeze(data_rVCVS_diff_region(:,:,i)))
end