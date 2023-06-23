% Plot the data across time for each channel.
time_axis = 0 : (1/fs) : (num_samples/fs); 
time_axis = time_axis(1:end-1);

parfor ch = 1:num_channels
    figure('Position',[200, 200, 2200, 600]);
    hold on;

    plot(time_axis, mu(ch, :), 'color', [0.6353, 0.0784, 0.1843]);
    xline(5,'linewidth', 2.5); 
    xline(20, 'linewidth', 2.5); 
    ylabel('Normalized freqband power (dB)');
    xlabel('Time[s]');
    title(sprintf("Channel %s", channel_label(ch)));
    
    plot(time_axis, CI_upper(ch, :));
    plot(time_axis, CI_lower(ch, :));
    
    patch(...
        [time_axis, fliplr(time_axis)], ...
        [CI_upper(ch, :), fliplr(CI_lower(ch, :))], ...
        [ 1, 0, 0], ...
        'FaceAlpha', 0.3);
    name = sprintf('%s_freqbandpowerovertime.png',deblank(channel_label(ch)));
    saveas(gcf,name); 
end