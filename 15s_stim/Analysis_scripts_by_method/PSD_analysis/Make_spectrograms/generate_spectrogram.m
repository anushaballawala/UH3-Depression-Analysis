function [] = generate_spectrogram(num_ch,outputdir,data,...
    ch_labels, t, freqs, lower_lim, upper_lim,freq_range)

for i = 1:num_ch  
      fig = figure(); 
        contourf(t,freqs(freq_range),squeeze(mean(data(:,i,freq_range,:),1)),'linecolor','none'); 
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        set(0,'defaultTextInterpreter','none'); 
        set(gca, 'YScale', 'log')
        caxis([lower_lim upper_lim])
        colormap redbluecmap 
        colorbar('westoutside');
        title(ch_labels(i)); 
    set(gca,'ydir','norm')
    set(gca,'ytick',freq_range,'yticklabel',round(freqs(freq_range)));    
    han=axes(fig,'visible','off');
    han.YLabel.Visible='on';
    ylabel(han,'log of Power in freq band');
    filename = sprintf('%s/%s.png',outputdir,ch_labels(i));
    saveas(gcf, filename); 
end 
end 