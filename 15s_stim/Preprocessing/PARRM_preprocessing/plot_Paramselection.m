
function [sld] = plot_Paramselection(ax1,ax2,data_onech,Period,sld)
figure
ax1 = subplot(1,2,1);
plot(mod(1:initWinSize*2-1,Period),diff(data_onech(1:initWinSize*2)),'o')
axis tight
ax2 = subplot(1,2,2);
plot(data_onech)
axis tight
h2 = uicontrol('style','slider','units','pixel','position',[20 20 300 20],'min',1,'max',round((length(data_onech)-1)/2),...
'callback',@(hObject, event) reset(ax1,ax2,data_onech,Period,hObject),'val',initWinSize);
h=zoom;
h.ActionPostCallback = @(obj,evd) zoomcallback(obj,data_onech,Period,h2);
h.Enable = 'on';
h3=uicontrol('style','text','position',[20 40 300 20],'String',"Window size ");