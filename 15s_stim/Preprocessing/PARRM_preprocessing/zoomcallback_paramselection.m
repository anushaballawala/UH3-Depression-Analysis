function zoomcallback_paramselection(obj,data_onech,Period,sld)
    allAxes = findobj(obj.Children,'Type','Axes');
    numClicked = find(gca==allAxes);
    if numClicked==2
        plot_Paramselection(allAxes(2),allAxes(1),data_onech,Period,sld)
    end
end