function multiselect
fig = uifigure('Position',[100 100 350 275]);

% Create Text Area
txt = uitextarea(fig,...
    'Position',[125 80 100 50]);

% Create List Box
lbox = uilistbox(fig,...
    'Position',[125 150 100 78],...
    'Multiselect','on',...
    'ValueChangedFcn',@selectionChanged);
lbox.Multiselect = 'on';

% ValueChangedFcn callback
function selectionChanged(src,event)
    txt.Value = src.Value;
    items = lbox.Value
end


end



