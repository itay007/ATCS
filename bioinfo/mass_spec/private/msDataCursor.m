function msDataCursor(hFig)
% MSDATACURSOR sets up the Data Cursor for mass spectrometry figures.

%   Copyright 2007 The MathWorks, Inc.


hDCM = datacursormode(hFig); % Does initial creation for Data Cursor Mode
set(hDCM,'UpdateFcn',@msDataCursorUpdateFcn,'Enable','on'); % Creates the context menu
hMenu = get(hDCM,'UIContextMenu'); % Get its context menu.
delete(findobj(hMenu,'Tag','DataCursorEditText'))
delete(findobj(hMenu,'Tag','DataCursorSelectText'));
set(hDCM,'Enable','off') %starts by default off

function txt = msDataCursorUpdateFcn(obj,event_obj) %#ok
% MSDATACURSORUPDATEFCN Function to update the Data Cursor text in all the
% mass spectrometry static customized figures. 


% Display 'M/Z' and 'Intensity'
pos = get(event_obj,'Position');
target = get(event_obj,'Target');
switch get(target,'Type')
    case 'image'
        switch getappdata(gcbf,'msheatmapType')
            case 'lcms'
               T = getappdata(gcbf,'msheatmapRetentionTimeVector');
               txt = {['M/Z: ',num2str(pos(1))],...
                    ['Time: ',num2str(T(pos(2)))]};
            case 'categorical'
                hIDs = getappdata(gcbf,'msheatmapCategoricalOriginalIDs');
                txt = {['M/Z: ',num2str(pos(1))],...
                    ['Spec ID: ',num2str(hIDs(pos(2)))]};
            case 'indexed'
                txt = {['M/Z: ',num2str(pos(1))],...
                    ['Spec ID: ',num2str(pos(2))]};
            otherwise
                txt = {'Datacursor unavailable'};
        end
    case 'line'
        switch get(target,'Tag')
            case 'msdotplot'
                txt = {['M/Z: ',num2str(pos(1))],...
                       ['Time: ',num2str(pos(2))],...
                       ['Intensity: ',num2str(pos(3))]};
            case 'mzmarker'
                txt = {['M/Z: ',num2str(pos(1))]};
            case 'mzcounts'
                txt = {['M/Z: ',num2str(pos(1))],...
                       ['Counts: ',num2str(pos(2))]};
            case 'mzdistance'
                txt = {['M/Z: ',num2str(pos(1))],...
                       ['Distance: ',num2str(pos(2))]}; 
            case 'mschroma'
                txt = {['Time: ',num2str(pos(1))],...
                       ['Intensity: ',num2str(pos(2))]};
            case 'spectrum'
                txt = {['M/Z: ',num2str(pos(1))],...
                       ['Intensity: ',num2str(pos(2))]};
            otherwise % we assume that any other line type has coordinates 
                      % 'M/Z' vs 'Intensity' 
                txt = {['M/Z: ',num2str(pos(1))],...
                       ['Intensity: ',num2str(pos(2))]};
        end
    otherwise % unrecognized hg object
        txt = [];
end
