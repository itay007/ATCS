function clearProteinPlotDataStatistics(figHandle)
%CLEARPROTEINPLOTDATASTATISTICS helper function to clear out the Data
% Statistics properties associated with the proteinplot figure. 


%   Copyright 2008 The MathWorks, Inc.

% This will be replaced at some point in the future with a cleaner method
% from the Data Statistics team.

% delete data stats if it is there
% % if ~isempty(findprop(handle(figHandle),'Data_Stats_GUI_Object'))
if ~isempty(bioinfoprivate.biofindprop(figHandle,'Data_Stats_GUI_Object'))
    ds = get(handle(figHandle),'Data_Stats_GUI_Object');
    if ~isempty(ds)
        ds.closeDataStats;
    end
end

% clear appdata for the figure
basicfitdatastat('bfitclearappdata', figHandle);

% clear appdata and Data Stats Properties for all the axes
axesList = findall(figHandle, 'type', 'axes');
for i = 1:length(axesList)
    basicfitdatastat('bfitclearappdata', axesList(i));
    tempProp = bioinfoprivate.biofindprop(axesList(i), 'bfit_AxesListeners');
    if ~isempty(tempProp)
        delete(tempProp)
    end
end

% remove properties (is there a way to get all the properties?
tempProp = bioinfoprivate.biofindprop(figHandle, 'Data_Stats_GUI_Object');
if ~isempty(tempProp)
    delete(tempProp)
end

tempProp = bioinfoprivate.biofindprop(figHandle, 'Data_Stats_Fig_Tag');
if ~isempty(tempProp)
    delete(tempProp)
end

tempProp = bioinfoprivate.biofindprop(figHandle, 'bfit_FigureListeners');
if ~isempty(tempProp)
    delete(tempProp)
end

tempProp = bioinfoprivate.biofindprop(figHandle, 'Basic_Fit_Data_Counter');
if ~isempty(tempProp)
    delete(tempProp)
end

tempProp = bioinfoprivate.biofindprop(figHandle, 'bfit_doublebuffer');
if ~isempty(tempProp)
    delete(tempProp)
end
