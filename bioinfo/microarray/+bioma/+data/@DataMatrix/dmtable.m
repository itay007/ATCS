function varargout = dmtable(obj, varargin)
%DMTABLE Create a two dimensional graphic table showing DataMatrix.
%
%   H = DMTABLE(DM) creates a two dimensional graphic table with data of a
%   DataMatrix DM and returns the uitable handle.
%
%   DMTABLE(DM, FH) specifies the parent handle of the uitable. The parent
%   can be a figure or uipanel handle.
%
%   Example:
%
%   % Create a DataMatrix d 
%   d = bioma.data.DataMatrix(rand(2,3), {'Row1', 'Row2'}, {'Col1', 'Col2', 'Col3'})
%   ht = dmtable(d)
%
%   See also DATAMATRIX.

%   Copyright 2009-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename)

%== Check input arguments
if nargin > 1
    arg = varargin{1};
    if ishandle(arg)
        hParent = arg;
    end
else
    figureName = bioinfoprivate.indexedFigureName('DataMatrix', 'dmtable');
    hParent = figure('Renderer',' zbuffer',...
        'Name', figureName,...
        'NumberTitle','off',...
        'Tag', 'DataMatrix',...
        'IntegerHandle','off',...
        'PaperPositionMode', 'auto',...
        'Visible', 'on');
end       

hTable = uitable('Parent', hParent);

%== Show data
set(hTable, 'Data', obj.Matrix)
set(hTable, 'ColumnName', obj.ColNames,...
            'RowName', obj.RowNames);
appdata = getAppData(hParent);
set(hParent, 'ResizeFcn', {@bioma.util.parentResizeCB,...
                           hTable, appdata.widthRatio, appdata.heightRatio});
setAppData(hParent,appdata); 
%== Output
if nargout >  0
    varargout{1} = hTable;
end
end % end of dmtable

% --------------
function setAppData(hfig,appdata)
setappdata(hfig,'DataMatrix',appdata);
end

function [appdata] = getAppData(hfig)

if isappdata(hfig,'DataMatrix')
    appdata = getappdata(hfig,'DataMatrix');
else
    appdata = guihandles(hfig);
    appdata.widthRatio = 1/20;
    appdata.heightRatio = 1/50;
end
end


