classdef HeatMap < hgsetget
    %HEATMAP A false color 2D image of the data values in a matrix.
    %   HeatMap displays the matrix data values as colors in a
    %   two-dimensional map. 
    %
    %   HeatMap properties:
    %       Standardize           - Dimension in which the values are standardized.
    %       Colormap              - Colormap used to display the heat map.
    %       DisplayRange          - The range of the values to be display.
    %       Symmetric             - Logical flag to scale color symmetrically about zero.
    %       ImputeFun             - Function name or handle for imputing missing values.
    %       ColumnLabels          - Cell array of strings to label columns (x-axis).
    %       RowLabels             - Cell array of strings to label rows (y-axis).
    %       ColumnLabelsRotate    - Column labels orientation.
    %       RowLabelsRotate       - Row labels orientation.
    %       ColumnLabelsLocation  - Location of column labels.
    %       RowLabelsLocation     - Location of row labels.
    %       Annotate              - Logical flag to display value text in heat map.
    %       AnnotPrecision        - Data value display precision.
    %       AnnotColor            - Annotation text color
    %       ColumnLabelsColor     - Structure of information for coloring column labels.
    %       RowLabelsColor        - Structure of information for coloring row labels.
    %       LabelsWithMarkers     - Logical flag to show colored markers for row/column labels.
    %
    %   HeatMap methods:
    %       HeatMap     - Create a HeatMap object.
    %       view        - Views a HeapMap object in a MATLAB figure.
    %       plot        - Render a heat map
    %       addXLabel   - Add HeatMap x-axis (column) label.
    %       addYLabel   - Add HeatMap y-axis (row) label.
    %       addTitle    - Add HeatMap Graph title.
    %   
    %   Example:
    %       data = gallery('invhess',20);
    %       hm = HeatMap(data);
    %
    %   See also CLUSTERGRAM.
    
    %   Copyright 2009-2012 The MathWorks, Inc.


    properties(SetObservable=true, AbortSet=true)
        %STANDARDIZE Dimension in which the values are standardized.
        %   The Standardize property is a string specifies the dimension in
        %   which the values are standardized. The standardized values are
        %   transformed so that the mean is 0 and the standard deviation is
        %   1 in the specified direction. The dimension can be 'COLUMN'
        %   (1), 'ROW' (2), or the default 'NONE' (3). 
        %
        %   Note: The default dimension of standardization for CLUSTERGRAM
        %   is 'ROW'. 
        %
        %    See also HEATMAP.
        Standardize = 'NONE';

        %SYMMETRIC Logical flag to scale color symmetrically about zero. 
        %   The Symmetric property is a logical indicating if the color
        %   scale of the heat map to be symmetric about zero. Default is
        %   TRUE.
        %
        %    See also HEATMAP.
        Symmetric = true;
        
        %DISPLAYRANGE The range of the values to display.
        %   The DisplayRange property is a positive value specifying the
        %   display range of standardized values. For example, if you
        %   specify REDGREENCMAP for the 'COLORMAP' property and the range
        %   value set to P, pure red represents values equal to or greater
        %   than P, and pure green represents values equal to or less
        %   than -P. The default value is the maximum absolute value in
        %   input data.
        %
        %   Note: The default range for CLUSTERGRAM is 3.
        %
        %    See also HEATMAP.
        DisplayRange = 3;
        
        %COLORMAP Colormap used to display the heat map.  
        %    The Colormap property specifies the colormap used to display
        %    the heat map. It can be the name of a colormap, the function
        %    handle of a function that returns a colormap, or an M-by-3
        %    array containing RGB values. Default is REDGREENCMAP. In this
        %    color map, black represents the mean, red color represents
        %    values above the mean, and green color represents values below
        %    the mean. 
        %
        %    See also HEATMAP, REDBLUECMAP, REDGREENCMAP.
        Colormap = redgreencmap(11);
        
        %IMPUTEFUN Function name or handle for imputing missing values.
        %   The ImputeFun property is a name of a function handle of a
        %   function that imputes missing values in the input data. It can
        %   also be a cell array with the first element being the function
        %   name or handle and other elements being the input
        %   property/value pairs for the function. The missing data points
        %   are colored gray in the heat map. 
        %
        %    See also HEATMAP.
        ImputeFun = [];
        
        %COLUMNLABELS Cell array of strings to label columns (x-axis).
        %    The ColumnLabels property can be a cell array of strings or
        %    vector of numbers specifying the labels for the columns. 
        %
        %    See also HEATMAP.
        ColumnLabels ={};
        
        %ROWLABELS Cell array of strings to label rows (y-axis).
        %    The RowLabels property can be a cell array of strings or
        %    vector of numbers specifying the labels for the rows. 
        %
        %    See also HEATMAP.
        RowLabels = {};

        %COLUMNLABELSROTATE Column labels orientation.
        %    The ColumnLabelsRotate property determines the orientation of
        %    the column (X-axis) labels. Specify the value of rotation in
        %    degrees (positive angles cause counterclockwise rotation).
        %    Default is 90 degrees.
        %
        %    See also HEATMAP.
        ColumnLabelsRotate = 90;
        
        %ROWLABELSROTATE Row labels orientation.
        %    The RowLabelsRotate property determines the orientation of the
        %    row (Y-axis) labels. Specify the value of rotation in degrees
        %    (positive angles cause counterclockwise rotation). Default is
        %    0 degree.
        %
        %    See also HEATMAP.
        RowLabelsRotate = 0;
        
        %COLUMNLABELSLOCATION Location of column (X-axis) labels.
        %    The ColumnLabelsLocation property controls the display
        %    location of the column (X-axis) labels. The location can be
        %    'top' or 'bottom' (default).
        %
        %    See also HEATMAP.
        ColumnLabelsLocation = 'bottom'; % top, bottom
        
        %ROWLABELSLOCATION Location of row (Y-axis) labels.
        %    The RowLabelsLocation property controls the display location
        %    of the row (Y-axis) labels. The location can be 'right' or
        %    'left' (default).
        %
        %    See also HEATMAP.
        RowLabelsLocation = 'left'; % left, right
        
        %ANNOTATE Logical flag to display data value text in the heat map.
        %   The Annotate property is a logical specifying whether to
        %   display the data values on the heat map. The default is FALSE.
        %
        %    See also HEATMAP.
        Annotate = 'off';
        
        %ANNOTPRECISION Data value display precision.
        %   The AnnotPrecision property is a number specifying number
        %   digits of precision in the display of data values.  The default
        %   number of digits of precision is 2.
        %
        %    See also HEATMAP.
        AnnotPrecision = 2;
        
        %ANNOTCOLOR Annotation text color.
        %   The AnnotColor property is a valid Handle Graphics color
        %   specifying the color of the annotation text. It can be a
        %   three-element RGB vector or one of the predefined names. The
        %   default is white.
        %
        %    See also HEATMAP.
        AnnotColor = 'w';
        
        %COLUMNLABELSCOLOR Structure of information for coloring column labels.
        %   The ColumnLabelsColor property is a structure of information
        %   for coloring the column (X-axis) labels. The structure should
        %   contain these fields:
        %       Labels
        %       Colors
        %   Labels field contains the column labels to be colored. The
        %   elements in Colors field can be strings or three-element
        %   vectors of RGB values specifying colors, which are used to
        %   color the labels. If this field is empty, default colors are
        %   selected.
        %
        %    See also HEATMAP.
        ColumnLabelsColor = [];
        
        %ROWLABELSCOLOR Structure array of information for coloring row labels.
        %   The RowLabelsColor property is a structure array of information
        %   for coloring the row (Y-axis) labels. The structure should
        %   contain these fields:
        %       Labels
        %       Colors
        %   Labels field contains the row labels to be colored. The
        %   elements in Colors field can be strings or three-element
        %   vectors of RGB values specifying colors, which are used to
        %   color the labels. If this field is empty, default colors are
        %   selected.
        %
        %    See also HEATMAP.
        RowLabelsColor = [];
        
        %LABELSWITHMARKERS Logical flag to show colored markers for row/column labels.
        %   The LabelsWithMarkers property is a logical specifying whether
        %   to display colored marker for row/column labels. Default is
        %   FALSE, that if the annotation information are specified only
        %   the row/column labels are colored.
        %
        %    See also HEATMAP.
        LabelsWithMarkers = false;
    end
    
    properties(GetAccess='protected', SetAccess='protected',...
                                      SetObservable=true,...
                                      AbortSet=true,...
                                      Dependent=true,...
                                      Hidden = true)
                
        %XABEL X-axis label.
        %    The XLabel property is a string for label the x-axis. By
        %    default there is none. 
        %
        %    See also HeatMap.
        XLabel = '';
        
        %YABEL Y-axis label.
        %    The YLabel property is a string for label the y-axis. By
        %    default there is none. 
        %
        %    See also HEATMAP.
        YLabel = '';

        %Title HeatMap title.
        %    The Title property is a string as the title for the HeatMap.
        %    By default there is none. 
        %
        %    See also HEATMAP.
        Title = '';
    end
    
    properties(GetAccess='protected', SetAccess='protected', Hidden =true)
        Data = []; %After standardized  
        OriginalData = [];
        
        MissingDataFlag = false;
        Scales = [];
        
        CopyOnly = false;
        
        FigureHandle = [];
        HMAxesHandle = [];

        %== The ratio between the width/height of HMaxes and the annotation
        %marker axes.
        AnnotMarkerRatio = 20;
        %== Limit to show tick labels.
        TickLimit = 200;
        
        PreOrderRowLabels = [];
        PreOrderColumnLabels = [];
        PreOrderData = [];
        PreOrderOriginalData = [];
        
        TitlePVPairs = {'String', ''};
        XLabelPVPairs = {'String', ''};
        YLabelPVPairs = {'String', ''};
    end    
    
 properties(GetAccess='protected', SetAccess='protected',...
                                      SetObservable=true,...
                                      AbortSet=true )             
        %COLORBAR Insert colorbar.
        %   The Colorbar property determines whether to show the colorbar
        %   of the heatmap. The property value can be 'ON' or 'OFF'
        %   (default).
        %
        %    See also HEATMAP.
        Colorbar = 'off';
 end
    
 
    %Constructor block
    methods
        function obj = HeatMap(data, varargin)
            %HEATMAP Create a HeatMap object.
            %
            %   HEATMAP(DATA) displays a heat map representing the values
            %   in the matrix DATA.
            %
            %   HEATMAP(..., 'COLUMNLABELS',LABELS) uses the contents of a
            %   cell array of strings or a numeric array, LABELS, as labels
            %   for the column (X-axis) labels in the heat map.
            %
            %   HEATMAP(...,'ROWLABELS',LABELS) uses the contents of cell
            %   array of strings, or a numeric array, LABELS as labels for
            %   the row (Y-axis) labels in heat map.
            %
            %   HEATMAP(...,'STANDARDIZE',DIM) specifies the direction in
            %   which the values are standardized. The standardized values
            %   are transformed so that the mean is 0 and the standard
            %   deviation is 1 in the specified direction. The direction
            %   can be 'column' (1), 'row' (2), or the default 'none' (3).
            %
            %   HEATMAP(...,'COLORMAP',CMAP) allows you to specify the
            %   colormap used to display the heat map. It can be the name
            %   of a colormap, the function handle of a function that
            %   returns a colormap, or an M-by-3 array containing RGB
            %   values. Default is REDGREENCMAP.
            %
            %   HEATMAP(...,'DISPLAYRANGE', P) sets the display range of
            %   standardized values. P must be a positive scalar.  For
            %   example, if you specify REDGREENCMAP for the 'COLORMAP'
            %   property, pure red represents values equal to or greater
            %   than P, and pure green represents values equal to or less
            %   than -P. The default value is the maximum absolute value of
            %   the input data array.
            %
            %   HEATMAP(...,'SYMMETRIC', FALSE) the color scale of the heat
            %   map is not symmetric about zero. Default is TRUE.
            %
            %   HEATMAP(...,'IMPUTEFUN',FUN) allows you to specify the
            %   name or function handle of a function that imputes missing
            %   data. FUN can also be a cell array with the first element
            %   being the function name or handle and other elements being
            %   the input property/value pairs for the function. The
            %   missing data points are colored gray in the heat map.
            %
            %   HEATMAP(...,'COLUMNLABELSCOLOR', S) is an optional
            %   structure for coloring the column (X-axis) labels. The
            %   structure should contain these fields:
            %      Labels
            %      Colors
            %   The Labels field contains the column labels to be colored.
            %   The Colors field specifies the colors. The elements in the
            %   Colors field can be strings or three-element vectors of RGB
            %   values specifying colors for the column labels. 
            %
            %   HEATMAP(...,'ROWLABELSCOLOR', S) is an optional structure
            %   for coloring the row (Y-axis) labels. The structure should
            %   contain these fields:
            %       Labels
            %       Colors
            %   The Labels field contains the row labels to be colored. The
            %   Colors field specifies the colors. The elements in the
            %   Colors field can be strings or three-element vectors of RGB
            %   values specifying colors for row labels. 
            %
            %   HEATMAP(...,'LABELSWITHMARKERS', TRUE) displays colored
            %   markers instead of colored text. Default is FALSE.
            %
            %   Examples:
            %       data = bioma.data.DataMatrix(gallery('invhess',20),...
            %                     'Rownames', 'feature',... 
            %                     'Colnames', 'sample');
            %       hmap = HeatMap(data);
            %
            %   See also HEATMAP.
            
            if nargin == 0 || isempty(data)
                return;
            end 
            %== Set display data
            obj.Data = data;

            %== Handle inputs, the second argument can be a switch for if
            %   to view the heat map. (undocumented)
            viewFlag = true;
            argCnt = 1;
            if nargin > 1
                arg = varargin{argCnt};
                if islogical(arg)
                    viewFlag = arg;
                    argCnt = argCnt + 1;
                end
            end
            
            %== All pvpairs shall have an object property counterpart.
            %   Parameter names may be case insensitive but NO partial
            %   match is accepted. Once parameter names are standardized,
            %   pvpairs are just passed to the class set method.
            
            pvPairs = varargin(argCnt:end);
            if ~isempty(pvPairs)
                if rem(numel(pvPairs),2)== 1
                    error(message('bioinfo:HeatMap:HeatMap:IncorrectNumberOfArguments', mfilename))
                end
                propertyNames = fieldnames(obj);
                for j=1:2:numel(pvPairs)
                    k = bioinfoprivate.pvpair(pvPairs{j},[],propertyNames,'HeatMap:HeatMap',false);
                    pvPairs{j} = propertyNames{k};
                end
                set(obj, pvPairs{:});
            end
            
           
            %== Check for missing data
            checkMissingData(obj);
            
            %== Standardize data 
            standardizeData(obj);
            
            %== Add listeners
            addHMPropertyListeners(obj);

            %==View Heatmap
            if viewFlag
                obj.view;
            end
        end    
    end
    
    %== class method block
    methods
        function ht = addXLabel(obj, label, varargin) 
            %ADDXLABEL  Add HeatMap X-axis (column) label.
            %
            %   ADDXLABEL(HM, LABELSTR) adds the string, LABELSTR, beside
            %   the X-axis (column) of the HeatMap object HM graph display.
            %
            %   ADDXLABEL(HM, LABELSTR, 'Property1', PropertyValue1,
            %   'Property2', PropertyValue2,...) sets the values of the
            %   specified properties of the X-axis label. 
            %
            %   H = ADDXLABEL(...) returns the handle to the text object
            %   used as the X-axis label. 
            %   
            %   See also HEATMAP, ADDTITLE, ADDYLABEL.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'HeatMap:addXLabel')
            ht = [];
            if nargin > 1
               set(obj, 'XLabel', {label, varargin{:}}); %#ok
            end
            
            if ~isempty(obj.HMAxesHandle) && ishandle(obj.HMAxesHandle)
                hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');
                ht = hmdata.HMAxisTitles(1);
            end
        end
        
        function ht = addYLabel(obj, label, varargin) 
            %ADDYLABEL  Add HeatMap Y-axis (row) label.
            %
            %   ADDYLABEL(HM, LABELSTR) adds the string, LABELSTR, beside
            %   the Y-axis (row) of the HeatMap object HM graph display.
            %
            %   ADDYLABEL(HM, LABELSTR, 'Property1', PropertyValue1,
            %   'Property2', PropertyValue2,...) sets the values of the
            %   specified properties of the label. 
            %
            %   H = ADDYLABEL(...) returns the handle to the text object
            %   used as the Y-axis label. 
            %   
            %   See also HEATMAP, ADDTITLE, ADDXLABEL.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'HeatMap:addYLabel')
            ht = [];
            if nargin > 1
               set(obj, 'YLabel', {label, varargin{:}}); %#ok
            end
            
            if ~isempty(obj.HMAxesHandle) && ishandle(obj.HMAxesHandle)
                hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');
                ht = hmdata.HMAxisTitles(2);
            end
        end
        
        function ht = addTitle(obj, label, varargin) 
            %ADDTITLE  Add HeatMap Graph title.
            %
            %   ADDTITLE(HM, TITLESTR) adds the string, TITLESTR, at the
            %   top of the HeatMap object HM graph display.
            %
            %   ADDTITLE(HM, TITLESTR, 'Property1', PropertyValue1,
            %   'Property2', PropertyValue2,...) sets the values of the
            %   specified properties of the title. 
            %
            %   H = ADDTITLE(...) returns the handle to the text object
            %   used as the title. 
            %   
            %   See also HEATMAP, ADDXLABEL, ADDYLABEL.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'HeatMap:addTitle')
            ht = [];
            if nargin > 1
               set(obj, 'Title', {label, varargin{:}}); %#ok
            end
            
            if ~isempty(obj.HMAxesHandle) && ishandle(obj.HMAxesHandle)
                hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');
                ht = hmdata.HMTitleText;
            end
        end
    end
    
    %== set and get method block
    methods
        %== Get
        function data = get.Data(obj)
           data = getDisplayData(obj); 
        end
        
        function labels = get.RowLabels(obj)
            labels = getDimensionLabels(obj, 1);
        end
        
        function labels = get.ColumnLabels(obj)
            labels = getDimensionLabels(obj, 2);
            labels = labels(:)';
        end
        
        function data = get.OriginalData(obj)
            data = getOriginalData(obj);
        end
        %== Set
        function set.Data(obj, data)
            setDisplayData(obj, data)
        end
        
        function set.OriginalData(obj, data)
            setOriginalData(obj, data)
        end
        
        function set.RowLabels(obj, labels)
            if isempty(labels)
                labels = {};
            else
                labels = validateRowLabels(obj, labels);
            end
            setDimensionLabels(obj, labels, 1)
        end
        
        function set.ColumnLabels(obj, labels)
            if isempty(labels)
                labels = {};
            else
                labels = validateColumnLabels(obj, labels);
            end
            setDimensionLabels(obj, labels, 2)
        end 
        
        function set.ColumnLabelsRotate(obj, ang)
            obj.ColumnLabelsRotate = validateRotation(ang, 'ColumnLabelsRotate');
        end

        function set.RowLabelsRotate(obj, ang)
            obj.RowLabelsRotate = validateRotation(ang, 'RowLabelsRotate');
        end
        
        function set.ColumnLabelsLocation(obj, loc)
            okloc = {'top', 'bottom'};
            if ischar(loc)
                [~,obj.ColumnLabelsLocation] = bioinfoprivate.optPartialMatch(loc, okloc,...
                    'ColumnLabelsLocation','HeatMap:set');
            else
                error(message('bioinfo:HeatMap:set:ColumnLabelsLocationFormatNotValid'));
            end
        end
        
        function set.RowLabelsLocation(obj, loc)
            okloc = {'left', 'right'};
            if ischar(loc)
                [~,obj.RowLabelsLocation] = bioinfoprivate.optPartialMatch(loc, okloc,...
                    'RowLabelsLocation','HeatMap:set');
            else
                error(message('bioinfo:HeatMap:set:RowLabelsLocationFormatNotValid'));
            end            
        end
        
        function  set.XLabel(obj, label)
            obj.XLabelPVPairs = setAxesTitles(label, 'XLabel');
        end
        
        function set.YLabel(obj, label)
            obj.YLabelPVPairs = setAxesTitles(label, 'YLabel');
        end
        
        function title = get.Title(obj)          
            if ~isempty(obj.HMAxesHandle) && ishandle(obj.HMAxesHandle)
                hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');             
                title = get(hmdata.HMTitleText, 'String');
                if isempty(title)
                    title = '';
                end
            else
                title = obj.TitlePVPairs{2};
            end
        end
        
        function label = get.XLabel(obj)            
            if ~isempty(obj.HMAxesHandle) && ishandle(obj.HMAxesHandle)
                hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');
                label = get(hmdata.HMAxisTitles(1), 'String');
                if isempty(label)
                    label = '';
                end
            else
                label = obj.XLabelPVPairs{2};
            end
        end
        
        function label = get.YLabel(obj)           
            if ~isempty(obj.HMAxesHandle) && ishandle(obj.HMAxesHandle)
                hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');
                label = get(hmdata.HMAxisTitles(2), 'String');
                if isempty(label)
                    label = '';
                end
            else
                label = obj.YLabelPVPairs{2};
            end
        end
        
        function set.Title(obj, label)
            obj.TitlePVPairs = setAxesTitles(label, 'Title');
        end

        function set.Standardize(obj, dim)
            okdim = {'COLUMN', 'ROW', 'NONE'};
            if isnumeric(dim) && isscalar(dim)
                if dim == 1 || dim == 2 || dim == 3
                    dim = okdim{dim};
                else
                    dim = okdim{3};
                end
            end
            if ischar(dim)
                [~,obj.Standardize] = bioinfoprivate.optPartialMatch(dim, okdim,...
                    'Standardize','HeatMap:set');
            else
                error(message('bioinfo:HeatMap:set:StandardizeFormatNotValid'));
            end
        end
    
        function set.Symmetric(obj, t)
            obj.Symmetric = bioinfoprivate.opttf(t,'Symmetric','HeatMap:set');
        end
        
        function set.DisplayRange(obj, d)
            if isnumeric(d) && isscalar(d) && isreal(d)
                obj.DisplayRange = abs(d);
            else
                error(message('bioinfo:HeatMap:set:InvalidDisplayRangeInput'))
            end
        end
        
        function set.Colormap(obj, cm)
            validcmap = true;
            colormap = [];
            if ischar(cm) || isa(cm,'function_handle')
                try
                    colormap = eval(cm);
                catch ME %#ok
                    validcmap = false;
                end
            elseif iscell(cm) && isa(cm{1}, 'function_handle')
                try
                    colormap = feval(cm{1}, cm{2:end});
                catch ME %#ok
                    validcmap = false;
                end
            elseif isnumeric(cm)
                colormap = cm;
            else
                validcmap = false;
            end
            
            if ~validcmap || ~ismatrix(colormap) ||size(colormap,2) ~= 3
                error(message('bioinfo:HeatMap:set:InvalidColormap'))
            end
            
            obj.Colormap = colormap;
        end
        
        function set.ImputeFun(obj, x)
            if ischar(x) || isa(x,'function_handle')
                obj.ImputeFun = {x};
            elseif iscell(x)
                obj.ImputeFun = x;
            else
                error(message('bioinfo:HeatMap:set:InvalidImputeFun'));
            end
        end
        
        function set.Annotate(obj, t)
            tf = bioinfoprivate.opttf(t,'Annotate','HeatMap:set');
            if tf
                obj.Annotate = 'on';
            else
                obj.Annotate = 'off';
            end
        end
        
        function set.AnnotPrecision(obj, t)
            if isscalar(t) && isnumeric(t)
                obj.AnnotPrecision = t;
            else
                error(message('bioinfo:HeatMap:set:InvalidAnnotPrecision'));
            end
        end
        
        function set.AnnotColor(obj, t)
            if ischar(t)
                obj.AnnotColor = t;
            elseif isnumeric(t) && isvector(t) && max(size(t))== 3 &&... 
                   max(t) <= 1 && min(t)>=0
                obj.AnnotColor = t(:)';
            else
                error(message('bioinfo:HeatMap:set:InvalidAnnotColor'));
            end 
        end
        
        function set.LabelsWithMarkers(obj, t)
            obj.LabelsWithMarkers = bioinfoprivate.opttf(t,'LabelsWithMarkers','HeatMap:set');
        end
        
        function set.ColumnLabelsColor(obj, x)
            obj.ColumnLabelsColor = validateAnnotationStruct(x, obj, 'ColumnLabelsColor');
        end
        
        function set.RowLabelsColor(obj, x)
            obj.RowLabelsColor = validateAnnotationStruct(x, obj, 'RowLabelsColor');
        end
        
        function set.Colorbar(obj, t)
            tf = bioinfoprivate.opttf(t,'Colorbar','HeatMap:set');
            if tf
                obj.Colorbar = 'on';
            else
                obj.Colorbar = 'off';
            end
        end
        
        function setdisp(obj)
            %SETDISP displays all properties and property values.
            %
            %   SETDISP(HM) Special display format of the property names
            %   and their possible values.
            %
            %   See also HGSETGET.SETDISP.
            
            propertyNames = fieldnames(obj);
            propDescrs = cell2struct(cell(size(propertyNames)), propertyNames, 1);
            propDescrs.Standardize = '[column | row | {none}]';
            propDescrs.Symmetric = '[true | false].';
            propDescrs.DisplayRange = 'Scalar.';
            propDescrs.Colormap = [];
            propDescrs.ImputeFun = 'string -or- function handle -or- cell array';
            propDescrs.ColumnLabels = 'Cell array of strings, or an empty cell array';
            propDescrs.RowLabels = 'Cell array of strings, or an empty cell array';
            propDescrs.ColumnLabelsRotate = [];
            propDescrs.RowLabelsRotate = [];
            propDescrs.ColumnLabelsLocation = '[ top | {bottom} ]';
            propDescrs.RowLabelsLocation = '[ {left} | right ]';
            propDescrs.Annotate = '[on | {off}]';
            propDescrs.AnnotPrecision = [];
            propDescrs.AnnotColor = [];
            propDescrs.ColumnLabelsColor = 'A structure array.';
            propDescrs.RowLabelsColor = 'A structure array.';
            propDescrs.LabelsWithMarkers = '[true | false].';
            
            disp(propDescrs);
        end
    end %method block
    
    %== Hidden method block
    methods(Hidden = true)
        function data = getDisplayData(obj)     
            data = obj.PreOrderData;
        end
        
        function setDisplayData(obj, data) 
            obj = setDisplayDataOnly(obj, data);
            obj.DisplayRange = max(abs(obj.PreOrderData(:)));
        end
        
        function obj = setDisplayDataOnly(obj, data)
            bioma.util.validateMatrix(data, 'DATA', 'HeatMap')
            if isa(data, 'bioma.data.DataMatrix')
                obj.PreOrderData = data.(':')(':');
                obj.PreOrderColumnLabels = colnames(data);
                obj.PreOrderRowLabels = rownames(data);
            else
                obj.PreOrderData = data;
            end
        end
        
        function data = getOriginalData(obj)
            % To be overwritten by subclass
            data = obj.PreOrderOriginalData;
        end
        
        function setOriginalData(obj, data)
            % Can be overwritten by subclass
             obj.PreOrderOriginalData = data;
        end
        
        function setDimensionLabels(obj, labels, dir)
            switch dir
                case 1
                    obj.PreOrderRowLabels = labels;
                case 2
                    obj.PreOrderColumnLabels = labels;
            end
        end
        
        function labels = getDimensionLabels(obj, dir)
            switch dir
                case 1 % Row
                    labels = obj.PreOrderRowLabels;
                case 2 % Column
                    labels = obj.PreOrderColumnLabels;
            end
        end
        
        function delete(obj)
            if ishandle(obj.HMAxesHandle)
                delete(obj.HMAxesHandle)
            end
            
            if ishandle(obj.FigureHandle)
                delete(obj.FigureHandle); 
            end
        end
        
        function updateHMAxesProp(obj, src, evt) %#ok
            if ishandle(obj.FigureHandle)
                plot(obj, obj.FigureHandle, src.Name);
            end
        end
        
        function updateGraph(obj, src, evt) %#ok
            switch src.Name
                case 'Standardize'
                    standardizeData(obj);
                    if ishandle(obj.FigureHandle)
                        view(obj);
                    end
                case 'ImputeFun'
                    imputeData(obj, obj.PreOrderData);
                    if ishandle(obj.FigureHandle)
                        view(obj);
                    end
            end
        end
        
        function standardizeData(obj)
            %== Standardize data
            k = find(strncmpi(obj.Standardize, {'COLUMN', 'ROW', 'NONE'}, numel(obj.Standardize)));
            if k ~= 3
                if isempty(obj.PreOrderOriginalData)
                    obj.PreOrderOriginalData = obj.PreOrderData;
                end
                data = obj.PreOrderData;
                %== check data to be 2 dimensional matrix
                bioma.util.validateMatrix(data, 'DATA', 'standardizeData')
                isaDM = isa(data, 'bioma.data.DataMatrix');
                %== DataMatrix
                if isaDM
                    sdata = data;
                    data = data.(':')(':');
                end
                center = mean(data,k);
                scale = std(data, 0,k);
                %== Standardized data
                tscale = scale;
                %=Check for zeros and set them to 1 so not to scale them.
                scale(tscale == 0) = 1;
                %== Center and scale the data
                data = bsxfun(@minus, data, center);
                if isaDM
                    sdata.(':')(':') = bsxfun(@rdivide, data, scale);
                else
                    sdata = bsxfun(@rdivide, data, scale);
                end
                obj.PreOrderData = sdata;
                obj.Scales = tscale;
            elseif k == 3 && ~isempty(obj.PreOrderOriginalData)
                obj.PreOrderData = obj.PreOrderOriginalData;
            end
        end
        
        function checkMissingData(obj)
            % Check if the data contains missing values. If imputation
            % function is provided, impute the data.
            if isempty(obj.PreOrderOriginalData)
                obj.PreOrderOriginalData = obj.PreOrderData;
            end
            
            nandata = isnan(obj.PreOrderData);
            if ~isempty(find(nandata, 1))
                obj.MissingDataFlag = true;
                
                if ~isempty(obj.ImputeFun)
                    imputeData(obj, data)
                end
            end
        end
    
        %------------------
        function imputeData(obj, data)
            % Impute the missing data and update Data property
            try
                if numel(obj.ImputeFun) > 1
                    obj.PreOrderData = feval(obj.ImputeFun{1}, data, obj.ImputeFun{2:end});
                else
                    obj.PreOrderData = feval(obj.ImputeFun{1}, data);
                end
            catch ME 
                error(message('bioinfo:HeatMap:set:FailImputeData'))
            end
            
            if size(obj.PreOrderData) ~= size(data)
                error(message('bioinfo:HeatMap:set:MisMatchImputeDataSize'))
            end
        end
        
       
        function display(obj)
            [nrows, ncols] = size(obj.Data);
            msg = sprintf('HeatMap object with %i rows and %i columns.\n',nrows, ncols);
            disp(msg);
        end
        
        function updateFigureModes(obj, dataCursorCB)
            %== Data cursor mode
            dcmObj = datacursormode(obj.FigureHandle);
            set(dcmObj, 'UpdateFcn', {dataCursorCB, dcmObj, obj},...
                        'Enable',    'on')
            % Remove unwanted contextmenu
            hdcmMenu = get(dcmObj,'UIContextMenu'); % Get its context menu.
            delete(findobj(hdcmMenu,'Tag','DataCursorEditText'))
            delete(findobj(hdcmMenu,'Tag','DataCursorSelectText'));
            set(dcmObj,'Enable','off') %starts by default off
            
            %== Turn off roate3D
            h = rotate3d(obj.FigureHandle);
            setAllowAxesRotate(h, obj.HMAxesHandle, false)
        end

         function positionAxes(obj, hHMAxes) %#ok
            %==Do nothing in HeatMap  
         end
         
        function [varargout] = subsref(obj,s)
            %SUBSREF Subscripted reference for a HeatMap object.
            %
            %   Indexing A(I,J) and cell indexing A{I,J} are not supported.
            %
            %   P = A.PROPERTYNAME returns a HeatMap property.     
            switch s(1).type
                case '()'
                    error(message('bioinfo:HeatMap:subsref:SubscriptNotSupport'));
                case '{}'
                    error(message('bioinfo:HeatMap:subsref:CellSubscript'));
                case '.'
                    try
                        bioma.util.validateMethodProps(obj, s(1).subs);
                        [varargout{1:nargout}] = builtin('subsref', obj, s);
                    catch ME
                        bioinfoprivate.bioclsrethrow(mfilename, 'subsref', ME)
                    end
            end
        end
        
        function [varargout] = subsasgn(obj,s,b)
            %SUBSASGN Subscripted assignment to a HeatMap object.
            %
            %   Subscript indexing A(I,J) = B and cell indexing A{I,J} = B
            %   is not supported.
            %
            %   A.PROPERTYNAME = P assigns to a HeatMap property.
            
            switch s(1).type
                case '()'
                    error(message('bioinfo:HeatMap:subsasgn:SubscriptNotSupport'));

                case '{}'
                    error(message('bioinfo:HeatMap:subsasgn:CellSubscript'));

                case '.'
                    % Return default subsasgn to this object
                    try
                        bioma.util.validateMethodProps(obj, s(1).subs);
                        [varargout{1:nargout}] = builtin('subsasgn',obj,s,b);
                    catch ME
                        bioinfoprivate.bioclsrethrow(mfilename, 'subsasgn', ME)
                    end
            end %switch
        end %subsasgn
        
        function a = eq(varargin)  %#ok
            error(message('bioinfo:HeatMap:eq:UndefinedMethod'));

        end
        function a = gt(varargin)  %#ok
            error(message('bioinfo:HeatMap:gt:UndefinedMethod'));

        end
        function a = ge(varargin)  %#ok
            error(message('bioinfo:HeatMap:ge:UndefinedMethod'));

        end
        function a = lt(varargin)  %#ok
            error(message('bioinfo:HeatMap:lt:UndefinedMethod'));

        end
        function a = le(varargin)  %#ok
            error(message('bioinfo:HeatMap:le:UndefinedMethod'));

        end
        function a = ne(varargin)  %#ok
            error(message('bioinfo:HeatMap:ne:UndefinedMethod'));

        end
        function a = notify(varargin)  %#ok
            error(message('bioinfo:HeatMap:notify:UndefinedMethod'));

        end
    end %hidden methods
end % HeatMap class

%--------Helper functions
function labels = validateRowLabels(obj, labels)
% Validate RowLabels
%== Convert numeric vector to strings
if isnumeric(labels) && isvector(labels)
    labels = cellstr(num2str(labels(:)));
end
if ~iscellstr(labels)
    error(message('bioinfo:HeatMap:set:InvalidRowLabelsFormat'))
end
if size(labels, 2) ~= 1
    labels = labels';
end
if  numel(labels) ~= size(obj.Data,1)
    error(message('bioinfo:HeatMap:set:RowLabelsSizeNotMatch'));
end
end
%----------
function labels = validateColumnLabels(obj, labels)
% Validate ColumnLabels
%== Convert numeric vector to strings
if isnumeric(labels) && isvector(labels)
    labels = cellstr(num2str(labels(:)));
end

if ~iscellstr(labels)
    error(message('bioinfo:HeatMap:set:InvalidColumnLabelsFormat'))
end

if size(labels, 2) ~= 1
    labels = labels';
end

if  numel(labels) ~= size(obj.Data,2)
    error(message('bioinfo:HeatMap:set:ColumnLabelsSizeNotMatch'));
end
end
%----------
function ang = validateRotation(ang, name)
% Validate rotation angle
try
    validateattributes(ang, {'numeric'}, {'scalar'}, 'set', name)
    return;
catch ME
    bioinfoprivate.bioclsrethrow(mfilename, 'set', ME);
end
end

%--------------
function annotStruct = validateAnnotationStruct(annotStruct, obj, varName)
% Checks for color marker structure field names
% varName is 'ColumnLabelsColor' or 'RowLabelsColor'
if isempty(annotStruct)
    return;
end

if ~isstruct(annotStruct)
    switch varName
        case 'ColumnLabelsColor'
            error(message('bioinfo:HeatMap:set:InvalidColumnLabelsColorInput'));
        case 'RowLabelsColor'
            error(message('bioinfo:HeatMap:set:InvalidRowLabelsColorInput'));
    end
end

if ~(isfield(annotStruct,'Labels') && isfield(annotStruct,'Colors'))
    switch varName
        case 'ColumnLabelsColor'
            error(message('bioinfo:HeatMap:set:InvalidColumnLabelsColorField'));
        case 'RowLabelsColor'
            error(message('bioinfo:HeatMap:set:InvalidRowLabelsColorField'));
    end    
end

switch varName
    case 'ColumnLabelsColor'
        tickLabels = obj.ColumnLabels;
        numTicks = size(obj.Data, 2);
    case 'RowLabelsColor'
        tickLabels = obj.RowLabels;
        numTicks = size(obj.Data, 1);
end

if numel(annotStruct) > 1
   annotStructTmp.Labels = {annotStruct.Labels};
   annotStructTmp.Colors = {annotStruct.Colors};
   annotStruct = annotStructTmp;
end

if isnumeric(annotStruct.Labels) && isvector(annotStruct.Labels)
    % Contains tick label indices
    tf = ismember(annotStruct.Labels, 1:numTicks);
    if ~all(tf)
        error(message('bioinfo:HeatMap:set:WrongAnnotLabelIndices', num2str( annotStruct.Labels( ~tf ) )));
    end
elseif iscellstr(annotStruct.Labels)
    tf = ismember(annotStruct.Labels, tickLabels);
    if ~all(tf)
        msg = '';
        notfoundLabels = annotStruct.Labels(~tf);
        for i = 1:numel(notfoundLabels)
            msg = [msg ',' notfoundLabels{i}]; %#ok
        end
        error(message('bioinfo:HeatMap:set:WrongAnnotLabels', msg));
    end
else
    error(message('bioinfo:HeatMap:set:InvalidAnnotLabel'));
end
end

function titlePVPair = setAxesTitles(label, titleName)
% Set Title, Xlabel, YLabel properties.
% 
if ischar(label)
    titlePVPair = {'String', label};
elseif iscell(label)
    [~, args] = axescheck(label{:});
    labelStr = args{1};
    if isempty(labelStr)
        labelStr = '';
    end
    titlePVPair = {'String', labelStr, args{2:end}};
elseif isempty(label)
    titlePVPair = {'String', ''};
else
    switch titleName
        case 'XLabel'
            error(message('bioinfo:HeatMap:set:InvalidXLabelInput'))
        case 'YLabel'
            error(message('bioinfo:HeatMap:set:InvalidYLabelInput'))
        case 'Title'
            error(message('bioinfo:HeatMap:set:InvalidTitleInput'))
    end
end       
end

function addHMPropertyListeners(obj)
%== Add listeners to the properties
addlistener(obj, 'ColumnLabels', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'RowLabels', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'ColumnLabelsRotate', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'RowLabelsRotate', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'Title', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'XLabel', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'YLabel', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'ColumnLabelsLocation', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'RowLabelsLocation', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'ColumnLabelsColor', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'RowLabelsColor', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'LabelsWithMarkers', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'Colormap', 'PostSet', @obj.updateHMAxesProp);
%= Annotation texts
addlistener(obj, 'Annotate', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'AnnotPrecision', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'AnnotColor', 'PostSet', @obj.updateHMAxesProp);
%= Display range change
addlistener(obj, 'DisplayRange', 'PostSet', @obj.updateHMAxesProp);
addlistener(obj, 'Symmetric', 'PostSet', @obj.updateHMAxesProp);
%= Data changes
addlistener(obj, 'Standardize', 'PostSet', @obj.updateGraph);
addlistener(obj, 'ImputeFun', 'PostSet', @obj.updateGraph);
%=Display colorbar
addlistener(obj, 'Colorbar', 'PostSet', @obj.updateHMAxesProp);
end
