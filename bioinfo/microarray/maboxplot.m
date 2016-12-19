function hOut = maboxplot(maStruct,theFieldName,varargin)
%MABOXPLOT displays a box plot of microarray data.
%
%   MABOXPLOT(DATA) displays a box plot of the values in the columns of
%   DATA. DATA can be a numeric array, a DataMatrix object or a structure
%   containing a field called Data. If DATA is a DataMatrix object, its
%   column names labels the box plot.
%
%   MABOXPLOT(DATA,COLUMN_NAMES) labels the box plot column names.
%
%   For microarray data structures that are block-based, MABOXPLOT creates
%   a box plot of a given field for each block.
%
%   MABOXPLOT(MASTRUCT,FIELDNAME) displays a box plot of field FIELDNAME
%   for each block in microarray data structure MASTRUCT.
%
%   MABOXPLOT(...,'TITLE',TITLE) allows you to specify the title of the
%   plot. The default title is FIELDNAME.
%
%   MABOXPLOT(...,'NOTCH',TF) draws notched boxes when TF is set to TRUE.
%   The default is to show square boxes and TF is FALSE.
%
%   MABOXPLOT(...,'SYMBOL',SYM) allows you to specify the symbol used for
%   outlier values. The default symbol is '+'.
%
%   MABOXPLOT(...,'ORIENTATION',ORIENT) allows you to specify the
%   orientation of the box plot. The choices are 'Vertical' and
%   'Horizontal'. The default is 'Vertical'.
%
%   MABOXPLOT(...,'WHISKERLENGTH',WHIS) allows you to specify the whisker
%   length for the box plot. WHIS defines the maximum length of the
%   whiskers as a function of the interquartile range (IQR). The default
%   length is 1.5. The whisker extends to the most extreme data value
%   within WHIS*IQR of the box. If WHIS is 0, then BOXPLOT displays all
%   data values outside the box using the plotting symbol SYM.
% 
%   MABOXPLOT(...,'BOXPLOT',BOXARGS) allows you to pass arguments to the
%   BOXPLOT, the function used to create the box plot. BOXARGS should be a
%   cell arrays of parameter/value pairs that can be passed to BOXPLOT. See
%   the help for BOXPLOT for more details of the available options.
%
%   H = MABOXPLOT(...) returns the handle of the box plot axes.
%
%   [H,HLINES] = MABOXPLOT(...) returns the handles of the lines used to
%   separate the different blocks in the image.
%
%   Examples:
%
%       load yeastdata
%       maboxplot(yeastvalues,times);
%       xlabel('Sample Times');
%
%       % Using a structure
%       geoStruct = getgeodata('GSM1768');
%       maboxplot(geoStruct,'title','GSM1768');
%
%       % For block-based data
%       madata = gprread('mouse_a1wt.gpr');
%       maboxplot(madata,'F635 Median - B635','TITLE','Cy5 Channel FG - BG');
%
%   See also AFFYPREPROCESSDEMO, BOXPLOT, MAGETFIELD, MAIMAGE, MAIRPLOT,
%   MALOGLOG, MALOWESS, MANORM, MAVOLCANOPLOT, MOUSEDEMO.

% Copyright 2003-2012 The MathWorks, Inc.


% check inputs
import bioinfoprivate.*;
bioinfochecknargin(nargin,1,mfilename);

emptyFieldName = false;
numArgs = nargin;

if numArgs < 2
    emptyFieldName = true;
end

% check for the case where the inputs are (DATA,param,value,...) with no
% fieldname.
okargs = {'title','notch','symbol','orientation','whiskerlength','boxplot'};
if numArgs > 2 && rem(numArgs,2)== 1
    if any(strncmpi(theFieldName,okargs,numel(theFieldName)))
        varargin = {theFieldName , varargin{:}};
        theFieldName = [];
        emptyFieldName = true;
        numArgs = numArgs +1;
    end
end

dataIsBlockBased = false;
if ~isstruct(maStruct)
    if isa(maStruct, 'bioma.data.DataMatrix')
        data = maStruct.(':')(':');
        
        if emptyFieldName
            ColumnNames = colnames(maStruct);
            if isempty(ColumnNames)
                ColumnNames = 1:maStruct.NCols;
            end
        else
            ColumnNames = theFieldName;
        end
    else
        data = maStruct;
        
        if emptyFieldName
            ColumnNames = 1:size(data,2);
        else
            ColumnNames = theFieldName;
        end
    end
else
    if isfield(maStruct,'FullPathName') && ...
            ~isempty(strfind(upper(maStruct.FullPathName),'CHP'))
        % Affymetrix data
        if emptyFieldName
            fieldname = 'Signal';
        else
            fieldname = theFieldName;
        end
        try
            data = log2([maStruct.ProbeSets.(fieldname)]);
        catch ME %#ok
            error(message('bioinfo:maboxplot:BadAffyField'))
        end
        ColumnNames = sprintf('Log2 of %s',fieldname);
        
    elseif checkmastruct(maStruct) % GenePix, Imagene etc.
        if emptyFieldName
            theFieldName = maStruct.ColumnNames{1};
        end
        if strcmpi(maStruct.Header.Type,'FeatureExtractor')
            maStruct = addIndicesAndBlocks(maStruct);
        end
        dataIsBlockBased = true;

    elseif isfield(maStruct,'Header')&& isfield(maStruct.Header,'Type')&& ...
            (~isempty(strfind(maStruct.Header.Type,'Gene Expression Omnibus'))...
            ||~isempty(strfind(maStruct.Header.Type,'Illumina')))
        data = maStruct.Data;
        ColumnNames = maStruct.ColumnNames;
        if iscell(data)
            numCols = size(data,2);
            numRows = size(data,1);
            isNum = true(numCols,1);
            for colCheck = 1:numCols
                if ~isnumeric(maStruct.Data{1,colCheck})
                    isNum(colCheck) = false;
                end
            end
            isNumIndices = find(isNum);
            data = reshape(cell2mat({maStruct.Data{:,isNumIndices}}),numRows,numel(isNumIndices));
            ColumnNames = ColumnNames(isNumIndices);
        end
        
    elseif isfield(maStruct,'Data')
        warning(message('bioinfo:maboxplot:UnknownFormat'));
        data =  maStruct.Data;
        ColumnNames = '';
    else
        error(message('bioinfo:maboxplot:NotSupported'));
    end
end

titleString = '';
showNotch = false;
symbol = 'r+';
vertical = true;
whiskerLength = 1.5;
boxargs = {};

if numArgs > 2
    if rem(numArgs,2)== 1
        error(message('bioinfo:maboxplot:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'title','notch','symbol','orientation','whiskerlength','boxplot'};
    
    for j=1:2:numArgs-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        
        switch(k)
            case 1 % title
                titleString = pval;
            case 2 % notch
                showNotch = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 3 % symbol
                if ~isempty(pval) && ischar(pval)
                    symbol = pval(1);
                else
                    warning(message('bioinfo:maboxplot:BadSymbol'));
                end
            case 4 %orientation
                if ~isempty(pval) && ischar(pval)
                    if(lower(pval(1)) == 'h')
                        vertical = false;
                    end
                else
                    warning(message('bioinfo:maboxplot:BadOrientation'));
                end
            case 5 %whiskerlength
                whiskerLength = pval;
            case 6 % boxplot
                if iscell(pval)
                    boxargs = pval;
                else
                    boxargs = {pval};
                end
        end
    end
end
if showNotch
    notch = 'on';
else
    notch = 'off';
end

if vertical 
    orien = 'vertical';
else
    orien = 'horizontal';
end

if dataIsBlockBased

    col = find(strcmpi(maStruct.ColumnNames,theFieldName));

    if isempty(col)
        fields = char(maStruct.ColumnNames);
        fields(:,end+1) = repmat(sprintf('\n'),size(fields,1),1);
        error(message('bioinfo:maboxplot:BadMAField', sprintf('%s',fields')));
    end

    theData = maStruct.Data(:,col);
    numPoints = numel(theData);
    numBlocks = maStruct.Shape.NumBlocks;

    % data should be sorted by blocks, but we should make sure
    if ~issorted(maStruct.Blocks,'rows')
        [sBlocks,perm] = sortrows(maStruct.Blocks);
        theData = theData(perm);
    end
    theData = reshape(theData,numPoints/numBlocks,numBlocks);

    uniqueBlocks = unique(maStruct.Blocks,'Rows');
    boxplot(theData,showNotch,symbol,vertical,whiskerLength);
    hAxis = gca;

    if size(uniqueBlocks,2) > 1
        tickLabels = cell(numBlocks,1);
        for count = 1:numBlocks
            tickLabels{count} = regexprep(['{' num2str(uniqueBlocks(count,:)) '}'],' .?',',');
        end
    else
        tickLabels = num2str((1:numBlocks)');
    end

    boxplot(theData,'Notch', notch,...
                    'Symbol',symbol,...
                    'Orientation', orien,...
                    'Whisker', whiskerLength,...
                    'Labels', tickLabels,...
                    boxargs{:});

    if vertical
        xlabel('Block');
        ylabel(theFieldName)
    else
        ylabel('Block');
        xlabel(theFieldName)
    end

    if isempty(titleString)
        titleString = theFieldName;
    end
    title(titleString);

    if nargout > 0
        hOut = hAxis;
    end

else
    if isnumeric(ColumnNames)
        ColumnNames = cellstr(num2str((ColumnNames(:))));
    end
    boxplot(data,'Notch', notch,...
            'Symbol',symbol,...
            'Orientation', orien,...
            'Whisker', whiskerLength,...
            'Labels', ColumnNames,...
            boxargs{:});

    hAxis = gca;
    if vertical
        xlabel('');
    else
        ylabel('');
    end

    if ~isempty(titleString)
        title(titleString);
    end

    if nargout > 0
        hOut = hAxis;
    end
end
end % maboxplot

