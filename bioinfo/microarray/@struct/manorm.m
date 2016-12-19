function [normVal,gnorm] = manorm(maStruct,theFieldName,varargin)
%MANORM normalizes microarray structure data
%
%   For microarray data structures that are block-based, MANORM
%   normalizes each block or print-tip region.
%
%   MANORM(MASTRUCT,FIELDNAME) normalizes field FIELDNAME for each block in
%   microarray data structure MASTRUCT using normalization method METHOD.
%   METHOD can be one of 'Mean',' Median', 'STD', 'MAD', or a function
%   handle to a normalization function.  The output is a matrix with each
%   column corresponding to the normalized data for each block.
%
%   MANORM(...,'STRUCTOUTPUT',true) returns the input structure with
%   an additional data field for the normalized data.
%
%   MANORM(...,'NEWCOLUMNNAME',COLNAME), when using STRUCTOUTPUT, allows
%   you to specify the name of the column that is appended to the list of
%   ColumnNames in the structure. The default behavior is to prefix 'Block
%   Normalized' to the FIELDNAME string.
%
%   See also MALOWESS, MAMADNORM, MAMEDIANNORM, MANORM, MASTDNORM.


% Copyright 2004-2006 The MathWorks, Inc.


% The normalization function is simple so to keep the function clean, use
% an @struct\manorm to handle struct inputs and deal with blocks.

% check that data really is block based.
dataIsBlockBased = false;
if isstruct(maStruct)
    if isfield(maStruct,'Header')&& isfield(maStruct.Header,'Type') && ...
            (~isempty(strfind(maStruct.Header.Type,'GenePix')) ...
            || ~isempty(strfind(maStruct.Header.Type,'SPOT'))...
            || ~isempty(strfind(maStruct.Header.Type,'ImaGene')))
        if nargin < 2
            theFieldName = maStruct.ColumnNames{1};
        end
        dataIsBlockBased = true;
    else
        if any(strcmp('Blocks',fieldnames(maStruct)))
            dataIsBlockBased = true;
        end
        if ~(isfield(maStruct,'Header')&& isfield(maStruct.Header,'Type') && ...
            (~isempty(strfind(maStruct.Header.Type,'Illumina')) ...
            || ~isempty(strfind(maStruct.Header.Type,'FeatureExtractor'))))
        warning(message('bioinfo:manorm:NotSupported'));
        end
    end
else
    badInputId = find(cellfun('isclass',varargin,'struct'))+2;
    if isempty(badInputId)
        badInputId = 2;
    end
    error(message('bioinfo:manorm:BadInput', badInputId));
end

if nargin < 2
    error(message('bioinfo:manorm:NoMAField'));
end
% set some defaults
structOutput = false;
manormArgs = {};

% Handle the fieldnames
if ischar(theFieldName)
    col = find(strcmpi(maStruct.ColumnNames,theFieldName));
    if isempty(col)
        fields = char(maStruct.ColumnNames);
        fields(:,end+1) = repmat(sprintf('\n'),size(fields,1),1);
        error(message('bioinfo:manorm:BadMAField', fields'));
    end
    newColName = sprintf('Block Normalized %s', theFieldName);
    theData = maStruct.Data(:,col(1));
    if numel(col)> 1
        warning(message('bioinfo:manorm:AmbiguousFieldName', theFieldName));
    end
elseif isnumeric(theFieldName)
    % treat it as column number -- force to be scalar
    col = theFieldName(1);
    try
        theData = maStruct.Data(:,col);
    catch allExceptions 
        error(message('bioinfo:manorm:BadMAColumnNumber'));
    end
    newColName = sprintf('Block Normalized Column %d', theFieldName);
else
    error(message('bioinfo:manorm:FieldnameNotString'));
end

% handle optional inputs
if nargin > 2
    if rem(nargin,2)== 1
        error(message('bioinfo:manorm:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'structoutput','newcolumnname'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            manormArgs = {manormArgs{:} , pname,pval};
        elseif length(k)>1
            error(message('bioinfo:manorm:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 %struct output
                    structOutput = bioinfoprivate.opttf(pval);
                    if isempty(structOutput)
                        error(message('bioinfo:manorm:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2 %New column name output
                    newColName = pval;
                    if ~ischar(newColName)
                        error(message('bioinfo:manorm:NewColNameNotString'));
                    end
            end
        end
    end
end

numPoints = numel(theData);
if dataIsBlockBased
    numBlocks = maStruct.Shape.NumBlocks;
else
    numBlocks = 1;
end
regularBlocks = true;
% This assumes that all blocks have the same number of points.
if dataIsBlockBased
    % check that blocks are ordered
    if issorted(maStruct.Blocks) && (rem(numPoints,numBlocks) == 0)
        theData = reshape(theData,numPoints/numBlocks,numBlocks);
    else % if not we will loop over the blocks
        if rem(numPoints,numBlocks) ~= 0 && ~structOutput
            warning(message('bioinfo:manorm:BlockSizesNotEqual'));
        end
        regularBlocks = false;
    end
end

if regularBlocks
    % call the main manorm function to do the actual normalization.
    [normVal,gnorm] = manorm(theData,manormArgs{:});
else
    % loop over the blocks
    blockNames = unique(maStruct.Blocks);
    for blockCount = 1:numel(blockNames)
        mask = (maStruct.Blocks == blockNames(blockCount));
        [normVal(mask),gnorm] = manorm(theData(mask),manormArgs{:}); %#ok
    end
end

% reshape the output if necessary
if structOutput
    maStruct.Data(:,end+1) = normVal(:);
    maStruct.ColumnNames{end+1} = newColName;
    normVal = maStruct;
end
