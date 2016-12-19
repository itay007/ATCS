function [dimIndices,numIndices,maxIndex] = getDimensionIndices(obj,dimIndices, dimN, errOnEmpty)
%GETDIMENSIONINDICES Process string, logical, or numeric DataMatrix row/column indices.
%
%   dimN =1 for rows and 2 for columns.

%   Copyright 2008-2012 The MathWorks, Inc. 


if nargin < 4
    errOnEmpty = true;
end

if dimN == 1 % Rows of DataMatrix
    objDNames = obj.RowNames;
    objDNum = obj.NRows;
    dimStr = 'Row';
    dimStrs ='rows';
elseif dimN == 2 % Columns of DataMatrix
    objDNames = obj.ColNames;
    objDNum = obj.NCols;
    dimStr = 'Column';
    dimStrs = 'columns';
else
    error(message('bioinfo:DataMatrix:getDimensionIndices:InvalidDimensionNumber'));
end

% Translate row names into indices
if ischar(dimIndices)
    if strcmp(dimIndices, ':') % already checked ischar
        % leave them alone
% %         dimIndices = 1:objDNum;
        numIndices = objDNum;
        maxIndex = objDNum;
    elseif size(dimIndices,1) == 1 || isempty(dimIndices)
        dimName = dimIndices;
        if isempty(objDNames) && isempty(dimName)
           dimIndices = 1:objDNum;
        else
           dimIndices = find(strcmp(dimName,objDNames));
        end
        
        if isempty(dimIndices) 
            if errOnEmpty
                error(message('bioinfo:DataMatrix:getDimensionIndicesUnrecognizedName', dimStrs, dimName));
            else
                numIndices = 0;
                maxIndex = [];
            end
        else
           numIndices = 1;
           maxIndex = dimIndices;
        end
    else
        error(message('bioinfo:DataMatrix:getDimensionIndices:InvalidName'));
    end
elseif iscellstr(dimIndices)
    dimNames = dimIndices;
% %     dimIndices = zeros(1,numel(dimIndices));
    dimIndices = cell(1,numel(dimIndices));
    for i = 1:numel(dimIndices)
        dimIndex = find(strcmp(dimNames{i},objDNames));
        if isempty(dimIndex)
            if errOnEmpty
                error(message('bioinfo:DataMatrix:getDimensionIndices:UnrecognizedName', dimNames{ i }, dimStrs));
            end    
        else
            dimIndices{i} = dimIndex(:)';
        end
    end
    dimIndices = cell2mat(dimIndices);
    numIndices = numel(dimIndices);
    maxIndex = max(dimIndices);
elseif isnumeric(dimIndices) || islogical(dimIndices)
    % leave the indices themselves alone
    if isnumeric(dimIndices)
        numIndices = numel(dimIndices);
        maxIndex = max(dimIndices);
    else
        numIndices = sum(dimIndices);
        maxIndex = find(dimIndices,1,'last');
    end
    
    if ~isempty(maxIndex)&& maxIndex > objDNum && ~isempty(obj)
        error(message('bioinfo:DataMatrix:getDimensionIndices:IndexOutOfBounds', dimStr));
    end
elseif iscellmat(dimIndices)
    dimIndices = dimIndices{:};
    numIndices = numel(dimIndices);
    maxIndex = max(dimIndices);
else
    error(message('bioinfo:DataMatrix:getDimensionIndices:InvalidSubscript'));
end
dimIndices = dimIndices(:);
if dimN == 2
    dimIndices = dimIndices';
end
end % DataMatrix/getDimensionIndices
