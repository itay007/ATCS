function [labelIndices,numIndices,maxIndex] = findLabelIndices(labels,labelIndices)
%FINDLABELINDICES Process string, logical, or numeric indices. 
% 
%   FINDLABELINDICES processes string, logical, or numeric indices in
%   labelIndices in labels - a cell array of strings. Return labelIndices,
%   numIndices and maxIndex. 

%   Copyright 2009 The MathWorks, Inc. 


labels = labels(:);
NLabels = numel(labels);

% Translate row labels into indices
if ischar(labelIndices)
    if strcmp(labelIndices, ':') % already checked ischar
        numIndices = NLabels;
        maxIndex = NLabels;
    elseif size(labelIndices,1) == 1 || isempty(labelIndices)
        theLabel = labelIndices;
        if isempty(labels) && isempty(theLabel)
           labelIndices = 1:NLabels;
        else
           labelIndices = find(strcmp(theLabel,labels));
        end
        
        if isempty(labelIndices) 
            numIndices = 0;
            maxIndex = [];
        else
           numIndices = 1;
           maxIndex = labelIndices;
        end
    else
        labelIndices = [];
        numIndices = 0;
        maxIndex = [];
    end
elseif iscellstr(labelIndices)
    theLabels = labelIndices;
    labelIndices = cell(1,numel(labelIndices));
    for i = 1:numel(labelIndices)
        labelIndex = find(strcmp(theLabels{i},labels));
        if ~isempty(labelIndex)
            labelIndices{i} = labelIndex(:)';
        end
% %         if isempty(labelIndex)
% %             labelIndices = [];
% %         else
% %             labelIndices{i} = labelIndex(:)';
% %         end
    end
    labelIndices = cell2mat(labelIndices);
    numIndices = numel(labelIndices);
    maxIndex = max(labelIndices);
elseif isnumeric(labelIndices) || islogical(labelIndices)
    % Leave the indices themselves alone
    if isnumeric(labelIndices)
        numIndices = numel(labelIndices);
        maxIndex = max(labelIndices);
    else
        numIndices = sum(labelIndices);
        maxIndex = find(labelIndices,1,'last');
    end
    
    if ~isempty(maxIndex)&& maxIndex > NLabels 
        error(message('bioinfo:findLabelIndices:IndexOutOfBounds'));

    end
elseif iscellmat(labelIndices)
    labelIndices = labelIndices{:};
    numIndices = numel(labelIndices);
    maxIndex = max(labelIndices);
else
    error(message('bioinfo:findLabelIndices:InvalidIndexType'));
end
labelIndices = labelIndices(:);

end % findLabelIndices
