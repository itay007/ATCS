function mastruct = affytomastruct(affystruct,fieldname)
% AFFYTOMASTRUCT turns affy probe data into the same format as other microarray data

%   Copyright 2004-2009 The MathWorks, Inc.


col = strcmpi(affystruct.ProbeColumnNames,fieldname);
Xcol = strcmpi(affystruct.ProbeColumnNames,'PosX');
Ycol = strcmpi(affystruct.ProbeColumnNames,'PosY');
X = affystruct.Probes(:,Xcol);
Y = affystruct.Probes(:,Ycol);

mastruct.Shape.NumBlocks = 1;
mastruct.Shape.BlockRange = [1 1];
mastruct.Blocks=ones(size(X));
mastruct.ColumnNames = {'X','Y',fieldname};
mastruct.Data = [X Y affystruct.Probes(:,col)];
numRows = max(Y+1);
numCols = max(X+1);

% convert file indexing into MATLAB ordering -- row major
mastruct.Indices = reshape(sub2ind([numRows, numCols],Y+1,X+1),numRows,numCols);
