function obj = vertcat(varargin)
%VERTCAT Vertical concatenation of DataMatrix objects.
%
%   DM = VERTCAT(DM1, DM2,...) vertically concatenates the DataMatrix
%   objects DM1 and DM2 into a DataMatrix object DM. DM1 and DM2 must have
%   same number of columns. The column names and order of DM are the same
%   as DM1. The column names of DM2 and any other DataMatrix input
%   arguments are not preserved. The row names for DM are the row names of
%   DM1 and DM2, and other DataMatrix input arguments. 
%
%   DM = VERTCAT(DM1, B,...) vertically concatenates the DataMatrix object
%   DM1 and MATLAB numeric or logical array B into a DataMatrix object DM.
%   DM1 and B must have same number of columns. The column names and order
%   of DM are the same as DM1. The rows names for DM are the row names of
%   DM1 and empty for the rows of B.
% 
%   Y = VERTCAT(X1,X2,X3,...) is called for the syntax '[X1; X2; X3; ...]'
%   when any of X1, X2, X3, etc. is a DataMatrix object.
%
%   See also DATAMATRIX/CAT, DATAMATRIX/HORZCAT.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Check input type
dmIdx = 0;
for i = 1:nargin
    b = varargin{i};
    if isa(b, 'bioma.data.DataMatrix') && dmIdx == 0
           obj = b;
           dmIdx = i;
    elseif ~isa(b, 'bioma.data.DataMatrix') && ~isnumeric(b) && ~islogical(b)
           error(message('bioinfo:DataMatrix:vertcat:InvalidInput'));      
   end
end

%== Loop and cancatenate 
for i = 1:nargin
    if i == dmIdx
        continue;
    end
    
    b = varargin{i};
    
    %= Handle empty DataMatrix cases
    if isempty(obj)
        if isa(b, 'bioma.data.DataMatrix') 
            obj = b;
        else
            obj = bioma.data.DataMatrix(b);
        end
        continue;
    elseif isempty(b)
        % Do nothing
        continue;
    elseif isa(b, 'bioma.data.DataMatrix') && obj.NCols ~= b.NCols
        error(message('bioinfo:DataMatrix:vertcat:SizeMismatch'));
    elseif ~isa(b, 'bioma.data.DataMatrix') && obj.NCols ~= size(b,2)
        error(message('bioinfo:DataMatrix:vertcat:SizeMismatch'));
    end
    
    %== Concatenate
    if isa(b, 'bioma.data.DataMatrix')
        obj.Matrix = [obj.Matrix; b.Matrix];
        oldNRows = obj.NRows;
        %== Update row number
        obj.NRows = size(obj.Matrix, 1);
        if isempty(obj.RowNames) && isempty(b.RowNames)
            % still empty
        elseif isempty(obj.RowNames) && ~isempty(b.RowNames)
            obj.RowNames = [repmat({''}, oldNRows, 1); b.RowNames];
        elseif ~isempty(obj.RowNames) && isempty(b.RowNames)
            obj.RowNames = [obj.RowNames; repmat({''}, b.NRows, 1)];
        else
            obj.RowNames = [obj.RowNames; b.RowNames];
        end
    else
        obj.Matrix = [obj.Matrix; b]; 
        obj.NRows = size(obj.Matrix, 1);
        if isempty(obj.RowNames)
            % still empty
        else
            obj.RowNames = [obj.RowNames; repmat({''}, size(b, 1),1)];
        end
    end
    
end % for

end % DataMatrix/vertcat mathod
