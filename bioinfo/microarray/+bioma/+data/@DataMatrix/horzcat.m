function obj = horzcat(varargin)
%HORZCAT Horizontal concatenation of DataMatrix objects.
%
%   DM = HORZCAT(DM1, DM2,...) horizontally concatenates the DataMatrix
%   objects DM1 and DM2 into a DataMatrix object DM. DM1 and DM2 must have
%   same number of rows. The row names and order of DM are the same as DM1.
%   The row names of DM2 and any other DataMatrix input arguments are not
%   preserved. The columns names for DM are the column names of DM1 and DM2
%   and other DataMatrix object input arguments. 
%
%   DM = HORZCAT(DM1, B,...) horizontally concatenates the DataMatrix
%   object DM1 and MATLAB numeric or logical array B into a DataMatrix
%   object DM. DM1 and B must have same number of rows. The row names and
%   order of DM are the same as DM1. The columns names for DM are the
%   column names of DM1 and empty for the columns of B.
% 
%   Y = HORZCAT(X1,X2,X3,...) is called for the syntax '[X1 X2 X3 ...]'
%   when any of X1, X2, X3, etc. is a DataMatrix object.
%
%   See also DATAMATRIX/CAT, DATAMATRIX/VERTCAT.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Check input type
dmIdx = 0;
for i = 1:nargin
    b = varargin{i};
    if isa(b, 'bioma.data.DataMatrix') && dmIdx == 0
           obj = b;
           dmIdx = i;
    elseif ~isa(b, 'bioma.data.DataMatrix') && ~isnumeric(b) && ~islogical(b)
           error(message('bioinfo:DataMatrix:horzcat:InvalidInput'));      
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
        elseif ~isa(b, 'bioma.data.DataMatrix') 
            obj = bioma.data.DataMatrix(b);
        end
        continue;
    elseif isempty(b)
        % Do nothing
        continue;
    elseif isa(b, 'bioma.data.DataMatrix') && obj.NRows ~= b.NRows
        error(message('bioinfo:DataMatrix:horzcat:SizeMismatch'));
    elseif ~isa(b, 'bioma.data.DataMatrix') && obj.NRows ~= size(b,1)
        error(message('bioinfo:DataMatrix:horzcat:SizeMismatch'));
    end

    %== Concatenate
    if isa(b, 'bioma.data.DataMatrix')
        obj.Matrix = [obj.Matrix b.Matrix];
        oldNCols = obj.NCols;
        %== Update column number
        obj.NCols = size(obj.Matrix, 2);
        if isempty(obj.ColNames) && isempty(b.ColNames)
            % still empty
        elseif isempty(obj.ColNames) && ~isempty(b.ColNames)
            obj.ColNames = [repmat({''}, 1, oldNCols) b.ColNames];
        elseif ~isempty(obj.ColNames) && isempty(b.ColNames)
            obj.ColNames = [obj.ColNames repmat({''}, 1, b.NCols)];
        else
            obj.ColNames = [obj.ColNames b.ColNames];
        end
    else
        obj.Matrix = [obj.Matrix b]; 
        obj.NCols = size(obj.Matrix, 2);
        if isempty(obj.ColNames)
            % still empty
        else
            obj.ColNames = [obj.ColNames repmat({''}, 1, size(b, 2))];
        end
    end
    
end % for

end % DataMatrix/horzcat mathod
