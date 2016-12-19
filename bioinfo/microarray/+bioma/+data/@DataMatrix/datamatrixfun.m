function [varargout] = datamatrixfun(fun,a,varargin)
%DATAMATRIXFUN Apply a function to each element of a DataMatrix object.
% 
%   A = DATAMATRIXFUN(FUN, B) applies the function specified by FUN to each
%   element of DataMatrix B, and returns the results in a DataMatrix A. A
%   is the same size as B and has the same dimension properties, and the
%   (I,J,...)th element of A is equal to FUN(B(I,J,...)). FUN is a function
%   handle to a function that takes one input argument and returns a scalar
%   value. FUN must return values of the same class each time it is called.
%   See the help for ARRAYFUN for more details about FUN.
% 
%   A = DATAMATRIXFUN(FUN, B, C,  ...) evaluates FUN using elements of
%   DataMatrix object B and DataMatrix or MATLAB arrays C,  ... as input
%   arguments.  B, C, ... must all have the same size. The (I,J,...)th
%   element of A is equal to FUN(B(I,J,...), C(I,J,...), ...). If more than
%   one input argument is DataMatrix class. A will have B's dimensional
%   properties with the same row names and column names.
%
%   [A, B, ...] = DATAMATRIXFUN(FUN, C, ...), where FUN is a function
%   handle to a function that returns multiple outputs, returns arrays A,
%   B, ..., each corresponding to one of the output arguments of FUN.
%   DMARRAYFUN calls FUN each time with as many outputs as there are in the
%   call to ARRAYFUN.  FUN may return output arguments having different
%   classes, but the class of each output must be the same each time FUN is
%   called.
% 
%   [A, ...] = DATAMATRIXFUN(FUN, B,  ..., 'UNIFORMOUTPUT', FALSE) allows
%   you to specify a function FUN that returns values of different sizes or
%   types. DMARRAYFUN returns a cell array (or multiple cell arrays), where
%   the (I,J,...)th cell contains the value FUN(C(I,J,...), ...).  When 
%   'UniformOutput' is false, the outputs can be of any type. Setting
%   'UniformOutput' to TRUE is equivalent to the default behavior.
% 
%   [A,...] = DATAMATRIXFUN(FUN, B,  ..., 'DATAMATRIXOUTPUT', FALSE)
%   returns the output(s) of FUN without change them to DataMatrix class
%   type. The type of outputs is determined by 'UNIFORMOUTPUT' option.
%   Setting 'DATAMATRIXOUT' to TRUE is equivalent to the default behavior.
%
%   [B,...] = DATAMATRIXFUN(FUN, A,  ..., 'ROWS', ROWNAMES) allows you to
%   apply FUN only to the DataMatrix rows in A specified by ROWNAMES.
%   ROWNAMES is a positive integer, a vector of positive integers, a row
%   name, a cell array containing one or more row names, or a logical
%   vector.
%
%   [B,...] = DATAMATRIXFUN(FUN, A,  ..., 'COLUMNS', COLNAMES) allows you
%   to apply FUN only to the DataMatrix columns in A specified by COLNAMES.
%   COLNAMES is a positive integer, a vector of positive integers, a column
%   name, a cell array containing one or more column names, or a logical
%   vector.
%
%   [B, ...] = DATAMATRIXFUN(FUN, A,  ..., 'ERRORHANDLER', EFUN), where
%   EFUN is a function handle, specifies the function for MATLAB to call if
%   the call to FUN fails.  The error handling function will be called with
%   the following input arguments:
%     -  a structure, with the fields:  "identifier", "message", and
%        "index", respectively containing the identifier of the error that
%        occurred, the text of the error message, and the linear index into
%        the input array(s) at which the error occurred. 
%     -  the set of input arguments at which the call to the function 
%        failed.
%    See the help for ARRAYFUN for more details about EFUN.
%
%   See also ARRAYFUN DATAMATRIX/DMBSXFUN

%   Copyright 2008-2012 The MathWorks, Inc.


%== Input checking
bioinfochecknargin(nargin,2,'DataMatrix:datamatrixfun')

%== The first input to the function handle FUN must be a DataMatrix object
if ~isa(a, 'bioma.data.DataMatrix')
    error(message('bioinfo:DataMatrix:datamatrixfun:FirstInputMustBeDataMatrix'))
end

%== Find the first param name/value
argCount = 0;
while argCount < nargin-2
    arg = varargin{argCount+1};
    if bioma.util.isString(arg)
        break;
    end
    argCount = argCount + 1;
end
inStruct = parse_inputs(a, varargin{argCount+1:end});
%== If output of FUN is not uniform output, not DataMatrix output.
if ~inStruct.UNOutput
    inputStruct.DMOutput = false;
end

aMatrix = a.Matrix(inStruct.Rows, inStruct.Columns);
try
    if argCount > 0
        for i = 1:argCount
            if isa(varargin{i}, 'bioma.data.DataMatrix')
                varargin{i} = varargin{i}.Matrix;
            end

            if isequal(size(a.Matrix), size(varargin{i}))
                dmatrix = varargin{i};
                varargin{i} = dmatrix(inStruct.Rows, inStruct.Columns);
            end

        end

        if isempty(inStruct.ErrFun)
            [varargout{1:nargout}] = arrayfun(fun,aMatrix,varargin{1:argCount},...
                'UniformOutput',inStruct.UNOutput);
        else
            [varargout{1:nargout}] = arrayfun(fun,adata,varargin{1:argCount},...
                'UniformOutput',inStruct.UNOutput,'ErrorHandler',inStruct.ErrFun);
        end
    else
        if isempty(inStruct.ErrFun)
            [varargout{1:nargout}] = arrayfun(fun,aMatrix,...
                'UniformOutput',inStruct.UNOutput);
        else
            [varargout{1:nargout}] = arrayfun(fun,adata,...
                'UniformOutput',inStruct.UNOutput,'ErrorHandler',inStruct.ErrFun);
        end

    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','datamatrixfun', ME);
end 

%== DataMatrix output only if UniformOutput is true
if inStruct.DMOutput && inStruct.UNOutput
    for i = 1:nargout
        bMatrix= varargout{i};
        rowNames =a.RowNames;
        if ~isempty(a.RowNames)
            rowNames = a.RowNames(inStruct.Rows);
        end
        
        colNames =a.ColNames;
        if ~isempty(a.ColNames)
            colNames = a.ColNames(inStruct.Columns);
        end
         b = bioma.data.DataMatrix(bMatrix, rowNames,colNames);

        varargout{i} = b;
    end
end
end % DataMatrix/datamatrixfun

function inputStruct = parse_inputs(obj, varargin)
% Check for the right number of inputs
if rem(nargin-1, 2)== 1
    error(message('bioinfo:DataMatrix:datamatrixfun:IncorrectNumberOfArgs'))
end

% The allowed inputs
okargs = {'uniformoutput', 'datamatrixoutput','rows','columns','errorhandler'};

% Defaults
inputStruct.UNOutput = true;
inputStruct.DMOutput = true;   % DataMatrix output flag
inputStruct.Rows = 1:obj.NRows;          % rownames 
inputStruct.Columns = 1:obj.NCols;       % columnnames
inputStruct.ErrFun = [];       % Error handler

for i=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{i}, varargin{i+1},...
                                      okargs, 'DataMatrix:datamatrixfun');
    switch(k)
        case 1 % uniformaoutput
            inputStruct.UNOutput =  bioinfoprivate.opttf(pval,okargs{k},...
                                                'DataMatrix:datamatrixfun');
        case 2 %datamatrixoutpt
            inputStruct.DMOutput = bioinfoprivate.opttf(pval,okargs{k},...
                                                'DataMatrix:datamatrixfun');
        case 3 % rows
            inputStruct.Rows = getDimensionIndices(obj,pval, 1, false);
            if isempty(inputStruct.Rows)
               inputStruct.Rows = 1:obj.NRows; 
            end
        case 4 % columns
            inputStruct.Columns = getDimensionIndices(obj,pval, 2, false);
            if isempty(inputStruct.Columns)
               inputStruct.Columns = 1:obj.NCols; 
            end
        case 5 % error handler
            inputStruct.ErrFun = pval;
    end %switch
end 
end % parse_inputs
