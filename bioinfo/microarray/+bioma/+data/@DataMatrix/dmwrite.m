function dmwrite(obj, filename, varargin)
%DMWRITE Write DataMatrix object to a text file.
%
%   DMWRITE(DM, FILENAME) writes a DataMatrix object, DM, to a text file,
%   FILENAME using the delimiter (\t) to separate DataMatrix elements. The
%   data is written starting at the first column of the first row in the
%   destination file. FILENAME must be a string.
%
%   DMWRITE(..., 'DELIMITER', DELIMITER) specifies delimiter character
%   DELIMITER. Default is '\t'.
%
%   DMWRITE(..., 'PRECISION',PS) specifies the numeric precision, PS, for
%   writing data to the file. PS can be the number of significant digits or
%   a C-style format string starting in %, such as '%6.5f'. Default is 5.
%
%   DMWRITE(..., 'HEADER', HEADERTEXT) specifies the first line of the
%   file. Default is the name in the Name property of the DataMatrix
%   object.
%
%   DMWRITE(..., 'ANNOTATED', FALSE) does not write row and column names to
%   the file. Default is TRUE.
%
%   DMWRITE(..., 'APPEND', TRUE) appends the DataMatrix object to the end
%   of the specified existing file. Default is FALSE.
%
%  Example:
%
%   % Create a DataMatrix object d 
%   d = bioma.data.DataMatrix(rand(2,3), {'Row1', 'Row2'}, {'Col1', 'Col2', 'Col3'})
%   % Write d to a text file
%   dmwrite(d,'mydm.txt')
%
%   See also DATAMATRIX.

%   Copyright 2009-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename)

%== Check input arguments
inPV = parse_inputs(varargin{:});

if ~ischar(filename),
    error(message('bioinfo:DataMatrix:dmwrite:FilenameMustBeString'));
end

if isempty(obj.RowNames) && isempty(obj.ColNames)
    inPV.Annotated = false;
end

if inPV.Append
    fid = fopen(filename, 'at');
else
    fid = fopen(filename,'wt');
end
[theDir, theFile, theExtension] = fileparts(filename);

if fid == (-1)
    if ~isempty(theDir)
        error(message('bioinfo:DataMatrix:dmwrite:CouldNotOpenFileinDir', [ theFile, theExtension ], theDir));
    else
        error(message('bioinfo:DataMatrix:dmwrite:CouldNotOpenFileinPwd', filename));
    end
end

try
    if isnumeric(inPV.Precision);
        format = sprintf('%%.%df%s', inPV.Precision, inPV.Delimiter);
    else
        format = sprintf('%s%s', inPV.Precision, inPV.Delimiter);
    end
    if inPV.Annotated
        if isempty(inPV.Header)
            inPV.Header = obj.Name;
        end
        
        if ~isempty(inPV.Header)
            fprintf(fid,'%s\n', inPV.Header);
        end
        %== Write column names
        if ~isempty(obj.ColNames)
            fwrite(fid,inPV.Delimiter, 'uchar');
            for iloop = 1:obj.NCols-1
                fprintf(fid, '%s%s', obj.ColNames{iloop}, inPV.Delimiter);
            end
            fprintf(fid, '%s%s',obj.ColNames{obj.NCols}, inPV.Newline);
        end
        %== Write one row at a time
        for iloop = 1:obj.NRows
            str1 = sprintf('%s%s', obj.RowNames{iloop}, inPV.Delimiter);
            str2 = sprintf(format, obj.Matrix(iloop, :));
            if inPV.isDlmChar
                fprintf(fid, '%s%s%s', str1, str2(1:end-1), inPV.Newline);
            else
                fprintf(fid, '%s%s%s', str1, str2, inPV.Newline);
            end
        end
    else % just the numbers
        %== Write one row at a time
        for iloop = 1:obj.NRows
            str = sprintf(format, obj.Matrix(iloop, :));
            if inPV.isDlmChar
                fprintf(fid, '%s%s', str(1:end-1), inPV.Newline);
            else
                fprintf(fid, '%s%s', str, inPV.Newline);
            end
        end
    end
catch ME
    fclose(fid);
    bioinfoprivate.bioclsrethrow('DataMatrix','dmwrite', ME)
end
%== Close the file
fclose(fid);
end

%---Helper functions--------
function inputStruct = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

%==Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:DataMatrix:dmwrite:IncorrectNumberOfArguments', mfilename))
end

% The allowed inputs
okargs = {'delimiter', 'precision', 'header', 'annotated', 'append'};

% Defaults
inputStruct.Delimiter = sprintf('\t');   % tab-delimited
inputStruct.Precision = 5;      % precision
inputStruct.Header = '';        % Header text - one line string
inputStruct.Annotated = true;   % annotated with row and column names
inputStruct.Append = false;     % append data to the end of the file
inputStruct.isDlmChar = false;  % delimiter character is not '\t' or space
inputStruct.Newline = sprintf('\n');  % newline '\n'

for j=1:2:nargin
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs,'DataMatrix:dmwrite');
    switch(k)
        case 1 % Delimiter
            tmp = sprintf(pval);
            if ischar(pval) && length(tmp) <= 1
                inputStruct.Delimiter = tmp;
                inputStruct.isDlmChar = ~strcmp(pval, '\t') && ~isempty(strtrim(pval));
            else
                error(message('bioinfo:DataMatrix:dmwrite:Invaliddelimiter', pval));
            end
        case 2 % Precision
            if (isnumeric(pval) && isscalar(pval)) || ischar(pval)
                inputStruct.Precision = pval;
            else
                error(message('bioinfo:DataMatrix:dmwrite:InvalidPrecision'));
            end
        case 3  % Header string
            if ~ischar(pval)
                error(message('bioinfo:DataMatrix:dmwrite:HeaderNotString'));
            end
            inputStruct.Header = pval;
        case 4 % annotated
            inputStruct.Annotated = bioinfoprivate.opttf(pval,okargs{k},'DataMatrix:dmwrite');
        case 5  % append
            inputStruct.Append = bioinfoprivate.opttf(pval,okargs{k},'DataMatrix:dmwrite');
    end
    
end
end
