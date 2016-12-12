function multialignwrite(filename,alignment,varargin)
%MULTIALIGNWRITE writes a multiple alignment file.
%
%   MULTIALIGNWRITE(FILENAME, ALIGNMENT) writes the multiple alignment
%   ALIGNMENT to file FILENAME in either the ClustalW ALN format (default)
%   or MSF format. If the file extension is .msf the file will be written
%   in MSF format otherwise it will be written in ALN format.
%
%   MULTIALIGNWRITE(..., 'FORMAT',FORMAT) allows you to specify the
%   format for the file. Valid options are 'ALN' (default) and 'MSF'.
%
%   MULTIALIGNWRITE(..., 'HEADER',HEADERTEXT) allows you to specify the
%   first line of the file. The default is 'MATLAB multiple sequence
%   alignment'.
%
%   MULTIALIGNWRITE(..., 'WRITECOUNT',TF) allows you to specify whether to
%   add the residue counts to the end of each line. Choices are true
%   (default) or false.
%
%  Example:
%
%   % Align seven cellular tumor antigen p53 sequences
%   p53 = fastaread('p53samples.txt')
%   ma = multialign(p53,'verbose',true)
%   % Write to a file called p53.aln
%   multialignwrite('p53.aln',ma)
%
%   See also FASTAWRITE, GETHMMALIGNMENT, MULTIALIGN, MULTIALIGNREAD,
%   SEQALIGNVIEWER, SEQCONSENSUS, SEQDISP, SEQPROFILE.

%   Copyright 2008-2012 The MathWorks, Inc.


% More information on the ClustalW/ALN format can be found here:
% http://www.ebi.ac.uk/help/formats.html#aln

% Some currently unsupported options

%
%   MULTIALIGNWRITE(..., 'COLUMNS',HEADER) allows you to specify the
%   first line of the file. The default is 'MATLAB multiple sequence
%   alignment'.
%
%   MULTIALIGNWRITE(..., 'NAMELENGTH',LEN) allows you to specify the
%   length of the sequence names that are written to the file. The default
%   is 15.
%
%   MULTIALIGNWRITE(..., 'STRONG',STRONGGROUPS) allows you to pass a cell
%   array of 'strong' groups. The default the 'strong' groups are:
%       {'STA', 'NEQK', 'NHQK', 'NDEQ', 'QHRK', 'MILV', 'MILF', 'HY',
%       'FYW'}
%
%   MULTIALIGNWRITE(..., 'WEAK',WEAKGROUPS) allows you to pass a cell
%   array of 'weak' groups. The default 'weak' groups are:
%   	{'CSA', 'ATV', 'SAG', 'STNK', 'STPA', 'SGND', 'SNDEQK', 'NDEQHK',
%   	'NEQHRK', 'FVLIM', 'HFY'}

bioinfochecknargin(nargin,2,mfilename)
[headerString,columns,nameLength,writeCount,strong,weak,format] = parse_inputs(varargin{:});

if ~ischar(filename),
    error(message('bioinfo:multialignwrite:FilenameMustBeString'));
end


fid = fopen(filename,'wt');
[theDir, theFile, theExtension] = fileparts(filename);

if fid == (-1)
    if ~isempty(theDir)
        error(message('bioinfo:multialignwrite:CouldNotOpenFileinDir', [ theFile, theExtension ], theDir));
    else
        error(message('bioinfo:multialignwrite:CouldNotOpenFileinPwd', filename));
    end
end
theExtension = strtok(theExtension,'.');
if isempty(format)
    if strcmpi('msf',theExtension)
        format = 'msf';
    else
        format = 'aln';
    end
else
    if ~strcmpi(format,theExtension)
        warning(message('bioinfo:multialignwrite:ExtensionFormatMismatch', format, theExtension));
    end
end
switch format
    case {'msf'}
        % write Header line
        fprintf(fid,'%s\n\n',headerString);
        
        alignLength = numel(alignment(1).Sequence);
        % decide if we are working with nucleotides or proteins
        if bioinfoprivate.isnt(alignment(1).Sequence)
            alphabet = 'N';
        else
            alphabet = 'P';
        end
        
        % check is required by some parsers but it is not clear than anyone
        % cares about the value so just fill in some number.
        fprintf(fid,'   MSF:  %d  Type: %s    Check:  1234   ..\n\n',alignLength,alphabet);
        
        % now write the names. Again check seems arbitrary as does weight.
        numSeqs = numel(alignment);
        
        
        seqHeaders = char({alignment.Header,blanks(nameLength+1)});
        seqHeaders = seqHeaders(:,1:nameLength+1);
        seqHeaders(:,end) = repmat(' ',numSeqs+1,1);
        seqHeaders(end,:)=[];
        for count = 1:numSeqs
            fprintf(fid,'Name: %s  Len:  %d  Check:  1111  Weight:  1.0\n',seqHeaders(count,:),numel(alignment(count).Sequence));
        end
        
        fprintf(fid,'\n//\n\n');
        writeCount = false;
        
        
    case {'aln'}
        
        fprintf(fid,'%s\n\n\n',headerString);
        
        alignLength = numel(alignment(1).Sequence);
        
        alignment(end+1).Sequence = multialignconserved(alignment,'strong',strong,'weak',weak);
        alignment(end).Header = '';
        
        numSeqs = numel(alignment);
        seqHeaders = char({alignment.Header,blanks(nameLength+1)});
        seqHeaders = seqHeaders(:,1:nameLength+1);
        seqHeaders(:,end) = repmat(' ',numSeqs+1,1);
        
end

numBlocks = alignLength/columns;

% Write out full blocks in loop
numFullBlocks = ceil(numBlocks-1);

pos = 1;

if writeCount
    resCount = 0;
else
    resCount =[];
end
seqCounts = cell(numSeqs,1);

for ndx = 1:numSeqs
    seqCounts{ndx} = 0;
end

for outer = 1:numFullBlocks
    
    for inner = 1:numSeqs
        seqToWrite = alignment(inner).Sequence(pos:pos+columns-1);
        seqCounts{inner} = seqCounts{inner} + numel(strrep(seqToWrite,'-',''));
        if writeCount
            resCount = seqCounts{inner};
        end
        if inner == numSeqs
            resCount = [];
        end
        fprintf(fid,'%s%s%4d\n',seqHeaders(inner,:),seqToWrite, resCount);
    end
    pos = pos+columns;
    fprintf(fid,'\n');
end

% Treat the last block as a special case
if writeCount
    resCount = alignLength;
end

for inner = 1:numSeqs
    seqToWrite = alignment(inner).Sequence(pos:end);
    seqCounts{inner} = seqCounts{inner} + numel(strrep(seqToWrite,'-',''));
    
    if writeCount
        resCount = seqCounts{inner};
    end
    if inner == numSeqs
        resCount = [];
    end
    fprintf(fid,'%s%s%4d\n',seqHeaders(inner,:),seqToWrite, resCount);
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [headerString,columns,nameLength,writeCount,strong,weak,format] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:multialignwrite:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'header','columns','namelength','writecount','strong','weak','format'};

% Set default values
headerString = 'MATLAB multiple sequence alignment';
nameLength = 15;
columns = 60;
writeCount = true;
strong = '';
weak = '';
format = '';
% Loop over the values
for j=1:2:nargin
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % headerstring
            if ~ischar(pval)
                error(message('bioinfo:multialignwrite:HeaderNotString'));
            end
            headerString = pval;
        case 2  % columns
            if ~isnumeric(pval) || ~isscalar(pval)
                error(message('bioinfo:multialignwrite:ColumnsNotScalar'));
            elseif (pval <= 0)
                error(message('bioinfo:multialignwrite:ColumnsNotPositive'));
            else
                columns = pval;
            end
        case 3  % nameLength
            if ~isnumeric(pval) || ~isscalar(pval)
                error(message('bioinfo:multialignwrite:NameLengthNotScalar'));
            elseif (pval <= 0)
                error(message('bioinfo:multialignwrite:NameLengthNotPositive'));
            else
                nameLength = pval;
            end
        case 4  % writeCount
            
            writeCount = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 5  % strong
            strong = pval;
        case 6  % weak
            weak = pval;
        case 7  % format
            if ~ischar(pval)
                error(message('bioinfo:multialignwrite:FormatNotString'));
            end
            okformats = {'aln','msf'};
            format = lower(pval);
            if ~ismember(format,okformats)
                error(message('bioinfo:multialignwrite:FormatNotSupported', format));
            end
    end
    
end
