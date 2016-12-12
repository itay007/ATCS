function data = soapread(filename, varargin)
%SOAPREAD reads SOAP format file.
%
%   DATA = SOAPREAD(FILENAME) reads a SOAP format file FILENAME, returning
%   the data in the file as a structure. FILENAME can also be a MATLAB
%   string that contains the text of a SOAP format file. The output
%   DATA will be an Nx1 structure, where N is the number of alignment
%   records stored in the SOAP file.
%
%   SOAPREAD(...,'BLOCKREAD', M) allows you to read in a single entry or
%   block of entries from a file containing multiple alignment records. If
%   M is a scalar then the M'th entry in the file is read. If M is a two
%   element vector then the block of entries starting at entry M(1) and
%   ending at entry M(2) will be read.  Use Inf for M(2) to read all
%   entries in the file starting at position M(1).
% 
%   SOAPREAD(...,'ALIGNDETAILS', false) will skip the additional field
%   called AlignDetails, which contains information on mismatches,
%   insertions and deletions in the alignment.  Default is true.
% 
%   Examples:
%
%       % Read the alignment records stored in 'sample01.soap'
%       data = soapread('sample01.soap');
%
%       % Read a block of entries from a file
%       data = soapread('sample01.soap','blockread', [5 10]);
%
%   See also BOWTIEREAD, SAMREAD, FASTQREAD

%   Copyright 2010-2012 The MathWorks, Inc.


% For a detailed description of the SOAP alignment format, see
% http://soap.genomics.org.cn/

bioinfochecknargin(nargin,1,mfilename);
if ~ischar(filename) || size(filename, 1)>1
    error(message('bioinfo:soapread:InvalidInput'));
end

% Parse inputs
[blockRead, blockSize, skip, aligndetails] = parse_inputs(varargin{:});

if  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
    fid = fopen(filename, 'r');
    c = onCleanup(@()fclose(fid));
    if blockRead
        if isinf(skip)
            ftext = cell(0,1); % SOAPREAD(...,'BLOCKREAD', Inf)
        else
            if skip>0
                textscan(fid,'%*s',skip,'Delimiter','\n'); % jump records
            end
            if isinf(blockSize)
                ftext = textscan(fid, '%s','Delimiter', '\n'); % read everything
            else
                ftext = textscan(fid,'%s', blockSize,'Delimiter','\n');
            end
            ftext = strtrim(ftext{1});
        end
        if isempty(ftext)
            error(message('bioinfo:soapread:StartTooBig', sprintf('%d',skip + 1)));
        end
    else
        try
            ftext = textscan(fid,'%s','Delimiter','\n');
            ftext = strtrim(ftext{1});
            
        catch theErr
            if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:soapread:FileTooBig'));
            else
                rethrow(theErr);
            end
        end
    end
else  % input must be a string
    %assume that if the string does not contain a delimiter, it is intended
    %to be a filename
    if all(filename~=sprintf('\t'))
        error(message('bioinfo:soapread:FileNotFound', filename));
    end
    ftext = textscan(filename,'%s','Delimiter','\n');
    ftext = strtrim(ftext{1});
end

try
    % Initialize output structure
    if aligndetails
        data(numel(ftext),1) = struct('QueryName', [],'Sequence', [],'Quality',[], 'NumHits',[], ...
            'PairedEndSourceFile',[], 'Length', [], 'Strand', [], 'ReferenceName', [], 'Position', [], ...
            'AlignDetails', []);
    else
        data(numel(ftext),1) = struct('QueryName', [],'Sequence', [],'Quality',[], 'NumHits',[], ...
            'PairedEndSourceFile',[], 'Length', [], 'Strand', [], 'ReferenceName', [], 'Position', []);
    end

    for row = 1:numel(ftext)
        % Read the text and parse each line
        [rtext pos] = textscan(ftext{row}, '%s%s%s%u32%s%u32%s%s%u32',1, 'Delimiter', '\t');
        if isempty(rtext{9}) % when textscan fails to read one of the fields they are returned empty
            error(message('bioinfo:soapread:IncompleteEntry', row))
        else
            % index into cell and extract data
            rtext([1 2 3 5 7 8]) = [rtext{[1 2 3 5 7 8]}];
        end
        if aligndetails
            % The rest of the string in the file is dependent on the HitType
            lastpart = ftext{row}(pos+1:end);
            if isempty(lastpart)
                error(message('bioinfo:soapread:IncompleteEntry', row))
            end
            data(row) = cell2struct([rtext'; {lastpart}], {'QueryName', 'Sequence', 'Quality', 'NumHits',...
                'PairedEndSourceFile', 'Length', 'Strand', 'ReferenceName', 'Position', 'AlignDetails'});
        else
            data(row) = cell2struct(rtext', {'QueryName', 'Sequence', 'Quality', 'NumHits',...
                'PairedEndSourceFile', 'Length', 'Strand', 'ReferenceName', 'Position'});
        end
     end %for row
catch ME
    if (strfind(ME.identifier, 'bioinfo:soapread:IncompleteEntry'))
        rethrow(ME)
    elseif (strfind(ME.identifier, 'MATLAB:nomem'))
        error(message('bioinfo:samread:FileTooBig'));

    else % Generic Error
        error(message('bioinfo:soapread:IncorrectDataFormat'));
    end
end


%--------------------------------------------------------------------------
function [blockRead, blockSize, skip, aligndetails] = parse_inputs(varargin)
% Parse input PV pairs.

% default value
blockRead = false;
aligndetails = true;
blockSize = inf;
skip = 0;

% get input arguments
if rem(nargin,2) ~= 0
    error(message('bioinfo:soapread:IncorrectNumberOfArguments', mfilename));
end
okargs = {'blockread', 'aligndetails'};
for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % range
            if ~isnumeric(pval) || numel(pval)> 2 || ...
                    isempty(pval) || any(pval<1) || any(rem(pval,1))
                error(message('bioinfo:soapread:BadBlockRead'))
            end
            blockRead = true;
            pval = double(pval);
            if numel(pval)==1
                skip = pval-1;
                blockSize = 1;
            else
                skip = min(pval)-1;
                blockSize = max(pval)-skip;
            end
        case 2 %aligndetails
            aligndetails = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end

