function [headers, mAlign] = multialignread(filename,varargin)
%MULTIALIGNREAD read a multiple sequence alignment file.
%
%   S = MULTIALIGNREAD(FILENAME) reads a multiple sequence alignment file
%   in which the sequences are determined by blocks with interleaved lines.
%   Common multiple alignment file types, such as ClustalW (.aln), GCG
%   (.msf), and PHYLIP can be read with MULTIALIGNREAD. Every line should
%   start with the sequence header followed by a number (optional, not used
%   by MULTIALIGNREAD) and the space formatted section of the multiple
%   sequence alignment. Sequences are divided into several blocks, same
%   number of blocks must appear for every sequence. The output S is a
%   structure array; where S.Header contains the header information and
%   S.Sequence the amino acid or nucleotide sequences. FILENAME can also be
%   a URL or MATLAB character array that contains the text of the data
%   file.
%
%   [HEADER, SEQS] = MULTIALIGNREAD(FILENAME) reads the file into separate
%   variables HEADER and SEQS.
%
%   MULTIALIGNREAD(...,'IGNOREGAPS', true) removes any gap symbol ('-' or
%   '.') from the sequences. Default is false.
%
%   Example:
%
%       % Reads a multiple alignment of the gag polyprotein of several
%       % HIV strains.
%       gagaa = multialignread('aagag.aln')
%       seqalignviewer(gagaa)
%
%   See also FASTAREAD, GETHMMALIGNMENT, MULTIALIGN, MULTIALIGNWRITE,
%   SEQALIGNVIEWER, SEQCONSENSUS, SEQDISP, SEQPROFILE.

%   Copyright 2003-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

% check input is char
% in a future version we may accept also cells
if ~ischar(filename)
    error(message('bioinfo:multialignread:InvalidInput'))
end

% default
ignoreGaps = false;

% get input arguments
if  nargin > 1
    if rem(nargin,2) ~= 1
        error(message('bioinfo:multialignread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'ignoregaps'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:multialignread:UnknownParameterName', pname));
            %elseif length(k)>1
            %    error('bioinfo:multialignread:faAmbiguousParameterName',...
            %        'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % ignore gaps
                    ignoreGaps = bioinfoprivate.opttf(pval);
                    if isempty(ignoreGaps)
                        error(message('bioinfo:multialignread:IgnoreGapsNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end


try % safety error-catch block
    if size(filename,1)>1  % is padded string
        for i=1:size(filename,1)
            ftext(i,1)=strread(filename(i,:),'%s','whitespace','','delimiter','\n'); %#ok<AGROW>
            ftext{i}(find(~isspace(ftext{i}),1,'last')+1:end)=[]; %#ok<AGROW>
        end
        % try then if it is an url
    elseif (strfind(filename(1:min(10,end)), '://'))
        if (~usejava('jvm'))
            error(message('bioinfo:multialignread:NoJava'))
        end
        try
            ftext = urlread(filename);
        catch allExceptions
            error(message('bioinfo:multialignread:CannotReadURL', filename));
        end
        ftext = textscan(ftext,'%s','delimiter','\n','whitespace','');
        ftext = ftext{1};
        % try then if it is a valid filename
    elseif  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
        fid = fopen(filename,'r');
        ftext = textscan(fid,'%s','delimiter','\n','whitespace','');
        ftext = ftext{1};
        fclose(fid);        
    else  % must be a string with '\n', convert to cell
        ftext = textscan(filename,'%s','delimiter','\n','whitespace','');
        ftext = ftext{1};
    end
    
    %erase empty lines from ftext
    ftext(strcmp('',ftext))=[];
    
    % get the preamble of every line (which should be the header)
    preamble = regexp(ftext,'\S+','once','match');
    
    % get a unique set of keys and indexes to keys for every entry in preamble
    [keys,i,h] = unique(preamble);
    
    % frequency of header occurrences (this gives a pattern of the potential
    % number of blocks)
    freqh = accumarray(h,1);
    
    % keys must contain at least one word character in order to be
    % considered potential headers, otherwise they could just be parts of
    % the consensus row.
    freqh(cellfun(@isempty,regexp(keys,'\w','once'))) = 1;
    
    % try to guess what is the potential number of blocks by detecting the mode
    % again, but ignoring all the headers which freqh==1 (i.e. appeared once)
    
    if all(freqh==1)
        modeFreqh = 1;
    else
        freqFreqh = accumarray(freqh(freqh~=1),1);
        [~,modeFreqh]=max(freqFreqh);
    end
    
    % logical index to the selected keys
    g = (modeFreqh == freqh);
    
    % h -> key index for every line in ftext (or every preamble)
    % i -> preamble index for every key
    % g -> logical index indicating the valid keys
    % g(h) -> logical index indicating the valid lines
    % i(g) -> index to preamble for the valid keys
    
    numSeqs = sum(g);
    numLines = sum(g(h));
    numBlocks = numLines/numSeqs;
    
    % if there are two or more valid keys (ie sequences) AND two or more blocks
    % AND the number of valid lines is a multiple of the number of blocks (or
    % valid keys) AND lines are ordered one after the other for every block
    % THEN we can assume it is a multiple alignment with several blocks
    if (numSeqs>1)  &&  (numBlocks>1)  &&  (rem(numBlocks,1)==0) && ...
            all(all(diff(reshape(find(g(h)),numSeqs,numBlocks))==1))
        recognizedFormat = true;
        ftext = ftext(g(h));
        headers = preamble(sort(i(g)));
    else
        recognizedFormat = false;
    end
    
    % if format has not been recognized, then try one block with several
    % sequences, the lines of interest will be those that have the same length
    if ~recognizedFormat
        % check to see if it is PHYLIP format -- which has number of
        % sequence and length values on the first line of the file and
        % headers only on the first line of each sequence.
        isPhylip = false;
        try
            phylipVals = sscanf(ftext{1},'%d');
            if numel(phylipVals) == 2
                phylipSeq = [ftext{2:phylipVals(1):end}];
                [~, phylipSeq] = strtok(phylipSeq,' ');
                phylipSeq = strrep(phylipSeq,' ','');
                if numel(phylipSeq) == phylipVals(2)
                    numSeqs = phylipVals(1);
                    headers = preamble(2:phylipVals(1)+1);
                    ftext(1)=[];
                    numLines = numel(ftext);
                    isPhylip = true;
                end      
            end
        catch allExceptions %#ok<NASGU>
            % Do nothing
        end
        if ~isPhylip
            lengths = cellfun('length',ftext);
            freqLengths = accumarray(lengths,1);
            [~,modeFreqLengths] = max(freqLengths);
            g = (lengths==modeFreqLengths);
            
            % check that the candidate sequence row are all equally spaces,
            % otherwise this most probably is not a multiple alignment
            if numel(unique(diff(find(g))))>1
                error(message('bioinfo:multialignread:NotEqualSpaces'))
            end
            
            % the selected lines must be contiguous, select the largest region of g
            % with ones
            cg = double(g);
            for j = 2:numel(cg)
                cg(j) = (cg(j-1)+cg(j))*cg(j);
            end
            [ma,lo]=max(cg);
            
            if ma<2
                error(message('bioinfo:multialignread:NotEnoughSequences'));
            end
            
            numSeqs = ma;
            numLines = ma;
            
            % clean input data accordingly with the selected valid lines
            ftext = ftext(lo-ma+1:lo);
            headers = preamble(lo-ma+1:lo);
            
            % in one block formats the top line with number can mess-up the
            % identification of valid lines, since it will also have the same length
            
            toKeep = cellfun('isempty',regexp(headers,'^\s*[0-9]+'));
            
            if any(~toKeep)
                ftext = ftext(toKeep);
                headers = headers(toKeep);
                numSeqs = sum(toKeep);
                numLines = sum(toKeep);
            end
        end
    end
    
    % remove the maximum width of the headers to every line at the beginning
    le = max(cellfun('length',headers));
    for j = 1:numLines
        ftext{j} = ftext{j}(le+1:end); %#ok<AGROW>
    end
    
    % remove all common spaces at the beginning
    le = min(cellfun('length',regexp(ftext,'\s*','match','once')));
    for j = 1:numLines
        ftext{j} = ftext{j}(le+1:end); %#ok<AGROW>
    end
    
    % some clustalw and plasmodb have numbers at the end of some lines
    ftext = regexprep(ftext,'[0-9]+$','');
    
    % some formats have a number here, remove it also
    le = max(cellfun('length',regexp(ftext,'\(*\s*[0-9]+\)*','match','once')));
    for j = 1:numLines
        ftext{j} = ftext{j}(le+1:end);
    end
    
    % after removing numbers there may be some remaining spaces, so we'll pad
    % with blanks and later if the blanks are common to all sequences they'll
    % be removed. We cannot just remove blanks because some formats use them to
    % mark gaps of the alignment
    padTo = max(cellfun('length',ftext));
    for j = 1:numLines
        ftext{j}(end+1:padTo) = ' ';
    end
    
    % concatenate every sequence into a single line
    mAlign = cell(1,numSeqs);
    for j = 1:numSeqs
        mAlign{j} = cat(2,ftext{j:numSeqs:numLines});
    end
    
    % remove all common spaces in every column and fill with trailing spaces
    temp = char(mAlign);
    mAlign = cellstr(temp(:,~all(temp==' ')));
    
    % check that all sequences at least have one symbol (this will catch the
    % conservation markup line)
    toKeep = ~cellfun('isempty',regexp(mAlign,'[A-Za-z]','once'));
    mAlign = mAlign(toKeep);
    headers = headers(toKeep);
    numSeqs = sum(toKeep);
    
    if isempty(headers) || isempty(mAlign)
        error(message('bioinfo:multialignread:EmptyCells'));
    end
    
    if ignoreGaps
        mAlign = regexprep(mAlign,'[-\.~\s]','');
    else % only uniformize gaps (we only allow - or .)
        mAlign = regexprep(mAlign,'[~\s]','-');
    end
    
    if sum(cellfun('length',regexp(mAlign,'[^A-Za-z\.\s*-]')))
        error(message('bioinfo:multialignread:InvalidCharacters')),
    end
    
    % in case of one output put everything into one structure (in headers)
    if nargout ~= 2
        [H(1:numSeqs).Header] = deal(headers{:});
        [H(1:numSeqs).Sequence] = deal(mAlign{:});
        headers = H;
    end
    
    if numel(headers)==0
        error(message('bioinfo:multialignread:ZeroSequencesFound'));
    end
    
catch err
    if ~strfind(lower(err.identifier),'bioinfo:multialignread')
        rethrow(err);
    end
    error(message('bioinfo:multialignread:IncorrectDataFormat'))
end
