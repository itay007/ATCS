function [parts, cuttingSites, lengths] = restrict(seq,pattern,varargin)
%RESTRICT Cuts a sequence at specified restriction sites.
%
%   FRAGMENTS = RESTRICT(SEQ,ENZYME) cuts a SEQ sequence into fragments at
%   the restriction sites of restriction enzyme ENZYME where ENZYME is the
%   name of a restriction enzyme from REBASE. The return values
%   are stored in a cell array of sequences.
%
%   FRAGMENTS = RESTRICT(SEQ,PATTERN,POSITION) cuts a SEQ sequence into
%   fragments at specified restriction sites specified by nucleotide
%   pattern PATTERN. POSITION defines the position on the PATTERN where the
%   sequence is cut. POSITION 0 corresponds to the 5' end of the PATTERN.
%   PATTERN can be a regular expression.
%
%   [FRAGMENTS, CUTTING_SITES] = RESTRICT(...) returns a numeric vector
%   with the indices representing the cutting sites. A 0 (zero) will be
%   added to the list, so numel(FRAGMENTS)==numel(CUTTING_SITES).
%   CUTTING_SITES+1 can be used to point to the first base of every
%   fragment respective to the original sequence.
%
%   [FRAGMENTS, CUTTING_SITES, LENGTHS] = RESTRICT(...) returns a numeric
%   vector with the lengths of every fragment.
%
%   RESTRICT(...,'PARTIALDIGEST',P) simulates a partial digest where each
%   restriction site in the sequence has a probability P of being cut.
%
%   REBASE, the Restriction Enzyme Database, is a collection of information
%   about restriction enzymes and related proteins. For more information on
%   REBASE, go to http://rebase.neb.com/rebase/rebase.html .
%
%   Examples:
%
%       seq = 'AGAGGGGTACGCGCTCTGAAAAGCGGGAACCTCGTGGCGCTTTATTAA'
%       fragmentsPattern = restrict(seq,'GCGC',3)
%       fragmentsEnzyme = restrict(seq,'HspAI')
%       fragmentsRegExp = restrict(seq,'GCG[^C]',3)
%
%       % Capture the cutting sites and fragment lengths with the fragments
%       [fragments, cut_sites, lengths] = restrict(seq,'HspAI')
%
%   See also CLEAVE, CLEAVELOOKUP, REGEXP, SEQ2REGEXP, SEQSHOWWORDS, REBASECUTS.

%   Copyright 2002-2012 The MathWorks, Inc.


% figure out whether we have a pattern or an enzyme name

if nargin < 2
    error(message('bioinfo:restrict:IncorrectNumberOfArguments', mfilename));
end

returnCuttingSites = false;
returnLengths = false;
if nargout > 1
    returnCuttingSites = true;
end
if nargout > 2
    returnLengths = true;
end

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

cut2 = [];
ratio = 1;
checkOpts = false;
offset = 1;
switch nargin
    case 2  % must be (sequence,enzyme)
        [pattern, position,cut2]=lookupEnzyme(pattern);
    case 3 % must be (sequence,pattern,position)
        position = varargin{1};
    case 4 % must be (sequence,enzyme,varargin_pairs)
        checkOpts = true;
        [pattern, position,cut2]=lookupEnzyme(pattern);
    case 5 % must be (sequence,pattern,position,varargin_pairs)
        checkOpts = true;
        offset = 2;
        position = varargin{1};
    otherwise
        error(message('bioinfo:restrict:IncorrectNumberOfArguments', mfilename));
end

% validate position (it may be corrupted if the user set it)
if ~isnumeric(position)
    error(message('bioinfo:restrict:NonNumericPosition'))
end
if numel(position)>2
    warning(message('bioinfo:restrict:RestrictOnlyTwoCuts'))
    position = position(1:2);
end
if numel(position)==2
    cut2 = max(position);
    position = min(position);
end

if checkOpts
    if rem(nargin,2) == 2-offset
        error(message('bioinfo:restrict:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'partialdigest',''};
    for j=offset:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:restrict:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:restrict:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % ratio
                    ratio = pval;
                    if ~isnumeric(ratio) || ratio <0 || ratio >1
                        error(message('bioinfo:restrict:BadPartialDigestRatio'));
                    end

            end
        end
    end

end

pattern = seq2regexp(pattern,'alphabet','nt');

% deal with the posibility that seq is numeric
if isnumeric(seq)
    seq = int2nt(seq);
elseif ischar(seq)
    seq = upper(seq);   % make sure it is upper case
end

% where are the restriction sites?
sites = regexp(seq,pattern);

position = position - 1;  % regexp finds position of first match

% check that we don't have any overhangs at either end
seqLength = length(seq);
cutSites = sites + position;
if ~isempty(cut2)
    cut2Sites = sites + cut2 -1;
else
    cut2Sites = cutSites;
end
mask = (cutSites < 1) | (cutSites >= seqLength) | (cut2Sites >= seqLength) ;
sites = sites(~mask);

% how many sites are left

numSites = length(sites);

if ratio < 1    % discard some sites if ratio is not 1
    mask = rand(1,numSites)>ratio;
    sites = sites(mask);
    numSites = length(sites);
end

% now for each site we need to split the sequence

% set up outputs
parts = cell(numSites+1,1);
if returnCuttingSites
    cuttingSites = zeros(numSites+1,1);
end
if returnLengths
    lengths = zeros(numSites+1,1);
end

% loop through the
pos = 1;
chunk = 1;
for count = 1:numSites
    chunkEnd = sites(count)+position;
    parts{chunk} =  seq(pos:chunkEnd);
    if returnCuttingSites
        cuttingSites(chunk,1) = pos - 1;
    end
    if returnLengths
        lengths(chunk,1) = chunkEnd - pos + 1;
    end
    pos = chunkEnd + 1;
    chunk = chunk + 1;
    if ~isempty(cut2)
        chunkEnd = sites(count)+cut2-1;
        parts{chunk} =  seq(pos:chunkEnd);
        if returnCuttingSites
            cuttingSites(chunk,1) = pos - 1;
        end
        if returnLengths
            lengths(chunk,1) = chunkEnd - pos + 1;
        end
        pos = chunkEnd + 1;
        chunk = chunk + 1;
    end
end

parts{chunk} = seq(pos:seqLength);
if returnCuttingSites
    cuttingSites(chunk,1) = pos - 1;
end
if returnLengths
    lengths(chunk,1) = seqLength - pos + 1;
end


%-----------------------------------------------------------------------
function [pattern, c1,c3] = lookupEnzyme(theEnzyme)
[name, pattern,  len, ncuts,blunt,c1,c2,c3,c4]=readrebase;%#ok
match = find(strcmpi(name,theEnzyme));

if isempty(match)
    error(message('bioinfo:restrict:CannotFindEnzyme', theEnzyme));
end

pattern = upper(pattern{match});
c1 = c1(match);
if (ncuts(match) > 2 && c3(match)~= c1)
    c3 = c3(match);
else
    c3 = [];
end
% rebase counts go -3 -2 -1 1 2 3 so we have to add 1 to negative values
c1(c1<0) = c1(c1<0)+1; 
c3(c3<0) = c3(c3<0)+1; 

%-----------------------------------------------------------------
% Reference:  REBASE The Restriction Enzyme Database
%
% The Restriction Enzyme data BASE
% A collection of information about restriction enzymes and related
% proteins. It contains published and unpublished references, recognition
% and cleavage sites, isoschizomers, commercial availability, methylation
% sensitivity, crystal and sequence data. DNA methyltransferases, homing
% endonucleases, nicking enzymes, specificity subunits and control proteins
% are also included. Putative DNA methyltransferases and restriction
% enzymes, as predicted from analysis of genomic sequences, are also
% listed. REBASE is updated daily and is constantly expanding.
%
% AUTHORS:
% Dr. Richard J. Roberts and Dana Macelis
%
% LATEST REVIEW:
% Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2007)
% REBASE--enzymes and genes for DNA restriction and modification.
% Nucl. Acids Res. 35: D269-D270. 
%
% OFFICIAL REBASE WEB SITE:
% http://rebase.neb.com
%
