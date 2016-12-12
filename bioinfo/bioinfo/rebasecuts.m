function [parts, cSites] = rebasecuts(seq,varargin)
%REBASECUTS Finds restriction enzymes that cut a sequence.
%
%   [ENZYMES, SITES] = REBASECUTS(SEQ) finds all the restriction enzymes
%   that cut the sequence SEQ. ENZYMES is a cell array with the names of
%   restriction enzymes from REBASE. SITES is a vector of cut
%   sites with the base number before every cut relative to the sequence.
%
%   REBASECUTS(SEQ, GROUP) limits the search to those enzymes contained in
%   GROUP. GROUP is a cell array with the names of valid restriction
%   enzymes.
%
%   REBASECUTS(SEQ, [Q,R]) limits the search to those enzymes that cut
%   after base position Q and before base position R relative to the
%   sequence.
%
%   REBASECUTS(SEQ, S) limits the search to those enzymes that cut just
%   after base position S.
%
%
%   REBASE, the Restriction Enzyme Database, is a collection of information
%   about restriction enzymes and related proteins. For more information on
%   REBASE, go to http://rebase.neb.com/rebase/rebase.html.
%
%   Examples:
%
%       seq = 'AGAGGGGTACGCGCTCTGAAAAGCGGGAACCTCGTGGCGCTTTATTAA'
%       % Look for all possible cleavage sites in the sequence seq.
%       [enzymes sites] = rebasecuts(seq)
%       % Find where restriction enzymes CfoI and Tru9I will cut.
%       [enzymes sites] = rebasecuts(seq, {'CfoI','Tru9I'})
%       % Search for any possible enzymes that cut after base 7.
%       enzymes  = rebasecuts(seq, 7)
%       % Get the subset of enzymes that cut between base 11 and 37
%       enzymes  = rebasecuts(seq, [11 37])
%
%   See also CLEAVE, CLEAVELOOKUP, PRIMERDEMO, REGEXP, SEQ2REGEXP,
%   SEQSHOWWORDS, RESTRICT.

%   Copyright 2003-2012 The MathWorks, Inc.

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

if nargin > 2
    error(message('bioinfo:rebasecuts:IncorrectNumberOfArguments', mfilename));
elseif(nargin == 2)
    position = varargin{1};
else
    position = true;
end

if isnumeric(seq)
    seq = int2nt(seq);
elseif ischar(seq)
    seq = upper(seq);   % make sure it is upper case
end

[name, patt,len, ncuts,blunt,c1,c2,c3,c4]=readrebase; %#ok
zymesind = [ ];
cSites = [ ];
% rebase counts go -3 -2 -1 1 2 3 so we have to add 1 to negative values
c1(c1<0) = c1(c1<0)+1; 
c3(c3<0) = c3(c3<0)+1; 

for i = 1:length(patt)
    strt = regexp(seq, seq2regexp(patt{i},'alphabet','nt'));
    if(~isempty(strt))
        zymesind = [zymesind; repmat(i,length(strt),1)]; %#ok<AGROW>
        cSites = [cSites; strt'+c1(i)-1]; %#ok<AGROW>
        if(c3(i) ~= 0)
            zymesind = [zymesind; repmat(i,length(strt),1)]; %#ok<AGROW>
            cSites = [cSites; strt'+c3(i)-1]; %#ok<AGROW>
        end
    end
end
parts = name(zymesind);

% filtering the enzymes depending on their cut sites
if(isnumeric(position)&& isreal(position))
    if(~all(position<length(seq)&position>0))
        error(message('bioinfo:rebasecuts:BadPositionValue'));
    end
    if(length(position)==1)
        ind = (cSites == position);
        parts = parts(ind);
        cSites = cSites(ind);
    elseif(length(position)==2)
        if(position(1) <= position(2))
            between = position(1):(position(2)-1);
            ind = ismember(cSites,between);
            parts = parts(ind);
            cSites = cSites(ind);
        else
            error(message('bioinfo:rebasecuts:BadPositionSet'));
        end
    end
end

% filtering the enzymes depending on the input subset
if(iscell(position))
    exactname = [];
    for i = 1:length(position)% check for spelling case
        found =  strcmpi(position(i),name);
        exactname = [exactname; name(found)]; %#ok<AGROW>
    end
    ind = ismember(parts,exactname);
    parts = parts(ind);
    cSites = cSites(ind);
end

% filtering the enzymes depending on single char input
if(ischar(position))
    found =  strcmpi(position,parts);
    parts = parts(found);
    cSites = cSites(found);
end


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


