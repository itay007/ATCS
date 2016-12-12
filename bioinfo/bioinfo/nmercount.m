function [outputCell,buckets] = nmercount(seq,len,c)
%NMERCOUNT counts n-mers in a sequence.
%
%   NMERCOUNT(SEQ,LEN) counts the number of n-mers, or patterns, of
%   length LEN in sequence SEQ.
%
%   NMERCOUNT(SEQ,LEN,C) only returns n-nmers with cardinality at least C.
%
%   Examples:
%
%       % Load the human insulin receptor precursor from GenPept
%       S = getgenpept('AAA59174','SequenceOnly',true)
%
%       % Counts all 4-mers and display the first 20
%       nmers = nmercount(S,4);
%       nmers(1:20,:)
%
%       % Counts all 3-mers with cardinality at least 3
%       nmercount(S,3,3)
%
%   See also BASECOUNT, CODONCOUNT, DIMERCOUNT.

%   Copyright 2002-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename);

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

% figure out how much big everything is.
seqLen = length(seq);
numberOfNmers = seqLen - len + 1;

% now stash all of the nmers in the array
for count = len:-1:1
    nmers(:,count) = seq(count:count+numberOfNmers-1); %#ok<AGROW>
end

% look for the unique nmers
[nmers,i,j] = unique(nmers,'rows'); 

% now use j to calculate multiplicity
buckets = accumarray(j,1);

% ignore nmers with cardinality less than C
if nargin>2
    h = buckets>=c;
    buckets = buckets(h);
    nmers = nmers(h,1:len);
end

% we want to order by the multiplicity
[buckets,perm] = sort(buckets,1,'descend');
nmers = nmers(perm,:);

% with one output we save this into a cell array.
if nargout <= 1
    outputSize = numel(buckets);
    outputCell = cell(outputSize,2);
    for count = 1:outputSize
        outputCell{count,1} = nmers(count,:);
        outputCell{count,2} = buckets(count);
    end
else % for two outputs the first output has only the nmers
    outputCell = nmers;
end
