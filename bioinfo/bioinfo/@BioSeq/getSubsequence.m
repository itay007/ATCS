function out = getSubsequence(obj, X, Z)
%GETSUBSEQUENCE retrieve the partial value of the 'Sequence' property.
%
%   OUT = GETSUBSEQUENCE(OBJ, X, Z) retrieves the positions indexed by Z of
%   the 'Sequence' property for the elements indexed by X in a BioSeq (or
%   derived) object. X must be a vector of positive integers, a logical
%   vector, or a cell array of strings corresponding to valid 'Header'
%   values. Z must be an array of positive integers or a logical array. OUT
%   is a cell array of strings.
%
%   Examples:
%
%   % Create a BioSeq (or derived) object.
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Retrieve the first 5 positions of each sequence.
%   getSubsequence(obj, 1:3, 1:5)
%
%   % Retrieve the first five positions of the sequence with header 'H2'.
%   getSubsequence(obj, {'H2'}, 1:5)
%   getSubsequence(obj, [false true false], 1:5)
%
%   See also BIOSEQ, BIOSEQ/GET, BIOSEQ/GETHEADER, BIOSEQ/GETSEQUENCE,
%   BIOSEQ/GETSUBSET, BIOSEQ/SETSUBSEQUENCE.

%   Copyright 2009-2012 The MathWorks, Inc.

bioinfochecknargin(nargin, 3, ['BioSeq:' mfilename])
seqs = getProperty(obj,'Sequence',X);

%=== Return empty if there are no sequences
if obj.getNSeqs == 0
	out = {};
	return
end

%=== check that Z is valid
if ~isnumeric(Z) && ~islogical(Z)
    error(message('bioinfo:BioSeq:getSubsequence:InvalidSubsequencePosition'))
end

if isempty(seqs)
	out = {};
	return
end

N = numel(seqs);
out = cell(N,1);
for i =1:N	
    out{i} = seqs{i}(Z);
end


