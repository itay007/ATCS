function obj = setSubsequence(obj, Y, X, Z)
%SETSUBSEQUENCE set the partial value of the 'Sequence' property.
%
%   SETSUBSEQUENCE(OBJ, Y, X, Z) sets the positions indexed by Z of the
%   elements indexed by X in the 'Sequence' property of a BioSeq (or
%   derived) object OBJ to the value Y. Y must be a cell array of strings
%   with size compatible with the number of elements indexed by X. X must
%   be a vector of positive integers, a logical vector, or a cell array of
%   strings corresponding to valid 'Header' values. Z must be an array of
%   positive integers or a logical array. 
%
%   NOTE: SETSUBSEQUENCE(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETSUBSEQUENCE(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioSeq object.
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Create copies of obj with the first five positions of the second sequence
%   % equal to 'NNNNN'.
%   out1 = setSubsequence(obj, {'NNNNN'}, 2, 1:5)
%   out2 = setSubsequence(obj, {'NNNNN'}, {'H2'}, 1:5)
%   out3 = setSubsequence(obj, {'NNNNN'}, [false true false], 1:5)
%
%   % Set the 10th position of each sequence in obj equal to 'X'.
%   obj = setSubsequence(obj, {'X', 'X', 'X'}, 1:3, 10)
%
%   See also BIOSEQ, BIOSEQ/GETSUBSEQUENCE, BIOSEQ/SET, BIOSEQ/SETHEADER,
%   BIOSEQ/SETSEQUENCE, BIOSEQ/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 4, ['BioSeq:' mfilename])
checkScalarInput(obj);

%=== check that Z is valid
if ~isnumeric(Z) && ~islogical(Z)
    error(message('bioinfo:BioSeq:setSubsequence:InvalidSubsequencePosition'))
end

[X,ex] = uniformizeIndexToNumericArray(obj,X,numel(Y));

if ~iscellstr(Y) 
    error(message('bioinfo:BioSeq:setSubsequence:InvalidSubsequence'));
end
if ~isempty(ex)
	Y = Y(ex); % expand to account for duplicate Headers
end
if numel(Y) ~= numel(X)
    error(message('bioinfo:BioSeq:setSubsequence:IncorrectSubsequenceSize'));
end

%=== modify the subsequences according to the specified values and positions
seqs = obj.getSequence;
for i = 1:numel(X)
	seqs{X(i)}(Z) = Y{i};
end
% Note: If the position(s) specified are greater than the current sequence
%       length(s), the sequence(s) is(are) grown in size. Blank spaces are
%       left between the new positions of the original sequence last position.

%=== set Sequence to the new value 
obj = setSequence(obj, seqs);


