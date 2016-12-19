function out = seq2regexp(seq,varargin)
%SEQ2REGEXP converts extended NT or AA symbols into a regular expression.
%
%   SEQ2REGEXP(SEQ) for a cell array calls seq2regexp with individual
%   elements of a cell to convert sequences SEQ to regular expression
%   sequences.
%
%   See also REGEXP, REGEXPI, RESTRICT, SEQWORDCOUNT.

% Copyright 2005 The MathWorks, Inc.


out = cell(size(seq));
for i = 1:numel(seq)
    out{i} = seq2regexp(seq{i},varargin{:});
end
