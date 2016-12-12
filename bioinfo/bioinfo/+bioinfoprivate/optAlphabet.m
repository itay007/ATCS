function [isAmino, isRNA] = optAlphabet(pval, okarg, matlabfile)
%OPTALPHABET determines the type of alphabet from the input options.
%
%   [ISAMINO, ISRNA] = OPTALPHABET(PVAL, OKARG, MATLABFILE) evaluates the
%   value PVAL of property OKARG for a known type of alphabet and errors if
%   the intended value is invalid. Valid alphabets are 'aa','nt','rna','dna'.

% Copyright 2009-2012 The MathWorks, Inc.


if ~ischar(pval)||~isvector(pval)
    pval = [];
end

isAmino = false;
isRNA = false;

k = find(strncmpi(pval,{'aa','nt','dna','rna','amino','protein'},numel(pval)));

if k
    switch k(1) % The only ambiguety occurs when both options are aa
        case 2
            % nop
        case 3
            % nop
        case 4
            isRNA = true;
        case {1,5,6}
            isAmino = true;
    end
else
    okarg(1) = upper(okarg(1));
    msg = getString(message('bioinfo:optAlphabet:UnknownAlphabet',upper(okarg),upper(matlabfile)));
    msgId = sprintf('bioinfo:%s:%sUnknownAlphabet',matlabfile,okarg);
    x = MException(msgId, msg);
    x.throwAsCaller;
end
 
end % optAlphabet function
