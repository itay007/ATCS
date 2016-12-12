function seq = cleansequence(seq)
%CLEANSEQUENCE remove nonletter characters, make sequence lowercase

%   Copyright 2003-2004 The MathWorks, Inc.


% remove the nonletter characters from the sequence
remove = ~isletter(seq);
seq(remove) = '';
seq = lower(seq);

