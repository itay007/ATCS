function tf = isColon(indices)
%ISCOLON Check if a set of indices is ':'.

%   Copyright 2008 The MathWorks, Inc. 


% Check ischar first.  isequal(58,':') alone is true,
% and strcmp({':'},':') is also true
tf = ischar(indices) && strcmp(indices,':');
