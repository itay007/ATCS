function tf = isString(s)
%ISSTRING True for a single string
%
%   TF = ISSTRING(S) returns true if S is a single string, or false
%   oherwise.

%   Copyright 2008 The MathWorks, Inc. 


tf = ischar(s) && ( (isvector(s) && (size(s,1) == 1)) || all(size(s)==0) );
end % isString
