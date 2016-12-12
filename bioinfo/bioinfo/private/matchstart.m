function tf = matchstart(string,pattern)
%MATCHSTART matches start of string with pattern, ignoring spaces

% Copyright 2003-2004 The MathWorks, Inc.


tf = ~isempty(regexp(string,['^(\s)*?',pattern],'once'));
