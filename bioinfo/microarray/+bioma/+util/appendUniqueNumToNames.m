function names = appendUniqueNumToNames(name, nummax, start)
%APPENDUNIQUENUMTONAMES Generate a cell string of unique names by
%  appending unique numbers. 
% 
%   NAMES = APPENDUNIQUENUMTONAMES(NAME, NUMMAX, START) Generate a cell
%   string of unique names by appending unique numbers from START to NUMMAX
%   to NAME. By default, START = 1.

%   Copyright 2009 The MathWorks, Inc. 


if nargin < 3
    start = 1;
end

if start > nummax
    start = 1;
end

names = cellstr([ repmat(name,nummax-start+1,1)  num2str((start:nummax)')]);
end % appendUniqueNumToName.m
