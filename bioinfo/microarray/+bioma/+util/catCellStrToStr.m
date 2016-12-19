function str = catCellStrToStr(cStr)
%CATCELLSTRTOSTR Concatenate the elements of a cell string into a string.
% 
% Concatenate the elements of a cell string into a string with ',' in
% between the words. This is useful for warning and error messages.

%   Copyright 2009 The MathWorks, Inc.


str = '';
for i=1:numel(cStr)
    str = [str cStr{i} ', ']; %#ok<AGROW>
end
str = str(1:end-2); % remove extra comma and space from the end
