function display(a)
%DISPLAY Display a DataMatrix object.
% 
%   DISPLAY(A) prints the DataMatrix A, including row names and column
%   names.  DISPLAY is called when a semicolon is not used to terminate a
%   statement.
%
%   See also DATAMATRIX, @DATAMATRIX/DISP.

%   Copyright 2008 The MathWorks, Inc. 


isLoose = strcmp(get(0,'FormatSpacing'),'loose');

objectname = inputname(1);
if isempty(objectname)
   objectname = 'ans';
end

if (isLoose), fprintf('\n'); end
fprintf('%s = \n', objectname);
disp(a)
