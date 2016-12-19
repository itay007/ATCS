function disp(a)
%DISP Display DataMatrix object.
% 
%   DISP(A) displays the DataMatrix A, including row names and column
%   names, without printing the DataMatrix object name. In all other ways
%   it's the same as leaving the semicolon off an expression.
%
%   See also DATAMATRIX, @DATAMATRIX/DISPLAY, FORMAT.

%   Copyright 2008 The MathWorks, Inc.


isLoose = strcmp(get(0,'FormatSpacing'),'loose');
isLong = ~isempty(strfind(get(0,'Format'),'long'));
doubleDigits = 5 + 10*isLong; % 5 or 15
singleDigits = 5 + 2*isLong; % 5 or 7
maxWidth = matlab.desktop.commandwindow.size; % 103 17
maxWidth = maxWidth(1);

printNewLine(isLoose);

if (a.NRows > 0) && (a.NCols > 0)
    colPad = repmat(' ', a.NRows+1, 4);
    if isempty(a.RowNames)
        amChars = [colPad  vCat(' ', num2str((1:a.NRows)'))];
    else
        amChars = [colPad  vCat(' ', char(a.RowNames))];
    end
    
    if isempty(a.ColNames)
        colNames = cellstr(num2str((1:a.NCols)'));
    else
        colNames = a.ColNames;
    end
    
    for i = 1:a.NCols
        col = a.Matrix(:, i);

        if isa(col, 'double')
            colChars = num2str(col, doubleDigits);
        elseif isa(col, 'single')
            colChars = num2str(col, singleDigits);
        else % integer
            colChars = num2str(col);
        end

        colChars = vCat(colNames{i}, colChars);
        if size(amChars,2)+size(colPad,2) + size(colChars,2) > maxWidth
            disp(amChars)
            fprintf('\n');
            printNewLine(isLoose)
            
            if isempty(a.RowNames)
                amChars = [colPad  vCat(' ', num2str((1:a.NRows)'))];
            else
                amChars = [colPad vCat(' ', char(a.RowNames))];
            end
        end

        amChars = [amChars colPad colChars];
    end %for
else
    amChars = sprintf('[empty %d-by-%d DataMatrix]', a.NRows, a.NCols);
end
disp(amChars)
printNewLine(isLoose)
end % disp method

function printNewLine(isLoose)
% Print new line if format spacing is loose
    if isLoose
        fprintf('\n');
    end
end

function chars = vCat(str1, str2)
% Vertically concatenate strings
    chars = char({str1;str2});
end
