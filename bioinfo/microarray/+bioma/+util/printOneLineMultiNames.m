function outStr = printOneLineMultiNames(names, colPad, dimNums, clsName)
%PRINTONELINEMULTINAMES Returns a string of names to be print on one line.
% 
%   STR = PRINTONELINEMULTINAMES(NAMES, PADDINGS, DIMNAMES, CLSNAME)
%   returns a string of NAMES to be print on one line. If NAMES has more
%   than 3 names, use .... PADDINGS is column indentation. DIMNUMS is 1x2
%   vector of dimension numbers, CLSNAME is the class name for empty object.

% Copyright 2009 The MathWorks, Inc.


if (dimNums(1) > 0)
    if isempty(names)
        outStr = [colPad  'NA'];
    else
        outStr = [colPad names{1}];
        if dimNums(1) >= 2
            outStr = [outStr  ', ' names{2}];
        end
        
        if dimNums(1) > 2
            if dimNums(1) == 3
                outStr = [outStr ', ' names{dimNums(1)}];
            else
                outStr = [outStr  ', ...,' names{dimNums(1)}];
                outStr = [outStr ' (' num2str(dimNums(1)) ' total)'];
            end
        end
    end
    
else
    outStr = sprintf('[empty %d-by-%d %s]', dimNums(1), dimNums(2), clsName);
end

end
