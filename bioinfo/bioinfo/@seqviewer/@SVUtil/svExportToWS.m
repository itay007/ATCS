function msg = svExportToWS(variable_name, variable_data, check_variable)
%SVEXPORTTOWS Call by SEQVIEWER Java components to export data to workspace.
%
%   MSG = SVEXPORTTOWS(VARNAME, DATA,CHHKVAR) called by SEQVIEWER Java
%   components to assign data, DATA, to a variable, VARNAME, in the base
%   workspace. CHKVAR is logical flag to check for duplicated variable
%   names in the workspace. MSG contains a string resulting from any
%   error condition, otherwise MSG is empty.

%   Copyright 2005-2012 The MathWorks, Inc.


msg = '';

if numel(variable_name) == 1

    % When exporting only one variable:
    if check_variable
        msg = checkExistingVariableInWorkspace(variable_name);
    end
    if isempty(msg)
        try
            assignin('base', variable_name{1}, variable_data{1});
        catch ME
            if strcmp(ME.identifier,'MATLAB:assigninInvalidVariable')
                msg = 'The Variable Name must be a valid MATLAB variable name.';
            else
                msg = ME.message;
            end
        end
    end
    
else
    
    % When exporting multiple variables (this only happens when exporting
    % translated sequences) 
    uniqueArray = unique(variable_name);
    if numel(variable_name) == numel(uniqueArray)
        % All variables have UNIQUE names
        % Note: we do not need to check for valid MATLAB variable names
        % since this was already done by TextFieldTableCellEditor in java.
        if check_variable
            msg = checkExistingVariableInWorkspace(variable_name);
        end
        if isempty(msg)
            for i = 1: numel(variable_name)
                try
                    assignin('base', variable_name{i}, variable_data{i});
                catch ME
                    msg = ME.message;
                    break;
                end
            end
        end  
    else % Variable do NOT have unique names:
        if numel(variable_name)==3 && numel(uniqueArray)==1
            % If all three variables have the same name an structure is
            % created with proper fiednames:
            tempstruct.Frame1 = variable_data{1};
            tempstruct.Frame2 = variable_data{2};
            tempstruct.Frame3 = variable_data{3};
            if check_variable
                msg = checkExistingVariableInWorkspace(variable_name(1));
            end
            if isempty(msg)
                try
                    assignin('base', variable_name{1}, tempstruct);
                catch ME
                    msg = ME.message;
                end
            end
        else % There are repeated names, but it is not possible to figure  
             % out what frames they belong to.
             msg = 'The Variable Names for each frame must be unique.';
        end
    end
    
end
end


function msg = checkExistingVariableInWorkspace(variable_name)
    %check for names already in the workspace
    n = numel(variable_name);
    flag = false(n,1);
    for i = 1:n
        if evalin('base',['exist(''',variable_name{i},''', ''var'');'])
            flag(i) = true;
        end
    end
    dupnames = unique(variable_name(flag));
    switch numel(dupnames)
        case 0
            msg = '';
        case 1
            queststr = ['"' dupnames{1} '"'...
                ' already exists. Do you want to overwrite it?'];
            msg = [queststr '_exist'];
        case 2
            queststr = ['"' dupnames{1} '" and "' dupnames{2} ...
                '" already exist. Do you want to overwrite them?'];
            msg = [queststr '_exist'];
        case 3
            queststr = [sprintf('"%s" , ', dupnames{1:end-2}), ...
                '"' dupnames{end-1} '" and "' dupnames{end} ...
                '" already exist. Do you want to overwrite them?'];
            msg = [queststr '_exist'];
        otherwise
            msg = 'svExportToWS error';
    end
end




