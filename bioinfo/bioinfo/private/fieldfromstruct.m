function result = fieldfromstruct(seqStruct,field)
%FIELDFROMSTRUCT extracts 'field' value(s) from a structure
%   
%   FIELDFROMSTRUCT(STRUCTURE) returns the field 'field' from STRUCTURE. If
%    a field with some other capitalization of the name exists, then this
%    field is returned but a warning is given.
%
%   If the input contains a nested structure then all matching nodes are
%   returned.

%   Copyright 2003-2009 The MathWorks, Inc.


result{numel(seqStruct)}='';

for i = 1:numel(seqStruct)
    try
        % if the struct has a (field) then we use it
        result{i} = seqStruct(i).(field);
        % call recursively to find the sequence
        if isstruct(result{i})
            result{i} = char(fieldfromstruct(result{i},field));
        end
        if isempty(result{i})
            result{i} = '';
        end
    catch allExceptions %#ok<NASGU>
        % if not we check for fields with different capitalization
        fields = fieldnames(seqStruct);
        matches = find(strcmpi(fields,field));
        if numel(matches) == 1
            result{i} = fieldfromstruct(seqStruct,fields{matches});
            warning(message('bioinfo:fieldfromstruct:StructSeqCapitalization', field, fields{ matches }));
        else
            result{i} = '';
        end
    end
end
if numel(seqStruct) == 1
    result = result{1};
end
