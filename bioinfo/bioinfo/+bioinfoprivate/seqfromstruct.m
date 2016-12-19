function seq = seqfromstruct(seqStruct)
%SEQFROMSTRUCT returns Sequence field of a structure.
%
%   SEQFROMSTRUCT(STRUCTURE) returns the field Sequence from STRUCTURE. If
%   no Sequence field exists but a field with some other capitalization of
%   Sequence does exist, then this field is returned but a warning is
%   given.
%
%   If the input contains a nested structure then the leaf node is
%   returned.
%
%   The function gives an MException if it encounters an array of structures.

%   Copyright 2003-2012 The MathWorks, Inc.


if numel(seqStruct) > 1
    msgId = sprintf('bioinfo:seqfromstruct:SeqStructArray');
    x = MException(msgId,getString(message(msgId)));
    x.throwAsCaller;
end
% if the struct has a field Sequence then we use it
try
    seq = seqStruct.Sequence;
catch ME %#ok<NASGU>
    % if not we check for fields with different capitalization
    fields = fieldnames(seqStruct);
    matches = find(strcmpi(fields,'sequence'));
    if numel(matches) == 1
        seq = seqStruct.(fields{matches});
        warning(message('bioinfo:seqfromstruct:StructSeqCapitalization',fields{matches}));
    elseif numel(matches) > 1
        msgId = 'bioinfo:seqfromstruct:StructMultiSeqFields';
        x = MException(msgId,getString(message(msgId)));
        x.throwAsCaller;
    else
        msgId = 'bioinfo:seqfromstruct:StructNoSequence';
        x = MException(msgId,getString(message(msgId)));
        x.throwAsCaller;
    end
end
if isstruct(seq)
    % call recursively to find the sequence
    try
        seq = bioinfoprivate.seqfromstruct(seq);
    catch ME
        if strcmp(ME.identifier,'bioinfo:seqfromstruct:StructNoSequence')
            msgId = 'bioinfo:seqfromstruct:BadSequenceField';
            x = MException(msgId,getString(message(msgId)));
            x.throwAsCaller;
        else
            ME.throwAsCaller;
        end
    end
end
