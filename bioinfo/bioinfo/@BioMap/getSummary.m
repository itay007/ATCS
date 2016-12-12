function dsout = getSummary(obj)
%GETSUMMARY Print summary of a BioMap object.
%   GETSUMMARY(OBJ) prints a summary of a BioMap object including the names
%   of the references, the number of sequences mapped to each reference and
%   the genomic range that the sequences cover in each reference.
%
%   DS = GETSUMMARY(OBJ) returns the summary information in a dataset. 
%
%   See also BIOMAP, BIOMAP/GETINFO, BIOMAP/GETSUBSET.

%   Copyright 2012 The MathWorks, Inc.

checkScalarInput(obj);

[gr(:,1), gr(:,2)] = getRange(obj);
if isempty(obj.DictionaryMapping)
    r = obj.Index.getField('Reference');
else
    r = obj.DictionaryMapping(obj.Index.getField('Reference'));
end
c = uint32(accumarray(r(:),1,[numel(obj.SequenceDictionary),1]));

ds = dataset({c,'Number_of_Sequences'},{gr,'Genomic_Range'},'ObsNames',obj.SequenceDictionary(:));

s.Name = obj.Name;
if obj.Index.InMemory
    s.Container_Type = 'Data is in memory.';
else
    s.Container_Type = 'Data is file indexed.';
end
s.Total_Number_of_Sequences = obj.NSeqs;
s.Number_of_References_in_Dictionary = numel(obj.SequenceDictionary);

ds.Properties.UserData = s;

if nargout == 0
    disp('BioMap summary:')
    disp(ds.Properties.UserData)
    disp(ds)
else
    dsout = ds;
end