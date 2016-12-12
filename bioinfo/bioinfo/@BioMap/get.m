function out = get(obj,varargin)
%GET retrieve a BioMap object property.
%
%   GET(OBJ) returns all properties of BioMap object OBJ in a scalar
%   structure, where each field name is a property of OBJ, and each
%   field contains the value of that property.
%
%   GET(OBJ, PROP) returns the value of one or more properties of BioMap
%   object OBJ. PROP is a MATLAB string with the property name or a 1-by-N
%   or N-by-1 cell array of strings containing multiple property
%   names. When the input PROP is a cell, the output is also a cell
%   containing the values for each requested property.
%
%   Examples:
%   % Create a BioMap object from a SAM formatted file:
%   obj = BioMap('ex1.sam')
%
%   % Retrieve the values of the 'Header' property.
%   get(obj, 'Header')
%
%   % Retrieve the values of the 'Header' and 'Sequence' properties.
%   v = get(obj, {'Header', 'Sequence'})
%
%   % Transform the BioMap object OBJ into a structure.
%   str = get(obj)
%
%   See also BIOMAP, BIOMAP/GETSUBSET, BIOREAD.

%   Copyright 2012 The MathWorks, Inc.

out = get@BioSeq(obj,varargin{:});

if nargin>1
    pname = varargin{1};
    if ischar(pname) && strcmp(pname,'Reference')
        if isempty(obj.DictionaryMapping)
            out = obj.SequenceDictionary(out);
        else
            out = obj.SequenceDictionary(obj.DictionaryMapping(out));
        end
    elseif iscellstr(pname)
        h = find(strcmp(pname,'Reference'));
        if ~isempty(h)
            if isempty(obj.DictionaryMapping)
                out{h} = obj.SequenceDictionary(out{h});
            else
                out{h} = obj.SequenceDictionary(obj.DictionaryMapping(out{h}));
            end
        end
    end
else % == out is a structure
    % The upper class get method does not know how to fix the
    % Reference property to have names instead of numeric indices:
    if isempty(obj.DictionaryMapping)
        out.Reference = obj.SequenceDictionary(out.Reference);
    else
        out.Reference = obj.SequenceDictionary(obj.DictionaryMapping(out.Reference));
    end
end