function b = set(obj, varargin)
%SET set the value of a BioSeq object property.
%
%   SET(OBJ, 'PropertyName', PROPERTYVALUE) sets the specified property
%   'PropertyName' to the specified value PROPERTYVALUE.
%
%   NOTE: SET(OBJ, 'PropertyName', PROPERTYVALUE) without assigning to a
%   variable does not modify OBJ's properties.  Use OBJ = SET(OBJ,
%   'PropertyName', PROPERTYVALUE) to modify OBJ.
%
%   SET(OBJ, 'Property1', Value1, 'Property2', Value2,...) sets
%   multiple property values of a BioSeq object OBJ with a single
%   statement.
%
%   SET(OBJ, 'PropertyName') displays the possible values for the specified
%   property of the BioSeq object OBJ.
%
%   SET(OBJ) displays all properties of the BioSeq OBJ and their possible
%   values.
%
%   STR = SET(OBJ) returns all property names and their possible values for
%   the object OBJ. STR is a structure whose field names are the property
%   names of OBJ, and whose values are cell arrays of possible property
%   values or empty cell arrays.
%
%   Examples:
%
%   % Create a BioSeq object. 
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Set property 'Name' to new value.
%   obj = set(obj, 'Name', 'This is a BioSeq object')
%
%   See also BIOSEQ, BIOSEQ/GET, BIOSEQ/SETHEADER, BIOSEQ/SETSEQUENCE,
%   BIOSEQ/SETSUBSEQUENCE, BIOSEQ/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

checkScalarInput(obj);

%=== Display allowed values if new value is not given (nargin < 3)
if nargin < 3 % set(OBJ) OR set(OBJ, 'propName')
    to = eval(class(obj));
    props = properties(to);
    propVals = struct;
    propDescrs = struct;
    for i = 1:numel(props)
        propVals.(props{i}) = to.(props{i});
        if strcmp(props{i},'NSeqs')
            propDescrs.NSeqs = 'Non negative integer.';
        else
            switch class(to.(props{i}))
                case 'cell'
                    propDescrs.(props{i}) = 'Cell array of strings.';
                case 'char'
                    propDescrs.(props{i}) = 'String.';
                case 'uint8'
                    propDescrs.(props{i}) = 'Array of uint8.';
                case 'uint16'
                    propDescrs.(props{i}) = 'Array of uint16.';
                case 'uint32'
                    propDescrs.(props{i}) = 'Array of uint32.';
            end
        end
    end
 
	%=== set(OBJ, 'PropertyName')
	if nargin == 2
		pname = varargin{1};
		
		if all(strcmp(pname, props)==0)
            error(message('bioinfo:BioSeq:set:InvalidProperty', pname, class(to)))
		else
			if nargout == 1
				b = propVals.(pname);
			else
				disp(propDescrs.(pname));
			end
		end
		
	%=== set(OBJ)
	else
		if nargout == 1
			b = propDescrs;
		else
			disp(propDescrs);
		end
	end
	
elseif mod(nargin,2) == 1 % must be set(OBJ, 'PropertyName1', value1, ...)
	
	%=== copy the oject and set individual properties
	b = obj;
	for k = 1:2:nargin-1
		b.(varargin{k})= varargin{k+1};
	end
else
    error(message('bioinfo:BioSeq:set:WrongNumberArgs'));
end

