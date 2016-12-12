function obj = setMappingQuality(obj,varargin)
%SETMAPPINGQUALITY set the 'MappingQuality' property of a BioMap object.
%
%   SETMAPPINGQUALITY(OBJ, Y) sets the MappingQuality property of a BioMap
%   object OBJ to Y. Y  must be a numeric array.
%
%   SETMAPPINGQUALITY(OBJ, Y, X) sets the 'MappingQuality' property of the
%   elements indexed by X in a BioMap object OBJ to the value Y. Y must be
%   a numeric array. X can be an array of positive integers, a logical
%   array, or a cell array of strings corresponding to valid 'Header'
%   values.
%
%   NOTE: SETMAPPINGQUALITY(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETMAPPINGQUALITY(OBJ,...) to
%   modify OBJ.
%   
%   Examples:
%
%   % Create a BioMap object with the first 10 elements of the structure s.
%   s = samread('ex1.sam');
%   obj = BioMap(s(1:10));
%
%   % Create copies of obj and set the 'MappingQuality' value of the second
%   % element to zero.
%   out1 = setMappingQuality(obj, 0, 2)
%   out2 = setMappingQuality(obj, 0, {s(2).QueryName})
%   out3 = setMappingQuality(obj, 0, [false true])
%
%   % Set the 'MappingQuality' values of the first and third elements in obj.
%   obj = setMappingQuality(obj, [0 0], [1 3])
%
%   % Set all 'MappingQuality' values in the object to 1.
%   y = ones(1, obj.NSeqs);
%   obj = setMappingQuality(obj, y)
%
%   See also BIOMAP, BIOMAP/GETMAPPINGQUALITY, BIOMAP/SET.

%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'MappingQuality',false,varargin{:});