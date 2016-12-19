function out = getMappingQuality(obj, varargin)
%GETMAPPINGQUALITY retrieve the 'MappingQuality' property of a BioMap object.
%
%   OUT = GETMAPPINGQUALITY(OBJ) retrieves the 'MappingQuality' property in
%   a BioMap object. OUT is an array of integers.
%
%   OUT = GETMAPPINGQUALITY(OBJ, X) retrieves the 'MappingQuality' property
%   for the elements indexed by X in a BioMap object. X must be a vector
%   of positive integers, a logical vector, or a cell array of strings
%   corresponding to valid 'Header' values. OUT is an array of integers.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the 'MappingQuality' value of the second element in the object.
%   getMappingQuality(obj, 2)
%   getMappingQuality(obj, {s(2).QueryName})
%   getMappingQuality(obj, [false true])
%
%   % Retrieve the 'MappingQuality' value of the first and third elements in the
%   % object.
%   getMappingQuality(obj, [1 3])
%   getMappingQuality(obj, {s([1 3]).QueryName})
%   getMappingQuality(obj, [true false true])
%
%   % Retrieve the 'MappingQuality' values of all elements in the object.
%   getMappingQuality(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETMAPPINGQUALITY.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'MappingQuality',varargin{:});
