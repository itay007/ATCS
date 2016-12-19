function out = getMatePosition(obj, varargin)
%GETMATEPOSITION retrieve the 'MatePosition' property of a BioMap object.
%
%   OUT = GETMATEPOSITION(OBJ) retrieves the 'MatePosition' property in
%   a BioMap object. OUT is an array of integers.
%
%   OUT = GETMATEPOSITION(OBJ, X) retrieves the 'MatePosition' property
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
%   % Retrieve the 'MatePosition' value of the second element in the object.
%   getMatePosition(obj, 2)
%   getMatePosition(obj, {s(2).QueryName})
%   getMatePosition(obj, [false true])
%
%   % Retrieve the 'MatePosition' value of the first and third elements in the
%   % object.
%   getMatePosition(obj, [1 3])
%   getMatePosition(obj, {s([1 3]).QueryName})
%   getMatePosition(obj, [true false true])
%
%   % Retrieve the 'MatePosition' values of all elements in the object.
%   getMatePosition(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETMATEPOSITION.

%   Copyright 2010-2012 The MathWorks, Inc.

out = getProperty(obj,'MatePosition',varargin{:});
