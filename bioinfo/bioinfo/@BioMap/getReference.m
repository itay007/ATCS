function out = getReference(obj, varargin)
%GETREFERENCE retrieve the 'Reference' property of a BioMap object.
%
%   OUT = GETREFERENCE(OBJ) retrieves the 'Reference' property in a BioMap
%   object. OUT is a cell array of strings.
%
%   OUT = GETREFERENCE(OBJ, X) retrieves the 'Reference' property for the
%   elements indexed by X in a BioMap object. X must be a vector
%   of positive integers, a logical vector, or a cell array of strings
%   corresponding to valid 'Header' values. OUT is a cell array of strings.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the 'Reference' value of the second element in the object.
%   getReference(obj, 2)
%   getReference(obj, {s(2).QueryName})
%   getReference(obj, [false true])
%
%   % Retrieve the 'Reference' value of the first and third elements in the
%   % object.
%   getReference(obj, [1 3])
%   getReference(obj, {s([1 3]).QueryName})
%   getReference(obj, [true false true])
%
%   % Retrieve the 'Reference' values of all elements in the object.
%   getReference(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETREFERENCE.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'Reference',varargin{:});

