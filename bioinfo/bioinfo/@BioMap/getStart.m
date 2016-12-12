function out = getStart(obj, varargin)
%GETSTART retrieve the 'Start' property of a BioMap object.
%
%   OUT = GETSTART(OBJ) retrieves the 'Start' property for the elements
%   in a BioMap object. OUT is an array of integers.
%
%   OUT = GETSTART(OBJ, X) retrieves the 'Start' property for the elements
%   indexed by X in a BioMap object. X must be a vector of positive
%   integers, a logical vector, or a cell array of strings corresponding to
%   valid 'Header' values. OUT is an array of integers.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the 'Start' value of the second element in the object.
%   getStart(obj, 2)
%   getStart(obj, {s(2).QueryName})
%   getStart(obj, [false true])
%
%   % Retrieve the 'Start' value of the first and third elements in the
%   % object.
%   getStart(obj, [1 3])
%   getStart(obj, {s([1 3]).QueryName})
%   getStart(obj, [true false true])
%
%   % Retrieve the 'Start' values of all elements in the object.
%   getStart(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETSTART.

%   Copyright 2009-2011 The MathWorks, Inc.

out = getProperty(obj,'Start',varargin{:});