function obj = setReference(obj,varargin)
%SETREFERENCE set the 'Reference' property of a BioMap object.
%
%   SETREFERENCE(OBJ, Y) sets the 'Reference' property of a BioMap object OBJ
%   to Y. Y  must be a cell array of strings.
%
%   SETREFERENCE(OBJ, Y, X) sets the 'Reference' property of the elements
%   indexed by X in a BioMap object OBJ to the value Y. Y must be a cell
%   array of strings. X can be an array of positive integers, a logical
%   array, or a cell array of strings corresponding to valid 'Header'
%   values.
%
%   NOTE: SETREFERENCE(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETREFERENCE(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioMap object with the first 10 elements of the structure s.
%   s = samread('ex1.sam');
%   obj = BioMap(s(1:10));
%
%   % Create copies of obj and set the 'Reference' value of the second
%   % element to '*'.
%   out1 = setReference(obj, {'*'}, 2)
%   out2 = setReference(obj, {'*'}, {s(2).QueryName})
%   out3 = setReference(obj, {'*'}, [false true])
%
%   % Set the 'Reference' values of the first and third elements in obj.
%   obj = setReference(obj, {'*', '*'}, [1 3])
%
%   See also BIOMAP, BIOMAP/GETREFERENCE, BIOMAP/SET.

%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Reference',true,varargin{:});