function obj = setMatePosition(obj,varargin)
%SETMATEPOSITION set the 'MatePosition' property of a BioMap object.
%
%   SETMATEPOSITION(OBJ, Y) sets the 'MatePosition' property of a BioMap
%   object OBJ to Y. Y  must be a numeric array.
%
%   SETMATEPOSITION(OBJ, Y, X) sets the 'MatePosition' property of the
%   elements indexed by X in a BioMap object OBJ to the value Y. Y must be
%   a numeric array. X can be an array of positive integers, a logical
%   array, or a cell array of strings corresponding to valid 'Header'
%   values.
%
%   NOTE: SETMATEPOSITION(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETMATEPOSITION(OBJ,...) to
%   modify OBJ.
%   
%   Examples:
%
%   % Create a BioMap object with the first 10 elements of the structure s.
%   s = samread('ex1.sam');
%   obj = BioMap(s(1:10));
%
%   % Create copies of obj and set the 'MatePosition' value of the second
%   % element to zero.
%   out1 = setMatePosition(obj, 0, 2)
%   out2 = setMatePosition(obj, 0, {s(2).QueryName})
%   out3 = setMatePosition(obj, 0, [false true])
%
%   % Set the 'MatePosition' values of the first and third elements in obj.
%   obj = setMatePosition(obj, [0 0], [1 3])
%
%   % Set all 'MatePosition' values in the object to 0.
%   y = zeros(1, obj.NSeqs);
%   obj = setMatePosition(obj, y)
%
%   See also BIOMAP, BIOMAP/GETMATEPOSITION, BIOMAP/SET.

%   Copyright 2010-2012 The MathWorks, Inc.

obj = setProperty(obj,'MatePosition',false,varargin{:});