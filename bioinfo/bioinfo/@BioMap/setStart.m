function obj = setStart(obj,varargin)
%SETSTART set the 'Start' property of a BioMap object.
%
%   SETSTART(OBJ, Y) sets the 'Start' property of a BioMap object OBJ
%   to Y. Y must be a numeric array.
%
%   SETSTART(OBJ, Y, X) sets the Start' property of the elements indexed by
%   X in a BioMap object OBJ to the value Y. Y must be a numeric array. X
%   can be an array of positive integers, a logical array, or a cell array
%   of strings corresponding to valid 'Header' values.
%
%   NOTE: SETSTART(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETSTART(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioMap object with the first 10 elements of the structure s.
%   s = samread('ex1.sam');
%   obj = BioMap(s(1:10));
%
%   % Create copies of obj and set the 'Start' value of the second element to zero.
%   out1 = setStart(obj, 0, 2)
%   out2 = setStart(obj, 0, {s(2).QueryName})
%   out3 = setStart(obj, 0, [false true])
%
%   % Set the 'Start' values of the first and third elements in obj.
%   obj = setStart(obj, [0 0], [1 3])
%
%   % Set all 'Start' values in the object to 1.
%   y = ones(1, obj.NSeqs);
%   obj = setStart(obj, y)
%
%   See also BIOMAP, BIOMAP/GETSTART, BIOMAP/SET.

%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Start',false,varargin{:});