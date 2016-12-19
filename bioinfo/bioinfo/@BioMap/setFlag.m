function obj = setFlag(obj,varargin)
%SETFLAG set the 'Flag' property of a BioMap object.
%
%   SETFLAG(OBJ, Y) sets the 'Flag' property of a BioMap object OBJ to Y. Y
%   must be a numeric array.
%
%   SETFLAG(OBJ, Y, X) sets the 'Flag' property of the elements indexed by X
%   in a BioMap object OBJ to the value Y. Y must be a numeric array. X can
%   be an array of positive integers, a logical array, or a cell array of
%   strings corresponding to valid 'Header' values.
%
%   NOTE: SETFLAG(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETFLAG(OBJ,...) to
%   modify OBJ.
%   
%   Examples:
%
%   % Create a BioMap object with the first 10 elements of the structure s.
%   s = samread('ex1.sam');
%   obj = BioMap(s(1:10));
%
%   % Create copies of obj and set the 'Flag' value of the second element to zero.
%   out1 = setFlag(obj, 0, 2)
%   out2 = setFlag(obj, 0, {s(2).QueryName})
%   out3 = setFlag(obj, 0, [false true])
%
%   % Set the first and third values of 'Flag' in obj.
%   obj = setFlag(obj, [0 0], [1 3])
%
%   % Set all 'Flag' values in the object to 1.
%   y = ones(1, obj.NSeqs);
%   obj = setFlag(obj, y)
%
%   See also BIOMAP, BIOMAP/GETFLAG, BIOMAP/SET.

%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Flag',false,varargin{:});