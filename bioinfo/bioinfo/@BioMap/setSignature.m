function obj = setSignature(obj,varargin)
%SETSIGNATURE set the 'Signature' property of a BioMap object.
%
%   SETSIGNATURE(OBJ, Y) sets the 'Signature' property of a BioMap object OBJ
%   to Y. Y  must be a cell array of strings.
%
%   SETSIGNATURE(OBJ, Y, X) sets the 'Signature' property of the elements
%   indexed by X in a BioMap object OBJ to the value Y. Y must be a cell
%   array of strings. X can be an array of positive integers, a logical
%   array, or a cell array of strings corresponding to valid 'Header'
%   values.
%
%   NOTE: SETSIGNATURE(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETSIGNATURE(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioMap object with the first 10 elements of the structure s.
%   s = samread('ex1.sam');
%   obj = BioMap(s(1:10));
%
%   % Create copies of obj and set the 'Signature' value of the second
%   % element to '*'.
%   out1 = setSignature(obj, {'*'}, 2)
%   out2 = setSignature(obj, {'*'}, {s(2).QueryName})
%   out3 = setSignature(obj, {'*'}, [false true])
%
%   % Set the 'Signature' values of the first and third elements in obj.
%   obj = setSignature(obj, {'*', '*'}, [1 3])
%
%   See also BIOMAP, BIOMAP/GETSIGNATURE, BIOMAP/SET.

%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Signature',true,varargin{:});