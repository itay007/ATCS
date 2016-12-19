function out = getSignature(obj, varargin)
%GETSIGNATURE retrieve the 'Signature' property of a BioMap object.
%
%   OUT = GETSIGNATURE(OBJ) retrieves the 'Signature' property in a BioMap
%   object. OUT is a cell array of strings.
%
%   OUT = GETSIGNATURE(OBJ, X) retrieves the 'Signature' property for the
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
%   % Retrieve the 'Signature' value of the second element in the object.
%   getSignature(obj, 2)
%   getSignature(obj, {s(2).QueryName})
%   getSignature(obj, [false true])
%
%   % Retrieve the 'Signature' value of the first and third elements in the
%   % object.
%   getSignature(obj, [1 3])
%   getSignature(obj, {s([1 3]).QueryName})
%   getSignature(obj, [true false true])
%
%   % Retrieve the 'Signature' values of all elements in the object.
%   getSignature(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETSIGNATURE.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'Signature',varargin{:});

