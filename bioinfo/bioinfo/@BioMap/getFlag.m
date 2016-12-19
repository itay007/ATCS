function out = getFlag(obj, varargin)
%GETFLAG retrieve the 'Flag' property of a BioMap object.
%
%   OUT = GETFLAG(OBJ) retrieves the 'Flag' property  in a BioMap object.
%   OUT is an array of integers. 
%
%   OUT = GETFLAG(OBJ, X) retrieves the 'Flag' property for the elements
%   indexed by X in a BioMap object. X must be a vector of positive
%   integers, a logical vector, or a cell array of strings corresponding to
%   valid 'Header' values. OUT is an array of integers. 
% 
% Note: Each flag must be interpreted bitwise (0 = false, 1 = true) to
% obtain the following information:
%     0x001: the read is paired in sequencing
%     0x002: the read is mapped in a proper pair
%     0x004: the read is unmapped
%     0x008: the mate is unmapped
%     0x010: strand of the read (0 = forward, 1 = reverse)
%     0x020: strand of the mate (0 = forward, 1 = reverse)
%     0x040: the read is the first in a pair
%     0x080: the read is the second in a pair
%     0x100: the alignment is not primary
%     0x200: the read fails platform/vendor quality checks
%     0x400: the read is a PCR or optical duplicate
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the 'Flag' value of the second element in the object.
%   getFlag(obj, 2)
%   getFlag(obj, {s(2).QueryName})
%   getFlag(obj, [false true])
%
%   % Retrieve the 'Flag' value of the first and third elements in the
%   % object.
%   getFlag(obj, [1 3])
%   getFlag(obj, {s([1 3]).QueryName})
%   getFlag(obj, [true false true])
%
%   % Retrieve the 'Flag' values of all elements in the object.
%   getFlag(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETFLAG.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'Flag',varargin{:});
