function out = getQuality(obj, varargin)
%GETQUALITY retrieve the 'Quality' property of a BioRead (or derived) object.
%
%   GETQUALITY(OBJ) retrieves the 'Quality' property in a BioRead (or
%   derived) object. OUT is a cell array of strings.
%
%   OUT = GETQUALITY(OBJ, X) retrieves the 'Quality' property for the
%   elements indexed by X in a BioRead (or derived) object. X must be a
%   vector of positive integers, a logical vector, or a cell array of
%   strings corresponding to valid 'Header' values. OUT is a cell array of
%   strings.
%
%   Examples:
%
%   % Create a BioRead object.
%   obj = BioRead('fastqfile', 'SRR005164_1_50.fastq')
%
%   % Retrieve the 'Quality' value of the second element in the object.
%   getQuality(obj, 2)
%   getQuality(obj, {'SRR005164.2'})
%   getQuality(obj, [false true])
%
%   % Retrieve the 'Quality' value of the first and third elements in the
%   % object.
%   getQuality(obj, [1 3])
%   getQuality(obj, {'SRR005164.1', 'SRR005164.3'})
%   getQuality(obj, [true false true])
%
%   % Retrieve the 'Quality' values of all elements in the object.
%   quals = getQuality(obj)
%
%   See also BIOREAD, BIOREAD/GET, BIOREAD/GETHEADER, BIOREAD/GETSUBSEQUENCE,
%   BIOREAD/GETSUBSET, BIOREAD/SETSEQUENCE.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'Quality',varargin{:});