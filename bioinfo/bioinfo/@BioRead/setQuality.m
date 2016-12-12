function obj = setQuality(obj,varargin)
%SETQUALITY set the 'Quality' property of a BioRead (or derived) object.
%
%   SETQUALITY(OBJ, Y) sets the 'Quality' property of a BioRead (or
%   derived) object OBJ to Y. Y  must be a cell array of strings.
%
%   SETQUALITY(OBJ, Y, X) sets the 'Quality' property  of the elements
%   indexed by X in a BioRead (or derived) object OBJ to the value Y. Y
%   must be a cell array of strings with as many elements as specified by
%   X. X must be a vector of positive integers, a logical vector, or a cell
%   array of strings corresponding to valid 'Header' values.
%   
%   NOTE: SETQUALITY(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETQUALITY(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioRead object.
%   obj = BioRead('fastqfile', 'SRR005164_1_50.fastq');
%
%   % Set the 'Quality' property value of the second element to a new value.
%   n(1:length(obj.Quality{2})) = 'N'; 
%   out1 = setQuality(obj, {n}, 2)
%   out2 = setQuality(obj, {n}, {'SRR005164.2'})
%   out3 = setQuality(obj, {n}, [false true false])
%
%   See also BIOREAD, BIOREAD/GETQUALITY, BIOREAD/SETHEADER, BIOREAD/SETSEQUENCE, 
%   BIOREAD/SETSUBSEQUENCE, BIOREAD/SETSUBSET.
  
%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Quality',true,varargin{:});
