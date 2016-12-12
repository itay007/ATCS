function obj = setHeader(obj,varargin)
%SETHEADER set the 'Header' property of a BioSeq (or derived) object.
%
%   SETHEADER(OBJ, Y) sets the 'Header' property of a BioSeq (or derived)
%   object OBJ to Y. Y  must be a cell array of strings with size equal to
%   the number of elements in the object OBJ.
%
%   SETHEADER(OBJ, Y, X) sets the 'Header' property of the elements indexed
%   by X in a BioSeq object OBJ to the value Y. Y must be a cell array of
%   strings with as many elements as specified by X. X must be a vector of
%   positive integers, a logical vector, or a cell array of strings
%   corresponding to valid 'Header' values.
%
%   NOTE: SETHEADER(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETHEADER(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioSeq (or derived) object. 
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Create copies of obj with the 'Header' value of the second element
%   % equal to 'N'. 
%   out1 = setHeader(obj, {'N'}, 2) 
%   out2 = setHeader(obj, {'N'}, {'H2'}) 
%   out3 = setHeader(obj, {'N'}, [false true false])
%
%   % Reset the 'Header' property.
%   obj = setHeader(obj, {})
%
%   % Set the 'Header' property to new values.
%   obj = setHeader(obj, {'N1', 'N2', 'N3'})
%   
%   See also BIOSEQ, BIOSEQ/GETHEADER, BIOSEQ/SET,
%   BIOSEQ/SETSEQUENCE, BIOSEQ/SETSUBSEQUENCE, BIOSEQ/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Header',true,varargin{:});
