function out = getHeader(obj, varargin)
%GETHEADER retrieve the 'Header' property of a BioSeq (or derived) object.
%
%   OUT = GETHEADER(OBJ) retrieves the 'Header' property in a BioSeq (or
%   derived) object OBJ. OUT is a cell array of strings.
%
%   OUT = GETHEADER(OBJ, X) retrieves the 'Header' property of elements
%   indexed by X in a BioSeq (or derived) object OBJ. X must be a vector of
%   positive integers, a logical vector, or a cell array of strings
%   corresponding to valid 'Header' values. OUT is a cell array of strings.
%
%   Examples:
%
%   % Create a BioSeq (or derived) object.
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%         'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Retrieve the 'Header' value of the second element in the object.
%   getHeader(obj, 2)
%   getHeader(obj, [false true false])
%
%   % Retrieve the 'Header' value of the first and third elements in the
%   % object.
%   getHeader(obj, [1 3])
%   getHeader(obj, [true false true])
%
%   % Retrieve the 'Header' values of all elements in the object.
%   getHeader(obj)
%   getHeader(obj, 1:3)
%
%   See also BIOSEQ, BIOSEQ/GET, BIOSEQ/GETSEQUENCE, BIOSEQ/GETSUBSEQUENCE,
%   BIOSEQ/GETSUBSET, BIOSEQ/SETHEADER.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'Header',varargin{:});
