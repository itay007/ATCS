function out = getSequence(obj, varargin)
%GETSEQUENCE retrieve the 'Sequence' property of a BioSeq (or derived) object.
%
%   OUT = GETSEQUENCE(OBJ) retrieves the 'Sequence' property in a BioSeq
%   (or derived) object OBJ. OUT is a cell array of strings.
%
%   OUT = GETSEQUENCE(OBJ, X) retrieves the 'Sequence' property for the
%   elements indexed by X in a BioSeq (or derived) object OBJ. X must be a
%   vector of positive integers, a logical vector, or a cell array of
%   strings corresponding to valid 'Header' values. OUT is a cell array of
%   strings.
%
%   Examples:
%
%   % Create a BioSeq (or derived) object. 
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Retrieve the 'Sequence' value of the second element in the object.
%   getSequence(obj, 2)
%   getSequence(obj, {'H2'})
%   getSequence(obj, [false true false])
%
%   % Retrieve the 'Sequence' value of the first and third elements in the
%   % object.
%   getSequence(obj, [1 3])
%   getSequence(obj, {'H1', 'H3'})
%   getSequence(obj, [true false true])
%
%   % Retrieve the 'Sequence' values of all elements in the object.
%   getSequence(obj)
%   getSequence(obj, 1:3)
%
%   See also BIOSEQ, BIOSEQ/GET, BIOSEQ/GETHEADER, BIOSEQ/GETSUBSEQUENCE,
%   BIOSEQ/GETSUBSET, BIOSEQ/SETSEQUENCE.

%   Copyright 2009-2012 The MathWorks, Inc.

out = getProperty(obj,'Sequence',varargin{:});
