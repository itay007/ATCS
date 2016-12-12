function obj = setSequence(obj,varargin)
%SETSEQUENCE set the 'Sequence' property of a BioSeq (or derived) object.
%
%   SETSEQUENCE(OBJ, Y) sets the 'Sequence' property of a BioSeq object OBJ
%   to Y. Y  must be a cell array of strings with size equal to the number
%   of sequences in the object OBJ.
%
%   SETSEQUENCE(OBJ, Y, X) sets the 'Sequence' property of the elements
%   indexed by X in a BioSeq object OBJ to the value Y. Y must be a cell
%   array of strings with as many elements as specified by X. X must be a
%   vector of positive integers, a logical vector, or a cell array of
%   strings corresponding to valid 'Header' values.
%
%   NOTE: SETSEQUENCE(OBJ,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ = SETSEQUENCE(OBJ,...) to
%   modify OBJ.
%
%   Examples:
%
%   % Create a BioSeq (or derived) object.
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Set the 'Sequence' property value of the second element to a new value 'NNNN'.
%   out1 = setSequence(obj, {'NNNN'}, 2)
%   out2 = setSequence(obj, {'NNNN'}, {'H2'})
%   out3 = setSequence(obj, {'NNNN'}, [false true false])
%
%   % Set the 'Sequence' property to new values.
%   out = setSequence(obj, {'N1', 'N2', 'N3'})
%   
%   See also BIOSEQ, BIOSEQ/GETSEQUENCE, BIOSEQ/SET, BIOSEQ/SETHEADER,
%   BIOSEQ/SETSUBSEQUENCE, BIOSEQ/SETSUBSET.
  
%   Copyright 2009-2012 The MathWorks, Inc.

obj = setProperty(obj,'Sequence',true,varargin{:});
