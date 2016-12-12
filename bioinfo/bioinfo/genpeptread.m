function data=genpeptread(gptext,varargin)
%GENPEPTREAD reads GenPept format data files.
%
%   DATA = GENPEPTREAD(FILE) reads in the GenPept formatted sequence from
%   FILE and creates a structure DATA containing fields corresponding to
%   the GenPept keywords. If the file contains information about multiple
%   sequences, then the information will be stored in an array of
%   structures.
%
%   FILE can also be a URL or a MATLAB character array that contains the
%   text of a GenPept format file.
%
%   Based on version 179.0 of GenBank
%
%   Examples:
%
%       % Download a GenPept file to your local drive.
%       getgenpept('AAA59174', 'TOFILE', 'HGENPEPTAAA59174.GPT')
%
%       % Then bring it into a MATLAB sequence.
%       data = genpeptread('HGENPEPTAAA59174.GPT')
%
%   See also FASTAREAD, GENBANKREAD, GETGENPEPT, PDBREAD, SEQVIEWER.

% Copyright 2002-2012 The MathWorks, Inc.


data = genbankread(gptext,varargin{:});
