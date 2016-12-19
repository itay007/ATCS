function varargout = getgenpept(accessnum,varargin)
%GETGENPEPT Retrieves sequence information from the NCBI GenPept database.
%
%   GPOUT = GETGENPEPT(ACCESSNUM) searches for the accession number in the
%   GenPept database and returns a structure containing information for
%   the sequence.
%
%   GPOUT = GETGENBANK(...,'PARTIALSEQ',SEQPARAMS) retrieves specified
%   subsequence for requested GenPept file.  SEQPARAMS is a two-element
%   array of integers containing start and end positions for the
%   subsequence [START_AA,END_AA]. START_AA must be an integer between 1
%   and END_AA. END_AA must be an integer between START_AA and length of
%   sequence.
%   
%   GPOUT = GETGENPEPT(...,'TOFILE',FILENAME) saves the data returned from
%   the NCBI database in the file FILENAME using the GenPept format.
%
%   FASTA = GETGENPEPT(...,'FILEFORMAT','FASTA') retrieves from NCBI a
%   FASTA data. The output structure FASTA contains only the fields
%   'Header' and 'Sequence'.
%
%   SEQ = GETGENPEPT(...,'SEQUENCEONLY',true) returns only the sequence
%   for the GenPept entry as a character array. When the SEQUENCEONLY and
%   TOFILE options are used together, the output file is in the Fasta
%   format. 
%
%   GETGENPEPT(...) without output arguments and without specifying the
%   option TOFILE displays information with hyperlinks to the URLS used to
%   search for and retrieve the data.
%
%   If an error occurs while retrieving the GenPept-formatted information,
%   try to run the query again at a later time.  Errors can occur due to
%   internet connectivity issues that are unrelated to the GenPept record.
%
%   Example:
%
%       S = getgenpept('AAA59174')
%
%       % Retrieve Furin-like repeats domain from AAA59174
%       FU = getgenpept('AAA59174','partialseq',[234,281]);
%
%   This retrieves the sequence for the human insulin receptor and stores
%   it in structure S.
%
%   See http://www.ncbi.nlm.nih.gov/About/disclaimer.html for more
%   information on GenPept data.
%
%   See also GENPEPTREAD, GETEMBL, GETGENBANK, GETPDB, SEQSTATSDEMO,
%   SEQVIEWER.

% Copyright 2002-2012 The MathWorks, Inc.


for n = 1:2:length(varargin);
    arg = lower(varargin{n});
    if find(strncmpi(arg,'fileformat',numel(arg)))
        pval = varargin{n+1};
        if find(strncmpi(pval,'GenBank',numel(pval)))
            error(message('bioinfo:getgenpept:BadProteinFileFormat'));
        end
    end
end

% Try to retrieve a record twice in order to circumvent problems due to
% instability of internet connectivity
try
    [varargout{1:nargout}] = getncbidata(accessnum,'database','protein','fileformat','GenPept',varargin{:});
catch ME %#ok<NASGU>
    pause(2); %Pausing for a second before trying query again
    [varargout{1:nargout}] = getncbidata(accessnum,'database','protein','fileformat','GenPept',varargin{:});
end

