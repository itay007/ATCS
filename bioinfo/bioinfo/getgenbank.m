function varargout = getgenbank(accessnum,varargin)
%GETGENBANK retrieves sequence information from the NCBI GenBank database.
%
%   GBOUT = GETGENBANK(ACCESSNUM) searches for the accession number in the
%   GenBank database and returns a structure containing information for
%   the sequence.
%
%   GBOUT = GETGENBANK(...,'PARTIALSEQ',SEQPARAMS) retrieves specified
%   subsequence for selected GenBank file.  SEQPARAMS is a two-element
%   array of integers containing the start and end positions for the
%   subsequence [START_BP,END_BP]. START_BP is an integer between 1 and
%   END_BP and END_BP is an integer between START_BP and length of
%   sequence.
%
%   GBOUT = GETGENBANK(...,'TOFILE',FILENAME) saves the data returned from
%   the NCBI database in the file FILENAME using the GenBank format.
%
%   FASTA = GETGENBANK(...,'FILEFORMAT','FASTA') retrieves from NCBI a
%   FASTA data. The output structure FASTA contains only the fields
%   'Header' and 'Sequence'.
%
%   SEQ = GETGENBANK(...,'SEQUENCEONLY',true) returns only the sequence
%   for the GenBank entry as a character array. When the SEQUENCEONLY and
%   TOFILE options are used together, the output file is in the Fasta
%   format. 
%
%   GETGENBANK(...) without output arguments and without specifying the
%   option TOFILE displays information with hyperlinks to the URLS used to
%   search for and retrieve the data.
%
%   If an error occurs while retrieving the GenBank-formatted information,
%   try to run the query again at a later time.  Errors can occur due to
%   internet connectivity issues that are unrelated to the GenBank record.
%
%   Example:
%
%       S = getgenbank('M10051')
%
%       % Retrieve just the coding sequence from M10051
%       CDS = getgenbank('M10051','partialseq',[139,4287]);
%
%   This retrieves the sequence from chromosome 19 that codes for the
%   human insulin receptor and stores it in structure S.
%
%   See http://www.ncbi.nlm.nih.gov/About/disclaimer.html for more
%   information on GenBank data.
%
%   See also GENBANKREAD, GETEMBL, GETGENPEPT, GETPDB, SEQSTATSDEMO,
%   SEQVIEWER.

% Copyright 2002-2012 The MathWorks, Inc.


for n = 1:2:length(varargin)
    arg = varargin{n};
    if find(strncmpi(arg,'fileformat',numel(arg)))
        pval = varargin{n+1};
        if find(strncmpi(pval,'GenPept',numel(pval)))
            error(message('bioinfo:getgenbank:BadNucleotideFileFormat'));
        end
    end
end

% Try to retrieve a record twice in order to circumvent problems due to
% instability of internet connectivity
try
    [varargout{1:nargout}] = getncbidata(accessnum,'fileformat','GenBank','database','nucleotide',varargin{:});
catch ME %#ok<NASGU>
    pause(2); %Pausing for a second before trying query again
    [varargout{1:nargout}] = getncbidata(accessnum,'fileformat','GenBank','database','nucleotide',varargin{:});
end
