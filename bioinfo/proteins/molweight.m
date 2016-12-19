function mw = molweight(str)
%MOLWEIGHT calculates molecular weight for a protein.
%
%   MOLWEIGHT(SEQUENCE) calculates the molecular weight of amino acid
%   sequence SEQUENCE.
%
%   Example:
%
%       rhodopsin = getgenpept('NP_000530')
%       rhodopsinMW = molweight(rhodopsin)
%
%   See also AACOUNT, ATOMICCOMP, ISOELECTRIC, OLIGOPROP, PROTEINPLOT, PROTEINPROPPLOT.

%   Copyright 2003-2008 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);
if isempty(str)
    mw = 0;
    return
end

% If the input is a structure then extract the Sequence data.
if isstruct(str)
    str = bioinfoprivate.seqfromstruct(str);
end

waterweight = 18.015;

if  ~bioinfoprivate.isaa(str)
    error(message('bioinfo:molweight:NotAASequence'));
end

if ~ischar(str)
    str = int2aa(str);
end
mw = sum(aminoacidinfo(str,'MolecularWeight')) - (waterweight * (length(str) - 1));

