function a = aaabsorbance(aa)
%AAABSORBANCE calculates the absorbance of an amino acid sequence.
%
%  AAABSORBANCE(AA) returns the absorbance (or optical density) of the
%  amino acid sequence AA.

% Author(s) : Gill S.C., von Hippel P.H.
% Reference : Anal. Biochem. 182:319-326(1989).


%   Copyright 2003-2004 The MathWorks, Inc.


if ~nargin
    a = [];
    return
elseif ~bioinfoprivate.isaa(aa)
    error(message('bioinfo:aaabsorbance:NotAnAminoAcid'));
end

e = aaextinctioncoefficient(aa);

a = e / sum(aamolecularweight(aa));
