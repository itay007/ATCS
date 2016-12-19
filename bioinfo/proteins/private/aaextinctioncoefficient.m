function e = aaextinctioncoefficient(aa)
%AAEXTINCTIONCOEFFICIENT calculates the extinction coeff. of amino acid sequence.
%
%  AAEXTINCTIONCOEFFICIENT(AA) returns the molar extinction coefficient of the
%  amino acid sequence AA.

% Author(s) : Gill S.C., von Hippel P.H.
% Reference : Anal. Biochem. 182:319-326(1989).


%   Copyright 2003-2004 The MathWorks, Inc.

if ~nargin
    e = [];
    return
elseif ~bioinfoprivate.isaa(aa)
    error(message('bioinfo:aaextinctioncoefficient:NotAnAminoAcid'));
end


num_cys = length(strfind(aa,'c'));
num_trp = length(strfind(aa,'w'));
num_tyr = length(strfind(aa,'y'));

e_cys = 120;
e_trp = 5690;
e_tyr = 1280;

e = e_cys * num_cys + e_trp * num_trp + e_tyr * num_tyr;
