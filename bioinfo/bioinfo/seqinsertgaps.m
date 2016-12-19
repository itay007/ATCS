function seqout = seqinsertgaps(SEQ,POS,N)
%SEQINSERTGAPS insert gaps into a sequence.
%
%   SEQINSERTGAPS(SEQ,POS) inserts gaps in the sequence SEQ before the
%   symbols indexed by the integers in the vector POS.
%
%   SEQINSERTGAPS(SEQ,GAPPEDSEQ) finds the gap positions in GAPPEDSEQ and
%   inserts them into SEQ.
%
%   SEQINSERTGAPS(SEQ,GAPPEDSEQ,N) indicates the symbol length relation
%   between the sequence SEQ and the gapped sequence GAPPEDSEQ, use N = 1
%   when both sequences have the same alphabet, and use N = 3 when SEQ
%   contains nucleotides representing codons and GAPPEDSEQ contains
%   aminoacids. Default is N = 3.
%
%   Example:
%
%    % Extract the coding region of the neuraminidase protein (NA) of two
%    % strains of the Influenza A virus (H5N1), align them as aminoacids
%    % and copy the gaps to the nucleotide sequences for calculating the
%    % synonymous and nonsynonymous substitutions rates:
%    hk01 = getgenbank('AF509094');
%    vt04 = getgenbank('DQ094287');
%    hk01_cds = featuresparse(hk01,'feature','CDS','Sequence',true);
%    vt04_cds = featuresparse(vt04,'feature','CDS','Sequence',true);
%    [sc,al]=nwalign(nt2aa(hk01_cds),nt2aa(vt04_cds),'extendgap',1);
%    hk01_aligned = seqinsertgaps(hk01_cds,al(1,:))
%    vt04_aligned = seqinsertgaps(vt04_cds,al(3,:))
%    [dn,ds] = dnds(hk01_aligned,vt04_aligned,'verbose',true)
%
%   See also BIRDFLUDEMO, DNDS, DNDSDEMO, DNDSML, INT2AA, INT2NT.

%   Copyright 2006-2008 The MathWorks, Inc.


% default for N
if nargin<3
    N = 3;
end

% If the inputs are structures then extract the sequence data from both:
if isstruct(SEQ)
    SEQ = bioinfoprivate.seqfromstruct(SEQ);
end
if isstruct(POS)
    POS = bioinfoprivate.seqfromstruct(POS);
end

if isnumeric(SEQ)
    error(message('bioinfo:seqinsertgaps:noNumeric'))
end

if ischar(POS)
    POS = find(POS(:)=='-');
elseif isnumeric(POS)
    POS = sort(POS(:)) + (0:numel(POS)-1)';
else
    error(message('bioinfo:seqinsertgaps:posInvalid'))
end

S = numel(SEQ) + numel(POS)*N;

if max(POS)*N-1>=S
    error(message('bioinfo:seqinsertgaps:seqTooShort'))
end

seqout = blanks(S);

for i = 0:N-1
    seqout(POS*N-i) = '-';
end

seqout(seqout==' ') = SEQ;


