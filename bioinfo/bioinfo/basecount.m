function bases = basecount(dna,varargin)
%BASECOUNT reports nucleotide counts for a sequence.
%
%   BASECOUNT(SEQ) counts the number of occurrences of each nucleotide in
%   the sequence and returns these numbers in a structure with the fields
%   A, C, G, T(U). Ambiguous nucleotide symbols and gaps (-) are ignored by
%   default. Other unrecognized characteres are also ignored, but a warning
%   message is displayed.
%
%   BASECOUNT(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
%   nucleotide symbols are present. Options are: 'Ignore' skips ambiguous
%   symbols, 'Bundle' counts and bundles them into the Ambiguous field of
%   the output structure, 'Prorate' counts and distributes them
%   proportionally in the appropriate fields for the four standard
%   nucleotide symbols, 'Individual' counts and reports them individually,
%   and 'Warn' ignores them and displays a warning message. Default is
%   'Ignore'.
%
%   BASECOUNT(...,'GAPS',true) adds a field to the output structure
%   with the gap count. Default is false, it ignores gap symbols.
%
%   BASECOUNT(...,'CHART',STYLE) creates a chart showing the relative
%   proportions of the nucleotides. Valid styles are 'Pie' and 'Bar'.
%
%   Example:
%
%       S = getgenbank('M10051')
%       basecount(S)
%
%   See also AACOUNT, BASELOOKUP, CODONCOUNT, CPGISLAND, DIMERCOUNT,
%   NMERCOUNT, NTDENSITY, SEQSTATSDEMO, SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.

pieChart = false;
barChart = false;
countGaps = false;
countAmbiguous = 1; % {'Ignore','Bundle','Individual','Warn','Prorate'}

% If the input is a structure then extract the Sequence data.
if isstruct(dna)
    dna = bioinfoprivate.seqfromstruct(dna);
end

if  nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:basecount:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'chart','gaps','ambiguous'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case 1  % graph
                k = bioinfoprivate.optPartialMatch(pval,{'Pie','Bar'}, okargs{k}, mfilename);
                if k==1
                    pieChart = true;
                else
                    barChart = true;
                end
            case 2 % gaps
                countGaps = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            case 3 % ambiguous
                countAmbiguous = bioinfoprivate.optPartialMatch(pval,{'Ignore','Bundle','Individual','Warn','Prorate'}, okargs{k}, mfilename);
        end
    end
end

    
[~, map] = nt2int('a'); 
maxNtInt = double(max(map));
if ischar(dna)
    try
        dnaint = nt2int(dna,'unknown',maxNtInt+1);
    catch allExceptions  %#ok<NASGU>
        % keep only the symbols that nt2int knows how to treat
        dnaint = nt2int(regexprep(dna,'[^\w-\*\?]','?'),'unknown',maxNtInt+1);
    end
else
    dnaint = dna;
    dnaint(dnaint == 0) = maxNtInt+1;
end

% * A C G T U R Y K M S  W  B  D  H  V  N  -(gap)
% 0 1 2 3 4 4 5 6 7 8 9 10 11 12 13 14 15 16

%'a  b c  d e f g  h i j k l m  n o p q r s t u  v  w  x y z  -'
%[1 11 2 12 0 0 3 13 0 0 7 0 8 15 0 0 0 5 9 4 4 14 10 15 6 0 16]);

seqLength = length(dna);

buckets = zeros(1,maxNtInt+1);

for count = 1:seqLength
    buckets(dnaint(count)) = buckets(dnaint(count)) + 1;
end

bases.A = buckets(nt2int('a'));
bases.C = buckets(nt2int('c'));
bases.G = buckets(nt2int('g'));
bases.T = buckets(nt2int('t'));

% For unknown symbols (bucket 17) a warning is always shown:
if buckets(end)
    if ischar(dna)
        unkn =  unique(dna(dnaint==maxNtInt+1));
        warning(message('bioinfo:basecount:UnknownSymbols', unkn));
    else
        warning(message('bioinfo:basecount:UnknownNumericSymbols'));
    end
end

% For ambiguous symbols (bucket 5 to 16):
switch countAmbiguous % {'Ignore','Bundle','Individual','Warn','Prorate'}
    case 1
        % do nothing
    case 2
        bases.Ambiguous = sum(buckets(5:end-2));
    case 3
        bases.R = buckets(nt2int('r'));
        bases.Y = buckets(nt2int('y'));
        bases.K = buckets(nt2int('k'));
        bases.M = buckets(nt2int('m'));
        bases.S = buckets(nt2int('s'));
        bases.W = buckets(nt2int('w'));
        bases.B = buckets(nt2int('b'));
        bases.D = buckets(nt2int('d'));
        bases.H = buckets(nt2int('h'));
        bases.V = buckets(nt2int('v'));
        bases.N = buckets(nt2int('n'));
    case 4
        buckets_with_amb = find(buckets(5:end-2))+4;
        if ~isempty(buckets_with_amb)
            warning(message('bioinfo:basecount:AmbiguousSymbols', int2nt( buckets_with_amb )));
        end
    case 5
        w = 1./[2;2;2;3;3;3;4];
        bases.A = bases.A + buckets(nt2int('rmwdhvn'))*w;
        bases.C = bases.C + buckets(nt2int('ymsbhvn'))*w;
        bases.G = bases.G + buckets(nt2int('rksbdvn'))*w;
        bases.T = bases.T + buckets(nt2int('ykwbdhn'))*w;
end

% Count gaps when required
if countGaps
   bases.Gaps = buckets(nt2int('-'));
end

if pieChart
    pielabels = fieldnames(bases);
    piecounts = structfun(@(x) x, bases);
    h = piecounts>0;
    pie(piecounts(h),pielabels(h))
elseif barChart
    barlabels = fieldnames(bases);
    barcounts = structfun(@(x) x, bases);
    h = barcounts>0;
    hg = bar(barcounts(h));
    set(get(hg,'parent'),'Xtick',1:sum(h),'XTickLabel',barlabels(h));
end
