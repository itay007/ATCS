function aminoAcids = aacount(peptide,varargin)
%AACOUNT reports amino acid counts for a sequence.
%
%   AACOUNT(SEQ) counts the number of occurrences of each amino acid in the
%   sequence and returns these numbers in a structure. Ambiguous amino acid
%   symbols (BZX), gaps (-), and end terminators (*) are ignored by default.
%   Other unrecognized characteres are also ignored, but a warning message
%   is displayed.
%
%   AACOUNT(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
%   amino acid symbols are present. Options are: 'Ignore' skips ambiguous
%   symbols, 'Bundle' counts and bundles them into the Ambiguous field of
%   the output structure, 'Prorate' counts and distributes them
%   proportionally in the appropriate fields for the standard amino acid
%   symbols, 'Individual' counts and reports them individually, and 'Warn'
%   ignores them and displays a warning message. Default is 'Ignore'.
%
%   AACOUNT(...,'GAPS',true) adds a field to the output structure with the
%   gap count. Default is false, it ignores gap symbols.
%
%   AACOUNT(...,'CHART',STYLE) creates a chart showing the relative
%   proportions of the amino acids. Valid styles are 'Pie' and 'Bar'.
%
%   Example:
%
%       S = getgenpept('AAA59174')
%       aacount(S)
%
%   See also AMINOLOOKUP, ATOMICCOMP, BASECOUNT, CODONCOUNT, DIMERCOUNT,
%   ISOELECTRIC, MOLWEIGHT, PROTEINPLOT, PROTEINPROPPLOT, SEQSTATSDEMO,
%   SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.

pieChart = false;
barChart = false;
countGaps = false;
countAmbiguous = 1; % {'Ignore','Bundle','Individual','Warn','Prorate'}

% If the input is a structure then extract the Sequence data.
if isstruct(peptide)
    peptide = bioinfoprivate.seqfromstruct(peptide);
end

if  nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:aacount:IncorrectNumberOfArguments', mfilename));
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

[~, map] = aa2int('a'); 
maxPeptideInt = double(max(map));
if ischar(peptide)
    try
        peptideint = aa2int(peptide,'unknown',maxPeptideInt+1);
    catch allExceptions %#ok<NASGU>
        % keep only the symbols that aa2int knows how to treat
        peptideint = aa2int(regexprep(dna,'[^\w-\*\?]','?'),'unknown',maxNtInt+1);
    end
else
    peptideint = peptide;
    peptideint(peptideint == 0) = maxPeptideInt+1;
end

%
%   A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *  -  ?
%   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 0

seqLength = length(peptide);

buckets = zeros(1,maxPeptideInt+1);

for count = 1:seqLength
    buckets(peptideint(count)) = buckets(peptideint(count)) + 1;
end

for count = 1:20
    aminoAcids.(upper(int2aa(count))) = buckets(count);
end

% For unknown symbols (bucket 26) a warning is always shown:
if buckets(end)
    if ischar(peptide)
        unkn =  unique(peptide(peptideint==maxPeptideInt+1));
        warning(message('bioinfo:aacount:UnknownSymbols', unkn));
    else
        warning(message('bioinfo:aacount:UnknownNumericSymbols'));
    end
end

% For ambiguous symbols (bucket 21 to 23):
switch countAmbiguous % {'Ignore','Bundle','Individual','Warn'}
    case 1
        % do nothing
    case 2
        aminoAcids.Ambiguous = sum(buckets(21:23));
    case 3
        aminoAcids.B = buckets(aa2int('b'));
        aminoAcids.Z = buckets(aa2int('z'));
        aminoAcids.X = buckets(aa2int('x'));
    case 4
        buckets_with_amb = find(buckets(21:23))+20;
        if ~isempty(buckets_with_amb)
            warning(message('bioinfo:aacount:AmbiguousSymbols', int2aa( buckets_with_amb )));
        end
    case 5
        wx = buckets(aa2int('x'))/20;
        wb = buckets(aa2int('b'))/2;
        wz = buckets(aa2int('z'))/2;
        for count = 1:20
            aa = upper(int2aa(count));
            aminoAcids.(aa) = aminoAcids.(aa)+wx;
        end
        aminoAcids.N = aminoAcids.N + wb;
        aminoAcids.D = aminoAcids.D + wb;
        aminoAcids.Q = aminoAcids.Q + wz;
        aminoAcids.E = aminoAcids.E + wz;
end

% Count gaps when required
if countGaps
   aminoAcids.Gaps = buckets(aa2int('-'));
end

if pieChart
    pielabels = fieldnames(aminoAcids);
    piecounts = structfun(@(x) x, aminoAcids);
    h = piecounts>0;
    pie(piecounts(h),pielabels(h))
elseif barChart
    barlabels = fieldnames(aminoAcids);
    barcounts =  structfun(@(x) x, aminoAcids);
    h = barcounts>=0;
    hg = bar(barcounts(h));
    set(get(hg,'parent'),'Xtick',1:sum(h),'XTickLabel',barlabels(h),'Xlim',[0 sum(h)+1]);
end
