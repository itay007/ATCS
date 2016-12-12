function [dimers,matrix] = dimercount(dna,varargin)
%DIMERCOUNT report nucleotide dimer counts for a sequence.
%
%   DIMERCOUNT(SEQ) counts the number of occurrences of nucleotide dimers
%   in the sequence and returns these numbers in a structure with the
%   fields AA, AC, AG,..., GT, TT. Dimers with ambiguous nucleotide symbols
%   are not counted by default. Gaps (-) are removed from the input
%   sequence. Dimers with other unrecognized characteres are not counted,
%   but a warning message is displayed.
%
%   [DIMERS, P] = DIMERCOUNT(SEQ) returns a 4x4 matrix of the relative
%   proportions of the dimers in SEQ, with the rows corresponding to A,C,G
%   and T in as the first element of the dimer and the columns
%   corresponding to A,C,G, and T in the second element.
%
%   DIMERCOUNT(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
%   nucleotide symbols are present in a dimer. Options are: 'Ignore' skips
%   dimers with ambiguous symbols, 'Bundle' counts and bundles them into
%   the Ambiguous field of the output structure, 'Prorate' counts and
%   prorates them into the other dimers with standard nucleotide symbols,
%   and 'Warn' ignores them and display a warning message. Default is
%   'Ignore'.
%
%   DIMERCOUNT(...,'CHART',STYLE) creates a chart showing the relative
%   proportions of the dimers. Valid styles are 'Pie', 'Bar' and 'HeatMap'.
%
%   Example:
%
%       dimercount('TAGCTGGCCAAGCGAGCTTG')
%
%   See also AACOUNT, BASECOUNT, BASELOOKUP, CODONCOUNT, NMERCOUNT,
%   NTDENSITY, SEQSTATSDEMO.

%   Copyright 2002-2012 The MathWorks, Inc.

pieChart = false;
barChart = false;
heatMapChart = false;
countAmbiguous = 1; % {'Ignore','Bundle','Warn','Prorate'}

% If the input is a structure then extract the Sequence data.
if isstruct(dna)
    dna = bioinfoprivate.seqfromstruct(dna);
end

% Remove any gap
if ischar(dna)
    if any(dna=='-')
        dna = strrep(dna,'-','');
    end
else
    if any(dna==nt2int('-'))
        dna(dna==nt2int('-')) =[];
    end
end

if  nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:dimercount:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'chart','ambiguous'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case 1  % graph
                k = bioinfoprivate.optPartialMatch(pval,{'Pie','Bar','HeatMap'}, okargs{k}, mfilename);
                if k==1
                    pieChart = true;
                elseif k==2
                    barChart = true;
                else
                    heatMapChart = true;
                end
            case 2 % ambiguous
                countAmbiguous = bioinfoprivate.optPartialMatch(pval,{'Ignore','Bundle','Warn','Prorate'}, okargs{k}, mfilename);
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

buckets = zeros(maxNtInt+1,maxNtInt+1);
for count = 1:seqLength-1
    buckets(dnaint(count),dnaint(count+1)) = buckets(dnaint(count),dnaint(count+1)) + 1;
end

aint = nt2int('a');cint = nt2int('c');
gint = nt2int('g');tint = nt2int('t');
dimers.AA = buckets(aint,aint);
dimers.AC = buckets(aint,cint);
dimers.AG = buckets(aint,gint);
dimers.AT = buckets(aint,tint);
dimers.CA = buckets(cint,aint);
dimers.CC = buckets(cint,cint);
dimers.CG = buckets(cint,gint);
dimers.CT = buckets(cint,tint);
dimers.GA = buckets(gint,aint);
dimers.GC = buckets(gint,cint);
dimers.GG = buckets(gint,gint);
dimers.GT = buckets(gint,tint);
dimers.TA = buckets(tint,aint);
dimers.TC = buckets(tint,cint);
dimers.TG = buckets(tint,gint);
dimers.TT = buckets(tint,tint);

% For unknown symbols (buckets [:,17] and [17,:]) a warning is always shown:
if any(buckets(:,end)) || any(buckets(end,:))
    if ischar(dna)
        unkn =  unique(dna(dnaint==maxNtInt+1));
        warning(message('bioinfo:dimercount:UnknownSymbols', unkn));
    else
        warning(message('bioinfo:dimercount:UnknownNumericSymbols'));
    end
end

% For ambiguous symbols (buckets [1:4,5:15] , [5:15,1:4] and [5:15,5:15]:
switch countAmbiguous % {'Ignore','Bundle','Warn','Prorate'}
    case 1
        % do nothing
    case 2
        dimers.Ambiguous = sum(sum(buckets(1:15,1:15)))-sum(sum(buckets(1:4,1:4)));
    case 3
        buckets_with_amb  = find(any(buckets(5:15,1:15),2) | any(buckets(1:15,5:15),1)')+4;
        if ~isempty(buckets_with_amb)
            warning(message('bioinfo:dimercount:AmbiguousSymbols', int2nt( buckets_with_amb(:)' )));
        end
    case 4
        w = 1./[1;2;2;2;3;3;3;4];
        buckets = [buckets(1:15,nt2int('armwdhvn'))*w ...
                   buckets(1:15,nt2int('cymsbhvn'))*w ...
                   buckets(1:15,nt2int('grksbdvn'))*w ...
                   buckets(1:15,nt2int('tykwbdhn'))*w]';
        buckets = [buckets(:,nt2int('armwdhvn'))*w ...
                   buckets(:,nt2int('cymsbhvn'))*w ...
                   buckets(:,nt2int('grksbdvn'))*w ...
                   buckets(:,nt2int('tykwbdhn'))*w]';
        dimers.AA = buckets(aint,aint);
        dimers.AC = buckets(aint,cint);
        dimers.AG = buckets(aint,gint);
        dimers.AT = buckets(aint,tint);
        dimers.CA = buckets(cint,aint);
        dimers.CC = buckets(cint,cint);
        dimers.CG = buckets(cint,gint);
        dimers.CT = buckets(cint,tint);
        dimers.GA = buckets(gint,aint);
        dimers.GC = buckets(gint,cint);
        dimers.GG = buckets(gint,gint);
        dimers.GT = buckets(gint,tint);
        dimers.TA = buckets(tint,aint);
        dimers.TC = buckets(tint,cint);
        dimers.TG = buckets(tint,gint);
        dimers.TT = buckets(tint,tint);
end

if pieChart
    pielabels = fieldnames(dimers);
    piecounts = structfun(@(x) x, dimers);
    h = piecounts>0;
    pie(piecounts(h),pielabels(h))
elseif barChart
    hg = bar3(buckets(1:4,1:4));
    set(get(hg(1),'parent'),'XTickLabel',{'A','C','G','T'},'YTickLabel',{'A','C','G','T'});
    ylabel('First Base');
    xlabel('Second Base');
elseif heatMapChart    
    imagesc(buckets(1:4,1:4));
    axis off;
    colormap(bone);
    colorbar;
    t = ['A','C','G','T'];
    for i = 1:4
        for j = 1:4
            text(i,j,[t(j),t(i)],...
                'color','r','horizontalAlignment','center');
        end
    end
end

if nargout>1
    matrix = buckets(1:4,1:4);
    matrix = matrix ./ sum(matrix(:));
end
    