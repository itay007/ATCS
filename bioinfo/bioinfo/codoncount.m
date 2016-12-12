function [outval,carray] = codoncount(dna,varargin)
%CODONCOUNT report codon counts for a sequence.
%
%   CODONCOUNT(SEQ) counts the number of occurrences of each codon in the
%   sequence and displays a formatted table of the result. Codons with
%   ambiguous nucleotide symbols are not counted by default. Gaps (-) are
%   removed from the input sequence. Codons with other unrecognized
%   characteres are not counted, but a warning message is displayed.
%
%   CODONS = CODONCOUNT(SEQ) returns these codon counts in a structure with
%   the fields AAA, AAC, AAG, ..., TTG, TTT. 
%
%   [CODONS, CARRAY] = CODONCOUNT(SEQ) returns a 4x4x4 array of the raw
%   count data for each codon. The three dimensions correspond to the three
%   positions in the codon. The index in each dimension corresponds to
%   nucleotides in the 'ACGT' order. For example the (2,3,4) element of the
%   array gives the number of 'CGT' codons in the input sequence.
%
%   CODONCOUNT(...,'FRAME',F) returns the codon count for reading frame
%   F, where F is 1, 2, or 3. Default is 1.
%
%   CODONCOUNT(...,'REVERSE',true) returns the codon count for the
%   reverse complement of SEQ.
%
%   CODONCOUNT(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
%   nucleotide symbols are present in a codon. Options are: 'Ignore' skips
%   codons with ambiguous symbols, 'Bundle' counts and bundles them into
%   the Ambiguous field of the output structure, 'Prorate' counts and
%   prorates them into the other codons with standard nucleotide symbols,
%   and 'Warn' ignores them and display a warning message. Default is
%   'Ignore'.
%
%   CODONCOUNT(...,'FIGURE',true) creates a figure showing a heat map
%   of the codon counts.
%
%   CODONCOUNT(...,'GENETICCODE',CODE) overlays a grid on the figure
%   grouping the synonymous codons according with the genetic code CODE.
%   Default is 'Standard' or 1. Set CODE to 'None' to create a heat map
%   without showing the grid.
%
%   Examples:
%
%       codons = codoncount('AAACGTTA')
%
%       r2codons = codoncount('AAACGTTA','Frame',2,'Reverse',true)
%
%   See also AACOUNT, BASECOUNT, BASELOOKUP, CODONBIAS, DIMERCOUNT,
%   NMERCOUNT, NTDENSITY, SEQRCOMPLEMENT, SEQSHOWORFS, SEQSTATSDEMO,
%   SEQWORDCOUNT.

%   Copyright 2002-2012 The MathWorks, Inc.

% check inputs

frame = 1;
reverseComp = false;
showFig = false;
countAmbiguous = 1; % {'Ignore','Bundle','Warn','Prorate'}
geneticCode = 1;

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
        error(message('bioinfo:codoncount:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'frame','reverse','figure','ambiguous','geneticcode'};
    % Loop over the values
    for j=1:2:nargin-2
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case 1  % frame
                if ~isnumeric(pval) || pval > 3 || pval < 1
                    error(message('bioinfo:codoncount:BadFrameNumber'))
                end
            frame = pval;
            case 2  % direction
                reverseComp = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 3  % figure
                showFig = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 4 % ambiguous
                countAmbiguous = bioinfoprivate.optPartialMatch(pval,{'Ignore','Bundle','Warn','Prorate'}, okargs{k}, mfilename);
            case 5 % geneticCode
                geneticCode = pval;
                if ischar(geneticCode)
                    if strncmpi(geneticCode,'none',numel(geneticCode))
                        geneticCode = 0;
                    end
                end
        end
    end
end

if reverseComp
    dna = seqrcomplement(dna);
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

% Use 3D indexing to count the codons
codons = zeros(maxNtInt+1,maxNtInt+1,maxNtInt+1);
for count = frame:3:seqLength-2
    codons(dnaint(count),dnaint(count+1),dnaint(count+2)) = codons(dnaint(count),dnaint(count+1),dnaint(count+2)) + 1;
end

% Place info in the output structure
labelcount = 0;
labels = cell(1,64);
for first = 1:4
    for second = 1:4
        for third = 1:4
            codon = int2nt([first second third]);
            labelcount = labelcount + 1;
            labels{labelcount} = codon;
            output.(codon) = codons(first,second,third);
        end
    end
end

% For unknown symbols (codons [:,:,17] [:,17,:] and [17,:,:]) a warning is always shown:
if any(any(codons(:,:,end))) || any(any(codons(:,end,:))) || any(any(codons(:,:,end)))
    if ischar(dna)
        unkn =  unique(dna(dnaint==maxNtInt+1));
        warning(message('bioinfo:codoncount:UnknownSymbols', unkn));
    else
        warning(message('bioinfo:codoncount:UnknownNumericSymbols'));
    end
end

% For ambiguous symbols:
switch countAmbiguous % {'Ignore','Bundle','Warn','Prorate'}
    case 1
        % do nothing
    case 2
        output.Ambiguous = sum(sum(sum(codons(1:15,1:15,1:15))))-sum(sum(sum(codons(1:4,1:4,1:4))));
    case 3
        codons_with_amb  = find(any(any(codons(5:15,1:15,1:15),2),3) | ...
                               any(any(codons(1:15,5:15,1:15),1),3)' | ...
                               squeeze(any(any(codons(1:15,1:15,5:15),1),2)))+4;
        if ~isempty(codons_with_amb)
            warning(message('bioinfo:codoncount:AmbiguousSymbols', int2nt( codons_with_amb(:)' )));
        end
    case 4
        w = 1./[1;2;2;2;3;3;3;4];
        codons(1:15,1:15,nt2int('a')) = sum(bsxfun(@times,codons(1:15,1:15,nt2int('armwdhvn')),reshape(w,1,1,8)),3);
        codons(1:15,1:15,nt2int('c')) = sum(bsxfun(@times,codons(1:15,1:15,nt2int('cymsbhvn')),reshape(w,1,1,8)),3);
        codons(1:15,1:15,nt2int('g')) = sum(bsxfun(@times,codons(1:15,1:15,nt2int('grksbdvn')),reshape(w,1,1,8)),3);
        codons(1:15,1:15,nt2int('t')) = sum(bsxfun(@times,codons(1:15,1:15,nt2int('tykwbdhn')),reshape(w,1,1,8)),3);
        codons = codons(1:15,1:15,1:4);
        codons(:,nt2int('a'),:) = sum(bsxfun(@times,codons(:,nt2int('armwdhvn'),:),w'),2);
        codons(:,nt2int('c'),:) = sum(bsxfun(@times,codons(:,nt2int('cymsbhvn'),:),w'),2);
        codons(:,nt2int('g'),:) = sum(bsxfun(@times,codons(:,nt2int('grksbdvn'),:),w'),2);
        codons(:,nt2int('t'),:) = sum(bsxfun(@times,codons(:,nt2int('tykwbdhn'),:),w'),2);
        codons = codons(:,1:4,:);
        codons(nt2int('a'),:,:) = sum(bsxfun(@times,codons(nt2int('armwdhvn'),:,:),w),1);
        codons(nt2int('c'),:,:) = sum(bsxfun(@times,codons(nt2int('cymsbhvn'),:,:),w),1);
        codons(nt2int('g'),:,:) = sum(bsxfun(@times,codons(nt2int('grksbdvn'),:,:),w),1);
        codons(nt2int('t'),:,:) = sum(bsxfun(@times,codons(nt2int('tykwbdhn'),:,:),w),1);
        codons = codons(1:4,:,:);
        
        % Place info in the output structure
        for first = 1:4
            for second = 1:4
                for third = 1:4
                    codon = int2nt([first second third]);
                    output.(codon) = codons(first,second,third);
                end
            end
        end
end

% Create pretty output when no outputs.
if nargout == 0
    numSpaces = ceil(log10(max(codons(:))));
    if any(rem(codons(:),1))
        formatString = sprintf('%%s - %%%d.2f     ',numSpaces);
    else
        formatString = sprintf('%%s - %%%dd     ',numSpaces);
    end
    outputStr = '';
    for outer = 0:3
        for middle = 0:3
            for inner = 1:4
                step = 16 * outer + 4 * middle + inner;
                outputStr = [outputStr, sprintf(formatString,labels{step},codons(outer+1,middle+1,inner))]; %#ok<AGROW>
            end
            outputStr = [outputStr,sprintf('\n')]; %#ok<AGROW>
        end
    end

    if isfield(output,'Ambiguous')
        outputStr = [outputStr, sprintf('Ambiguous - %d\n',output.Ambiguous)];
    end
    disp(outputStr);
else
    outval = output;
end

if nargout > 1 || showFig
    carray = codons(1:4,1:4,1:4);
end

if showFig
    ord2d = [ 1  2  5  6 17 18 21 22
              3  4  7  8 19 20 23 24
              9 10 13 14 25 26 29 30
             11 12 15 16 27 28 31 32
             33 34 37 38 49 50 53 54
             35 36 39 40 51 52 55 56
             41 42 45 46 57 58 61 62
             43 44 47 48 59 60 63 64];
    carray2(:,:,1) = squeeze(carray(1,:,:))';
    carray2(:,:,2) = squeeze(carray(2,:,:))';
    carray2(:,:,3) = squeeze(carray(3,:,:))';
    carray2(:,:,4) = squeeze(carray(4,:,:))';
    im = carray2(ord2d);
    labels = labels(ord2d);
    if geneticCode~=0
       aal = reshape(char(nt2aa(labels,'geneticcode',geneticCode,'alt',false)),8,8);
       gc = geneticcode(geneticCode);
       as = gc.Starts;
    end
    imagesc(im);
    axis off;
    colormap(bone);
    colorbar;
    for i = 1:8
        for j = 1:8
            if any(geneticCode~=0) && aal(i,j) == '*' 
               text(j,i,labels{i,j},'color',[.8 0 0],...
                 'horizontalAlignment','center');
            elseif any(geneticCode~=0) && any(strcmp(labels{i,j},as))
               text(j,i,labels{i,j},'color',[.1 .6 .1],...
                 'horizontalAlignment','center');
            else
               text(j,i,labels{i,j},'color',[.1 .1 1],...
                  'horizontalAlignment','center');
            end
        end
    end
    if any(geneticCode~=0)
      xvl = reshape(ones(16,1)*(1:7)+.5,2,56);
      yvl = repmat([(1:8)-.5;(1:8)+.5],1,7);
      xhl = yvl;
      yhl = xvl;
      yvl(:,diff(double(aal),[],2)==0) = NaN;
      yhl(:,diff(double(aal),[],1)'==0) = NaN;
      xhl = [reshape(xhl,16,7);nan(1,7)];
      yhl = [reshape(yhl,16,7);nan(1,7)];
      xvl = [reshape(xvl,16,7);nan(1,7)];
      yvl = [reshape(yvl,16,7);nan(1,7)];
      hold on
      plot([xhl(:);xvl(:);.5;.5;8.5;8.5;.5],[yhl(:);yvl(:);.5;8.5;8.5;.5;.5],...
           'color',[.1 .1 1],'Linewidth',2);
      hold off 
      text(8.5,9,['Genetic Code: ' gc.Name],'color',[.1 .1 1],...
          'horizontalAlignment','right','fontSize',8);
    end
end

