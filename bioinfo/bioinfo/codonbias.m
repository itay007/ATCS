function codon_bias = codonbias(sequence, varargin)
%CODONBIAS reports codon usage per amino acid for a DNA sequence.
%
%   CODONBIAS(SEQ) translates a DNA sequence in frame and returns the codon
%   frequencies for each amino acid as a structure.
%
%   CODONBIAS(...,'GENETICCODE',CODE) uses an alternative translation table
%   for the sequence. The default code is 'Standard' or 1.
%
%   CODONBIAS(...,'FRAME',F) returns the codon bias for reading frame
%   F, where F is 1, 2, or 3. Default is 1.
%
%   CODONBIAS(... ,'REVERSE',true) returns the codon usage bias for the
%   reverse complement of SEQ.
%
%   CODONBIAS(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
%   nucleotide symbols are present in a codon. Options are: 'Ignore' skips
%   codons with ambiguous symbols, 'Prorate' counts and prorates them into
%   the other codons with standard nucleotide symbols, and 'Warn' ignores
%   them and display a warning message. Default is 'Ignore'.
%
%   CODONBIAS(...,'PIE',true) creates a figure of 20 pie charts for each
%   amino acid filled with codon frequencies
%
%   Example:
%
%       S = getgenbank('m10051')
%       cb = codonbias(S.Sequence, 'PIE', true)
%
%   See also AMINOLOOKUP, CODONCOUNT, GENETICCODE, NT2AA.

%   Copyright 2003-2012 The MathWorks, Inc.


% determine the format of the sequence

if isstruct(sequence)
    seq = bioinfoprivate.seqfromstruct(sequence);
elseif(~bioinfoprivate.isnt(sequence))
    error(message('bioinfo:codonbias:IncorrectSequenceType'));
else
    seq = sequence;
end

%defaults
countAmbiguous = 1; % ignore
countAmbiguousStrs = {'Ignore','Bundle','Prorate'};
geneticCode = 1;
piech = false;
reverseComp= false;
frame = 1;

%%% check arguments
if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:codonbias:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'geneticcode','pie','frame','reverse','ambiguous'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case 1  % geneticcode
                if (length(pval)==1) && isnumeric(pval) && isreal(pval) && (pval>0)
                    geneticCode = pval;
                elseif ischar(pval) && isrow(pval)
                    geneticCode = pval;
                else
                    error(message('bioinfo:codonbias:InvalidGeneticCode'));
                end
            case 2  % pie
                piech = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 3  % frame
                if ~isnumeric(pval) || pval > 3 || pval < 1
                     error(message('bioinfo:codonbias:BadFrameNumber'))
                end
                frame = pval;
            case 4  % reverse
                reverseComp = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 5  % ambiguous
                countAmbiguous = bioinfoprivate.optPartialMatch(pval,{'Ignore','Warn','Prorate'}, okargs{k}, mfilename);
        end
    end
end



% count codons
codontally = codoncount(seq, 'FRAME', frame, 'REVERSE',reverseComp,'Ambiguous',countAmbiguousStrs{countAmbiguous});

% send warn when ambiguous symbols were found and it was asked for
if (countAmbiguous==2) && codontally.Ambiguous>0
    warning(message('bioinfo:codonbias:AmbiguousSeqCharacters'));
    codontally = rmfield(codontally, {'Ambiguous'});
end

%get freq and codons for each amino
all_codons = fieldnames(codontally);
intAA = double(aa2int(nt2aa(cell2mat(all_codons'), 'GENETICCODE', geneticCode)))';
codon_count = cell2mat(struct2cell(codontally));
total = accumarray(intAA,codon_count,[26 1]);
freq = codon_count ./ total(intAA);
aminos = cellstr(aminolookup(int2aa(1:25)'));

% put into struct
for i = 1:length(aminos)
    codon_bias.(aminos{i}).Codon = all_codons(intAA == i)';
    codon_bias.(aminos{i}).Freq = freq(intAA == i)';
end

% get rid of ambiguous fields
codon_bias = rmfield(codon_bias, {'Asx', 'Glx','Xaa', 'GAP'});

% make a pie chart
if(piech)
    aminos([25 23 22 21])= [];
    %make pie chart
    figure('Position', [170   183   765   649]);
    for  i = 1:length(aminos)-1 % no stop codons
        index =  find(intAA == i);
        subplot(4,5,i);
        X = freq(index);
        if(sum(X) > 0)
            z = (X ~= 0);
            Y = all_codons(index);
            pie(X(z),Y(z));
            title(aminos{i});
        else
            pie(1,{'none'})
            title(aminos{i});
        end
    end
    brighten(.7);
end


