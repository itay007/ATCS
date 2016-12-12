function out = localalign(seq1, seq2, varargin)
%LOCALALIGN returns local alignments between two sequences. 
%
% LOCALALIGN(SEQ1, SEQ2) finds the optimal local alignment between two
% sequences, SEQ1 and SEQ2. By default, LOCALALIGN returns the
% highest-scoring local alignment and related information. To retrieve
% multiple local alignments, limit the number of alignments by using the
% option NUMALN, MINSCORE or PERCENT (see below). When using these options,
% the alignments are returned in decreasing order of scores. The scale
% factor used to calculate the score is provided by the scoring matrix info
% (see below). If the scale factor is not defined either by the scoring
% matrix or the scale option, then LOCALALIGN returns the raw score. For
% optimal performances, use the shorter sequence as the first input SEQ1.
%
% OUT = LOCALALIGN(SEQ1, SEQ2) returns a MATLAB structure with the
% information relative to the local alignment(s) between two sequences SEQ1
% and SEQ2. OUT is a structure array where each element has the following
% fields:
%     Score     - score of the alignment.
%     Start     - starting point indices indicating the starting points of the
%                 alignment in the two sequences.
%     Stop      - ending point indices indicating the ending points of the 
%                 alignment in the two sequences.
%     Alignment - 3xN char array representing the alignment between the two
%                 sequences.
%
% LOCALALIGN(...,'DOALIGNMENT', TF) lets you specify whether to return the
% actual alignment(s) or not in the field 'Alignment' of the output
% structure. Default is true.
%
% LOCALALIGN(..., 'NUMALN', M) specifies the number of local alignments to
% return. By default only one alignment with the best score is returned.
% This option and the MINSCORE and PERCENT options are mutually exclusive.
% 
% LOCALALIGN(..., 'MINSCORE', S) specifies the minimum score of the local
% alignments to return. This option and the NUMALN and PERCENT options are
% mutually exclusive. 
% 
% LOCALALIGN(..., 'PERCENT', P) returns local alignments whose score is
% within P percent from the highest score. This option and the NUMALN and
% MINSCORE options are mutually exclusive.
% 
% LOCALALIGN(...,'ALPHABET', alpha) specifies whether the sequences are
% amino acids ('AA') or nucleotides ('NT'). Default is 'AA'.
% 
% LOCALALIGN(...,'SCORINGMATRIX', matrix) specifies the scoring matrix to
% be used for  the alignment(s). The default is BLOSUM50 for AA or NUC44 for
% NT.
% 
% LOCALALIGN(...,'SCALE', scale) indicates the scale factor of the scoring
% matrix to return the score using arbitrary units. If the scoring matrix
% info also provides a scale factor, then both are used.
% 
% LOCALALIGN(...,'GAPOPEN', penalty) specifies the penalty for opening a
% gap in the alignment. The default gap open penalty is 8.
% 
% 
%   Examples:
%       
%       % Retrieve the best 3 local alignments.
%       out = localalign('VSPAGMASGYDPGKA','IPGKATREYDVSPAG', 'numaln', 3)
%       
%       % Retrieve information relative to local alignments with score higher than 8.
%       out = localalign('VSPAGMASGYDPGKA','IPGKATREYDVSPAG', 'minscore', 8, 'doalignment', false)
%
%       % Retrieve all local alignments with score within 15% from the maximum score. 
%       out = localalign('VSPAGMASGYDPGKA','IPGKATREYDVSPAG', 'percent', 15)
%
%       % Retrieve the best local alignment using BLOSUM30 as scoring matrix.
%       out = localalign('VSPAGMASGYDPGKA','IPGKATREYDVSPAG','scoringmatrix', 'blosum30')
%
%       % Retrieve the best 3 local alignments and specify scoring matrix
%       % and gap opening penalty.
%       seq1 = 'CCAATCTACTACTGCTTGCAGTAC';
%       seq2 = 'AGTCCGAGGGCTACTCTACTGAAC';
%       sm = [10 -9 -9 -9; -9 10 -9 -9; -9 -9 10 -9; -9 -9 -9 10];
%       out = localalign(seq1, seq2, 'alpha', 'nt', 'gapopen', 20, ...
%             'scoringmatrix', sm, 'numaln', 3)
%             
%      
%   See also ALIGNDEMO, BLOSUM, MULTIALIGN, NWALIGN, PAM, SEQDOTPLOT,
%   SHOWALIGNMENT, SWALIGN.

%  References:
%  G. Barton. An efficient algorithm to locate all locally optimal
%  alignments between two sequences allowing for gaps. CABIOS (1993)
%  9:729-734.
%
%   Copyright 2009-2012 The MathWorks, Inc.

%=== check inputs
bioinfochecknargin(nargin, 2, mfilename)
[scoringMatrix, gapopen, isAminoAcid, scale, doAlignment, ...
	setNumAln, setMinScore, setPct, limitAln] = parse_inputs(varargin{:});

%=== extract the Sequence data from structure
if isstruct(seq1)
	seq1 = bioinfoprivate.seqfromstruct(seq1);
end
if isstruct(seq2)
	seq2 = bioinfoprivate.seqfromstruct(seq2);
end

%=== handle properly "?" characters typically found in pdb files
if isAminoAcid
	if ischar(seq1)
		seq1 = strrep(seq1,'?','X');
	else
		seq1(seq1 == 26) = 23;
	end
	if ischar(seq2)
		seq2 = strrep(seq2,'?','X');
	else
		seq2(seq2 == 26) = 23;
	end
end

%=== check input sequences for consistency in their alphabet
if isAminoAcid && ~(bioinfoprivate.isaa(seq1) && bioinfoprivate.isaa(seq2))
	error(message('bioinfo:localalign:InvalidAminoAcidSequences'));
elseif ~isAminoAcid && ~(bioinfoprivate.isnt(seq1) && bioinfoprivate.isnt(seq2))
	error(message('bioinfo:localalign:InvalidNucleotideSequences'));
end

%=== use numerical arrays for easy indexing
if ischar(seq1)
	seq1 = upper(seq1); %the output alignment will be all uppercase
	if isAminoAcid
		intseq1 = aa2int(seq1);
	else
		intseq1 = nt2int(seq1);
	end
else
	intseq1 = uint8(seq1);
	if isAminoAcid
		seq1 = int2aa(intseq1);
	else
		seq1 = int2nt(intseq1);
	end
end
if ischar(seq2)
	seq2 = upper(seq2); %the output alignment will be all uppercase
	if isAminoAcid
		intseq2 = aa2int(seq2);
	else
		intseq2 = nt2int(seq2);
	end
else
	intseq2 = uint8(seq2);
	if isAminoAcid
		seq2 = int2aa(intseq2);
	else
		seq2 = int2nt(intseq2);
	end
end

m = length(seq1);
n = length(seq2);
if ~n||~m
	error(message('bioinfo:localalign:InvalidLengthSequences'));
end

% if m > n  
% 	warning('bioinfo:localalign:SwapInputs', ...
% 		'For optimal performances, swap the two input sequences so that the shorter sequence is the first input.');
% end

%=== set the scoring matrix
if isempty(scoringMatrix)
	if isAminoAcid
		[scoringMatrix, scoringMatrixInfo] = blosum50;
	else
		[scoringMatrix, scoringMatrixInfo] = nuc44;
	end
else
	if ~isnumeric(scoringMatrix)
		try
			[scoringMatrix, scoringMatrixInfo] = feval(scoringMatrix);
		catch allExceptions
			error(message('bioinfo:localalign:InvalidScoringMatrix'));
		end
		
	end
end

%=== get the scale from scoringMatrixInfo, if it exists
if exist('scoringMatrixInfo','var') && isfield(scoringMatrixInfo,'Scale')
	scale = scale * scoringMatrixInfo.Scale;
end

%=== scale the minscore, if needed
if setMinScore
	limitAln = limitAln / scale;
end

%=== make sure scoringMatrix can handle unknown, ambiguous or gaps

% possible values are
% B  Z  X  * -
% 21 22 23 24 25

scoringMatrixSize = size(scoringMatrix,1);
highestVal = max([intseq1, intseq2]);
if highestVal > scoringMatrixSize
	
	%=== map to the symbol corresponding to 'Any', if the matrix allows it
	if isAminoAcid
		anyVal = aa2int('X');
	else
		anyVal = nt2int('N');
	end
	if scoringMatrixSize >= anyVal
		intseq1(intseq1>scoringMatrixSize) = anyVal;
		intseq2(intseq2>scoringMatrixSize) = anyVal;
	else
		error(message('bioinfo:localalign:InvalidSymbolsInInputSequences'));
	end
end

%=== call the core functions
if doAlignment
	[score,  allPath(:,1), allPath(:,2)] = localbartonsimplegapmex(intseq1,intseq2, ...
		gapopen, scoringMatrix, doAlignment, setNumAln, setMinScore, setPct, limitAln);
	
	%=== warn if no alignment was found
	if isempty(allPath) || isempty(allPath{1,1})
        % mex function returns {} when there are not alignment above
        % setMinScore and {[],[]} whne there are not positive alignment, we
        % catch both cases and give a general warning (G561420)
		allScore = zeros(0,1);
		allStart = zeros(0,2); 
		allStop = zeros(0,2);
		allAln = cell(0,1);
		
		warning(message('bioinfo:localalign:EmptyAlignment'))
	else
		%=== remove empty entries (occurs when numAln is set > than actual number of possible aln)
		isNotEmptyEntry = sum(cellfun(@isempty, allPath), 2) == 0;
		allPath = allPath(isNotEmptyEntry,:);
		
		%=== initialize output variables
		numAln = size(allPath,1); % recompute numAln after removing empty entries
		allScore = zeros(numAln,1);
		allStart = zeros(numAln,2);
		allStop = zeros(numAln,2);
		allAln = cell(numAln,1);
		
		%=== reconstruct alignments from the path information
		for a = 1:numAln
			
			alnLen = size(allPath{a},2);
			path = zeros(alnLen,2);
			path(:,1) = allPath{a,1}';
			path(:,2) = allPath{a,2}';
			path = path(sum(path,2)>0,:); % ignore zero entries padding the arrays
			path = flipud(path); % reverse order, starting from smaller sequence positions
			
			%=== set the size of the alignment
			aln = repmat(('- -')',1,size(path,1));
			
			%=== add sequences to alignment
			aln(1,path(:,1)>0) = seq1(path(path(:,1)>0,1));
			aln(3,path(:,2)>0) = seq2(path(path(:,2)>0,2));
			
			%=== find positions where there are no gaps
			h = find(all(path>0,2));
			
			if isAminoAcid
				noGaps1 = aa2int(aln(1,h));
				noGaps2 = aa2int(aln(3,h));
			else
				noGaps1 = nt2int(aln(1,h));
				noGaps2 = nt2int(aln(3,h));
			end
			
			%=== remove symbols that cannot be scored
			htodel = max([noGaps1;noGaps2]) > scoringMatrixSize;
			h(htodel) = [];
			noGaps1(htodel) = [];
			noGaps2(htodel) = [];
			
			%=== score pairs with no gap
			value = scoringMatrix(sub2ind(size(scoringMatrix),double(noGaps1),double(noGaps2)));
			
			%=== insert symbols of the match string into the alignment
			aln(2,h(value >= 0)) = ':';
			aln(2,h(noGaps1 == noGaps2)) = '|';
			
			%=== assign output
			allScore(a) = score(a) * scale;
			allStart(a,1:2) = path(1,:)';
			allStop(a,1:2) = path(end,:)';
			allAln{a} = aln;
			
		end
	end
		
else % do not compute alignments, only starts and stops
	[allScore, allStart, allStop] = localbartonsimplegapmex(intseq1,intseq2, ...
		gapopen, scoringMatrix, doAlignment, setNumAln, setMinScore, setPct, limitAln);
	 allStart = reshape(allStart, size(allScore, 2), 2)+1; 
	 allStop = reshape(allStop, size(allScore, 2), 2);
	 allScore = allScore' * scale;
	
end

out.Score = allScore;
out.Start = allStart;
out.Stop = allStop;

if doAlignment
	out.Alignment = allAln;
end



%--------------------------------------------------------------------------

function [scoringMatrix, gapopen, isAminoAcid, scale, doAlignment, ...
	setNumAln, setMinScore, setPct, limitAln] = ...
	parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
gapopen = -8;
isAminoAcid = true;
scale = 1;
scoringMatrix = [];
doAlignment = true;
setNumAln = false;
setMinScore = false;
setPct = false;

%=== check for the right number of inputs
if rem(nargin,2) == 1
	error(message('bioinfo:localalign:IncorrectNumberOfArguments', mfilename));
end

%=== allowed parameters
okargs = {'scoringmatrix','gapopen','alphabet','scale', 'doalignment', ...
	'numaln', 'minscore', 'percent'};

%=== parse inputs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	k = find(strncmpi(pname, okargs,numel(pname)));
	if isempty(k)
		error(message('bioinfo:localalign:UnknownParameterName', pname));
	elseif length(k)>1
		error(message('bioinfo:localalign:AmbiguousParameterName', pname));
	else
		switch(k)
			case 1  % scoring matrix
				if isnumeric(pval)
					scoringMatrix = pval;
				else
					scoringMatrix = lower(pval);
				end
			case 2 % gap open penalty
				gapopen = -pval;
			case 3 % alphabet
    			isAminoAcid = bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename);
			case 4 % scale
				scale = pval;
			case 5 % doalignment
				doAlignment = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
			case 6 % numaln
				if isnumeric(pval) && isscalar(pval) && pval > 0 && rem(pval,1) == 0 && pval <= 2^12
					limitAln = pval;
					setNumAln = true;
				else
					error(message('bioinfo:localalign:InvalidNumAln'));
				end
			case 7 % minscore
				if isnumeric(pval) && isscalar(pval) && pval >= 0 
					limitAln = pval;
					setMinScore = true;
				else
					error(message('bioinfo:localalign:InvalidMinScore'));
				end
			case 8 % percent
				if isnumeric(pval) && isscalar(pval) && pval >=0 && pval <=100
					limitAln = pval / 100;  
					setPct = true;
				else
					error(message('bioinfo:localalign:InvalidPercent'));
				end
		end
	end
end

%=== make sure that numaln, minscore and percent are mutually exclusive 
limitCheck = [setNumAln setMinScore setPct] == true;
if sum(limitCheck) == 0 % no criterion has been set
	setNumAln = true; % use NUMALN
	limitAln = 1;
elseif sum(limitCheck) > 1 % more than one criterion has been set
	error(message('bioinfo:localalign:InvalidParCombination'));
end





