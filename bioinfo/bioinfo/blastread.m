function data = blastread(blasttext)
%BLASTREAD reads NCBI BLAST format report files.
%
%   DATA = BLASTREAD(FILE) reads in a NCBI formatted BLAST report from
%   FILE and creates a structure, DATA, containing fields corresponding to
%   the BLAST keywords.
%
%   NOTE: FILE should be created with the NCBI 'Use new formatter' option.
%   FILE can also be a URL or a MATLAB character array that contains the
%   text output of a NCBI BLAST report.
%
%   DATA contains these fields:
%       RID
%       Algorithm
%       Query
%       Database
%       Hits.Name
%       Hits.Length
%       Hits.HSPs.Score
%       Hits.HSPs.Expect
%       Hits.HSPs.Identities
%       Hits.HSPs.Positives  (peptide sequences)
%       Hits.HSPs.Gaps
%       Hits.HSPs.Frame      (translated searches)
%       Hits.HSPs.Strand     (nucleotide sequences)
%       Hits.HSPs.Alignment  (3xn: Query- R1, Alignment- R2, Subject- R3)
%       Hits.HSPs.QueryIndices
%       Hits.HSPs.SubjectIndices
%       Statistics
%
%   Examples:
%
%       % Create a BLAST request with a GenPept accession number.
%       RID = blastncbi('AAA59174', 'blastp', 'expect', 1e-10)
%
%       % Then pass the RID to GETBLAST to download the report and save it
%       % to a text file.
%       getblast(RID,'TOFILE','AAA59174_BLAST.rpt')
%
%       % Using the saved file, read the results into a MATLAB structure.
%       results = blastread('AAA59174_BLAST.rpt')
%
%   See also BLASTFORMAT, BLASTLOCAL, BLASTNCBI, BLASTREADLOCAL, GETBLAST.

% Copyright 2003-2009 The MathWorks, Inc.



if ~ischar(blasttext) && ~iscellstr(blasttext)
	error(message('bioinfo:blastread:InvalidInput'));
end

if iscellstr(blasttext)
	blasttext = char(blasttext);
else %it is char, check if it has a url address or a file name
	if isvector(blasttext) && ~isempty(strfind(blasttext(1:min(10,end)), '://'))
		% must be a URL
		if (~usejava('jvm'))
			error(message('bioinfo:blastread:NoJava'))
		end
		try
			blasttext = urlread(blasttext);
		catch theErr
			error(message('bioinfo:blastread:CannotReadURL', blasttext));
		end
		% clean up any &amp s
		blasttext = strrep(blasttext,'&amp;','&');
	elseif isvector(blasttext) && exist(blasttext,'file')
		% the file exists
		fid = fopen(blasttext, 'rt');
		blasttext = fread(fid,'*char')';
		fclose(fid);
	end
end

% check whether error has occurred on the NCBI side (e.g. 500 Server Error)
if ~isempty(regexpi(blasttext, 'ERROR:'))
	error(message('bioinfo:blastread:Error'))
end

% If the input is a string of BLAST data then a score section must be present
if isempty(regexpi(blasttext,'Score =.*?','once'));
	if isempty(regexpi(blasttext,'No hits found','once'))
		error(message('bioinfo:blastread:NoValidData'))
	else
		warning(message('bioinfo:blastread:BlastSearchNull'));
		data = [];
		return;
	end
end

% Remove html tags
tags = find(blasttext == '<' | blasttext == '>');
h = find(diff(blasttext(tags)=='>')==1);
g = zeros(size(blasttext),'single');
g(tags(h)) = 1;
g(tags(h+1)+1) = g(tags(h+1)+1)-1;
blasttext = blasttext(~cumsum(g));

% Get the header info and cleanup
blastheader = blasttext(1:strfind(blasttext,'Sequences producing'));
blastheader = regexprep(blastheader,'\s+',' '); %remove extra white space

% Extract search information from blast header
data.RID = regexpi(blastheader, '(?<=RID\s*:\s*)[\w-\.]+','match','once'); % for backward compatibility with RID format used before April 16 2007
data.Algorithm = regexpi(blastheader,'(T*)BLAST[N|P|X]\s+\d+\.\d+\.\d+\s+[\[\(]\w{3}-\d{2}-\d{4}[\)\]]','match','once'); % allow both () and [] around date
if isempty(data.Algorithm) % as of May 2008, the algorithm revision format has changed
	data.Algorithm = regexpi(blastheader,'(T*)BLAST[N|P|X]\s+\d+\.\d+\.\d+\+?','match','once');
end
data.Query = regexpi(blastheader,'(?<=Query= ).*?(?=( Database|\s*Score| Distance| E-value))','match','once');
data.Database = regexpi(blastheader,'(?<=Database: ).*?(?= [\d,]*? sequences)','match','once');

% Get the hits section
align_start = strfind(lower(blasttext),'alignments');
blasttext = blasttext(align_start(size(align_start,2)):size(blasttext,2));

% Remove junk at the end of the report
blasttext = blasttext(1:regexpi(blasttext,'S2: (.*?) bits\)', 'end'));

% Save off the statistical data
stats_start = strfind(blasttext,'Database:');
if isempty(stats_start)
	stats_start = strfind(lower(blasttext),'lambda');
end
statstext = blasttext(stats_start(size(stats_start,1)):size(blasttext,2));
statstext = regexprep(statstext,'\r\s','\r'); %remove any extra line feeds
statsend = regexpi(statstext, 'S2: (.*?) bits\)', 'end');
statstext = statstext(1:statsend);

% Get alignment header info
% tags = sort([find(blasttext == '>'), strfind(blasttext,'Query')]);
% h = find(diff(blasttext(tags)=='Q')== 1);
% g = zeros(size(blasttext),'single');
% g(tags(h)) = 1;
% g(tags(h+1)) = g(tags(h+1))-1;
% alignheaders = [blasttext(logical(cumsum(g))) ' >'];

% Get alignment header info (including hits for which the alignment is skipped)
alignheaders = regexp([blasttext ' >'], '>(.*?)(?=Query|>)', 'match'); 
alignheaders = [alignheaders{:} ' >'];

% Remove and uniformize extra white space
clean = '';
chunk = 10000;
filesize = size(alignheaders,2);
for i = 1:chunk:filesize
	clean = [clean regexprep(alignheaders(i:min(i+chunk-1,filesize)),'\s+',' ')]; %remove extra white space
end

% Correct abbreviated e-scores
clean = regexprep(clean,' (e[-+])',' 1$1');

use_frame = false;
use_positives = false;
use_strand = false;
use_gaps = false;
skipWarnFlag = false; % keep track whether a warning has already been thrown for skipped aln


% Frame is included in translated searches
if strfind(clean,'Frame =')
	use_frame = true;
end

% Positives are included with peptide sequences
if strfind(clean,'Positives =')
	use_positives = true;
end

% Strand is included with nucleotide sequences
if ~isempty(strfind(clean,'Strand='))||~isempty(strfind(clean,'Strand ='))
	use_strand = true;
end

% Gaps
if ~isempty(strfind(clean,'Gaps ='))
	use_gaps = true;
end

%hit_lines = regexpi(clean,'>.*?(?=(>|database:|lambda\s+k\s+h))','match');
% The new regexp will also match an alignment header that may contain an
% '>' before the word 'score'
hit_lines = regexpi(clean,'>.*?(?=(score))score.*?(?=(>|database:|lambda\s+k\s+h))','match');
names_lengths = regexp(hit_lines,'>.*?(?= Score)','match')';
num_hits = size(hit_lines,2);

for i = 1:num_hits
	
	name = regexp(names_lengths{i},'(?<=>).*?(?= Length)','match','once');
	hit_length = regexp(names_lengths{i},'(?<=Length ?= ?)\d*','match','once');
	data.Hits(i).Name = char(name{1});
	data.Hits(i).Length = str2double(char(hit_length{1}));
	align_start = strfind(blasttext,strtok(data.Hits(i).Name));
	
	if i == num_hits
		align_stop = size(blasttext,2);
	else
		next_name = char(regexp(names_lengths{i+1},'(?<=>).*?(?= Length)','match','once'));
		align_stop = strfind(blasttext,strtok(next_name));
		
		% if two entries have same name (e.g. one has skipped alignment)
		if length(align_start) > 1
			if length(align_stop) > 1 % must be the 1st of two entries
				align_start = align_start(1);
				align_stop = align_stop(2);
			else % must be the 2nd of the two entries
				align_start = align_start(2); % stop_align obtained from different next_name
			end
			
		end
		
	end
	hsp_align = blasttext(align_start:align_stop);
	hsp_lines = regexp(hsp_align,'Score.*?(?=( Score|.$))','match');
	num_hsps = size(hsp_lines,2);
	
	for j = 1:num_hsps
		score = regexp(hsp_lines{j},'Score =\s{1,2}[e+\d\.]*','match');
		expect = regexp(hsp_lines{j},'Expect =\s{1,2}[e-\d\.]*','match');
		ident = regexp(hsp_lines{j},'Identities =\s\d*/\d* \(\d*%\)','match');
		
		if isempty(ident) % most likely the sequence no longer exists in the database
			data.Hits(i).HSPs(j).Score = sscanf(char(score), '%*s = %f');
			data.Hits(i).HSPs(j).Expect = sscanf(char(expect), '%*s = %f');
			data.Hits(i).HSPs(j).Identities = struct('Match', [], 'Possible', [], 'Percent', []);
			
			if ~skipWarnFlag
				skipWarnFlag = true;
				warning(message('bioinfo:blastread:SkippedAln'));
			end
			
			if use_positives
				data.Hits(i).HSPs(j).Positives = struct('Match', [], 'Possible', [], 'Percent', []);
			end
			
			if use_gaps
				data.Hits(i).HSPs(j).Gaps = struct('Match', [], 'Possible', [], 'Percent', []);
			end
			
			if use_frame
				data.Hits(i).HSPs(j).Frame = [];
			end
			
			if use_strand
				data.Hits(i).HSPs(j).Strand = [];
			end
			
			data.Hits(i).HSPs(j).Alignment = [];
			data.Hits(i).HSPs(j).QueryIndices = [];
			data.Hits(i).HSPs(j).SubjectIndices = [];
			
			
		else
			
			if use_positives
				% Get the 3 values for positives
				posit = regexp(hsp_lines{j},'Positives = \d*/\d* \(\d*%\)','match');
			end
			if use_frame
				frame = regexp(hsp_lines{j},'Frame = [+-]\d(\s*/\s*[+-]\d)?','match');
			end
			if use_strand
				% Get the strand orientation
				strand = regexp(hsp_lines{j},'Strand\s*=\s*(Plus|Minus)\s*/\s*(Plus|Minus)','match');
			end
			query_indices = regexp(hsp_lines{j},'Query  (\d+)\s+[a-zA-Z\-\*]+  (\d+)','tokens');
			subj_indices = regexp(hsp_lines{j},'Sbjct  (\d+)\s+[a-zA-Z\-\*]+  (\d+)','tokens');
			alignments = regexp(hsp_lines{j},'Query  \d+\s+([a-zA-Z\-\*]+)  \d+\s([a-zA-Z\+\*|\s]+)\sSbjct  \d+\s+([a-zA-Z\-\*]+)  \d+','tokens');
			
			data.Hits(i).HSPs(j).Score = sscanf(char(score), '%*s = %f');
			data.Hits(i).HSPs(j).Expect = sscanf(char(expect), '%*s = %f');
			ident_vals = sscanf(char(ident), '%*s = %d/%d (%d%)');
			data.Hits(i).HSPs(j).Identities.Match = ident_vals(1);
			data.Hits(i).HSPs(j).Identities.Possible = ident_vals(2);
			data.Hits(i).HSPs(j).Identities.Percent = ident_vals(3);
			
			if use_positives
				posit_vals = sscanf(char(posit), '%*s = %d/%d (%d%)');
				data.Hits(i).HSPs(j).Positives.Match = posit_vals(1);
				data.Hits(i).HSPs(j).Positives.Possible = posit_vals(2);
				data.Hits(i).HSPs(j).Positives.Percent = posit_vals(3);
			end
			
			if use_gaps
				gaps = regexp(hsp_lines{j},'Gaps = \d*/\d* \(\d*%\)','match');
				gap_vals = sscanf(char(gaps), '%*s = %d/%d (%d%)');
				data.Hits(i).HSPs(j).Gaps.Match = gap_vals(1);
				data.Hits(i).HSPs(j).Gaps.Possible = gap_vals(2);
				data.Hits(i).HSPs(j).Gaps.Percent = gap_vals(3);
			end
			
			if use_frame
				data.Hits(i).HSPs(j).Frame = sscanf(char(frame), 'Frame = %s %s %s');
			end
			
			if use_strand
				data.Hits(i).HSPs(j).Strand = sscanf(char(strand), 'Strand = %s %s %s');
			end
			
			hspAlign = '';
			seq_lines = size(alignments,2);
			
			for k = 1:seq_lines
				max_length = numel(char(alignments{k}(2)));
				lead_pos = 1 + max_length - numel(char(alignments{k}(1)));
				alignments{k}(2) = {alignments{k}{2}(lead_pos:max_length)};
				hspAlign = [hspAlign char(alignments{k})];
			end
			data.Hits(i).HSPs(j).Alignment = hspAlign;
			Q_start = query_indices{1}(1);
			Q_stop = query_indices{seq_lines}(2);
			data.Hits(i).HSPs(j).QueryIndices = str2double([Q_start Q_stop]);
			S_start = subj_indices{1}(1);
			S_stop = subj_indices{seq_lines}(2);
			data.Hits(i).HSPs(j).SubjectIndices = str2double([S_start S_stop]);
		end
	end
end

% Append statistical summary
data.Statistics = deblank(statstext);
