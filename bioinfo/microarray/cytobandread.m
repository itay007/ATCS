function S=cytobandread(filename)
%CYTOBANDREAD reads cytogenetic banding information.
%
%   S = CYTOBANDREAD(FILENAME) reads the G-banding data file FILENAME,
%   returning the data as a structure. The data file can be an NCBI
%   ideogram text file or a UCSC Genome Browser cytoband text file. The
%   return structure contains these fields:
%       ChromLabels
%       BandStartBPs
%       BandEndBPs
%       BandLabels
%       GieStains
% 
%   Note: For Homo sapiens, GieStains contains these values: gneg, gpos25,
%   gpos50, gpos75, gpos100, acen (centromere), stalk (repeat area, like
%   the p arm of Chromosome 22), and gvar.
% 
%   You can download cytogenetic banding information file from the NCBI or
%   the UCSC Genome Browser ftp site. For example, the current ideogram for
%   Homo sapiens can be downloaded from:
%   ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/mapview/ideogram.gz
%   or
%   ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBandIdeo.txt.gz
%
%   Examples:
%
%       % Read the ideogram of Homo sapiens.
%       hs_cytobands = cytobandread('hs_cytoBand.txt')
%
%   See also AFFYSNPCNVDEMO, BACACGHDEMO, CGHCBS, CHROMOSOMEPLOT.

%   Copyright 2007-2012 The MathWorks, Inc.


S =[];
bioinfochecknargin(nargin,1,mfilename);

if ~ischar(filename) || (iscellstr(filename)&& length(filename) > 1) 
    error(message('bioinfo:cytobandread:FileNameNotAString'));
end

try
    fid = fopen(filename, 'rt');
    % pass the empty lines
    pos = ftell(fid);
    while isempty(fgetl(fid))
        pos = ftell(fid);
    end
    fseek(fid, pos, -1);

    % Check for chr at the beginning of the first line
    tline = fgetl(fid);

    if ~isempty(regexp(tline, '^\#chromosome*', 'once'))% NCBI file with column header
        %-----------NCBI ideogram file -----------%
        c = textscan(fid, '%s%s%s%*d%*d%d%d%s%d%*d', 'Delimiter',' \b\t');
        fclose(fid);
        chromLabels = c{1};
        armNames = c{2};
        bandNames = c{3};       
        chromStarts = c{4};
        chromEnds = c{5};
        gieStainLabels = c{6};
        gieDensity = c{7};
        
        % Remove items without gie stain labels
		idx = cellfun('isempty', gieStainLabels);
        bandNames = bandNames(~idx);
        chromLabels = chromLabels(~idx);
        armNames = armNames(~idx);
        chromStarts = chromStarts(~idx);
        chromEnds = chromEnds(~idx);
        gieDensity = gieDensity(~idx);
        gieStainLabels = gieStainLabels(~idx);

        % Remove the extra arm labels. NCBI ideogram contain extra band
        % position information that is nor needed for ideogramplot.
        tmp = cellfun(@(x) isempty(regexp(x, '[a-z]', 'once')),...
                           bandNames, 'UniformOutput',false);
        idx = cell2mat(tmp);
        bandNames = bandNames(idx);
        chromLabels = chromLabels(idx);
        armNames = armNames(idx);
        chromStarts = chromStarts(idx);
        chromEnds = chromEnds(idx);
        gieDensity = gieDensity(idx);
        gieStainLabels = gieStainLabels(idx);

        % Concatenate arm and band labels
        bandLabels = strcat(armNames, bandNames);

        % Concatenate giestainlabel with density
        tmp = num2CellStr(gieDensity)';
        idx = strmatch('0', tmp);
        tmp(idx) = {''};
        gieStainLabels = strcat(gieStainLabels, tmp);
        
    elseif ~isempty(regexp(tline, '^chr.*', 'once')) % UCSC files with chr1 on line one
        %---- UCSC cytoband.txt file ---------------%
        fseek(fid, pos, -1);
        c = textscan(fid, '%s%d%d%s%s');
        fclose(fid);

        chromLabels = c{1};
        chromStarts = c{2};
        chromEnds = c{3};
        bandLabels = c{4};
        gieStainLabels = c{5};
        chromLabels = regexp(chromLabels, '(?<=^chr).*', 'match', 'once');
    else
        fclose(fid);
        error(message('bioinfo:cytobandread:UnknowCytoBandFileFormat', filename))
    end
    
    chromosome = convert2Chromosome(chromLabels);

    % Sort according to chromosome number - UCSC's order need to be
    % correct.
    [~, idx] = sort(chromosome);

    S.ChromLabels = chromLabels(idx);
    S.BandStartBPs = chromStarts(idx);
    S.BandEndBPs = chromEnds(idx);
    S.BandLabels = bandLabels(idx);
    S.GieStains = gieStainLabels(idx);
catch theException
    if strcmpi(theException.identifier,'bioinfo:cytobandread:UnknowCytoBandFileFormat')
		rethrow(theException)
	else
		error(message('bioinfo:cytobandread:CannotReadInput', filename))
    end
end
end % cytobandread function
%-------------------- Helper functions -----------------------------%
function chrs = convert2Chromosome(chromLabels)
% Convert chromosome label 'chr2' to a number 2

chrs = str2double(chromLabels);

chrs(strcmpi('X', chromLabels)) = max(chrs)+1;
chrs(strcmpi('Y', chromLabels)) = max(chrs)+1;
end % convert2Chromosome function
%-----------------------------------------------------------%
function c = num2CellStr(n)
% num2StrCell Convert vector of numbers to cell array of strings
N = numel(n);
c = cell(1,N);
for i=1:N
  c{i} = sprintf('%d', n(i));
end
end % num2CellStr function
